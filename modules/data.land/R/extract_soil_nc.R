#' Extract soil data from gssurgo
#'
#' Samples soil variabity within an area specified either as a radius around a
#' single lat/lon point OR as a spatial polygon.
#'
#' The samples are drawn from all soil types present in the area of interest,
#' as determined by the gSSURGO soil map at 30 meter pixel resolution.
#' Each soil component is sampled proportionally to its prevalence in the AOI,
#' with the restriction that each soil type is sampled at least once.
#' Note especially that this means sites with multiple soil types will always
#' have multiple output files, **even when** `size = 1`.
#'
#' @param outdir Output directory for writing down the netcdf file
#' @param lat Latitude of center point (single numeric value)
#' @param lon Longitude of center point (single numeric value)
#' @param aoi Spatial polygon defining area of interest.
#'  When this is specified, all three of `lat`, `lon`, and `radius` are ignored
#' @param size Ensemble size
#' @param radius distance in meters to sample around center point
#' @param depths depths in m of boundaries between soil layers
#'
#' @return list of paths to the generated soil netcdf files
#'
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @examples
#' \dontrun{
#'    outdir  <- "~/paleon/envTest"
#'    lat     <- 40
#'    lon     <- -80
#'    PEcAn.data.land::extract_soil_gssurgo(outdir, lat, lon)
#' }
#' @author Hamze Dokoohaki, Akash
#' @export
#'
extract_soil_gssurgo <- function(outdir,
                                 lat, lon, aoi = NULL,
                                 size = 1, radius = 500,
                                 depths = c(0.15, 0.30, 0.60)) {
  if (missing(aoi)) {
    aoi <- data.frame(lon = lon, lat = lat) |>
      terra::vect(crs = "epsg:4326") |>
      terra::buffer(radius)
  }

  all.soil.ens <- list()
  PEcAn.logger::logger.info("Querying gSSURGO Web Coverage Service for map unit keys")
  mu_raster <- soilDB::mukey.wcs(aoi = aoi, db = "gSSURGO", res = 30)

  # Extract unique mukeys and their pixel counts for area weighting
  mu_weights <- mu_raster %>%
    terra::values(dataframe = TRUE, na.rm = TRUE) %>%
    table() %>%
    as.data.frame() %>%
    dplyr::mutate(
      mukey = as.character(mukey),
      mu_wt = .data$Freq / sum(.data$Freq)
    )
  if (nrow(mu_weights) == 0) {
    PEcAn.logger::logger.severe("No mapunit keys were found for this site.")
  }

  # look up component soils (aka soil types) found in these map units
  components <-
    soilDB::get_component_from_SDA(
      paste("mukey IN", soilDB::format_SQL_in_statement(mu_weights$mukey)),
      duplicates = TRUE,
      childs = FALSE
    ) %>%
    dplyr::select("mukey", "cokey", "comppct_r") %>%
    dplyr::mutate(dplyr::across(c("mukey", "cokey"), as.character))
  component_weights <- mu_weights %>%
    dplyr::left_join(components, by = "mukey") %>%
    dplyr::mutate(comp_wt = (.data$comppct_r / 100) * .data$mu_wt)
  if (anyNA(component_weights)) {
    # This will only(?) happen when no components found for a map unit
    ok <- stats::complete.cases(component_weights)
    PEcAn.logger::logger.warn(
      "Dropping map units with no data found:",
      sQuote(component_weights$mukey[!ok])
    )
    component_weights <- component_weights[ok, ]
  }
  component_weights <- component_weights %>%
    dplyr::group_by(.data$cokey) %>%
    # consolidate any components that appear in multiple map units
    dplyr::summarize(
      comp_wt = sum(.data$comp_wt, na.rm = TRUE),
      .groups = "drop"
    )

  # Get horizon-by-horizon properties of each component
  hz_data <- tryCatch({
    soilDB::get_chorizon_from_SDA(
      paste("c.cokey IN",
            soilDB::format_SQL_in_statement(component_weights$cokey)),
      duplicates = FALSE,
      childs = TRUE,
      nullFragsAreZero = TRUE,
      droplevels = TRUE
    )
  }, error = function(e) {
    PEcAn.logger::logger.error("Failed to fetch soil horizon data:", e$message)
    return(NULL)
  })

  # Calculate properties of each requested layer
  # as the weighted average of all horizons appearing in it
  depths_cm <- depths * 100
  all_soil_data <- list()
  for (i in seq_along(depths_cm)) {
    top_depth <- if (i == 1) 0 else depths_cm[i - 1]
    bottom_depth <- depths_cm[i]

    depth_data <- hz_data %>%
      dplyr::filter(.data$hzdept_r < bottom_depth
                    & .data$hzdepb_r > top_depth) %>%
      dplyr::mutate(
        hz_top_adj = pmax(.data$hzdept_r, top_depth),
        hz_bot_adj = pmin(.data$hzdepb_r, bottom_depth),
        hz_thickness = .data$hz_bot_adj - .data$hz_top_adj
      ) %>%
      dplyr::group_by(.data$cokey) %>%
      dplyr::summarize(
        depth_layer = depths[i],
        hzdept_r = top_depth,
        hzdepb_r = bottom_depth,
        dplyr::across(
          c("sandtotal_r", "silttotal_r", "claytotal_r",
            "om_r", "dbthirdbar_r", "total_frags_pct"),
          \(x) stats::weighted.mean(x, .data$hz_thickness, na.rm = TRUE)
        ),
        .groups = "drop"
      )
    all_soil_data[[i]] <- depth_data
  }

  # Transform to match original code format
  soilprop <- do.call(rbind, all_soil_data)
  soilprop.new <- soilprop %>%
    dplyr::select(
      fraction_of_sand_in_soil = "sandtotal_r",
      fraction_of_silt_in_soil = "silttotal_r",
      fraction_of_clay_in_soil = "claytotal_r",
      depth_layer = "depth_layer",
      soil_depth = "hzdept_r",
      soil_depth_bottom = "hzdepb_r",
      organic_matter_pct = "om_r",
      bulk_density = "dbthirdbar_r",
      coarse_fragment_pct = "total_frags_pct",
      cokey = "cokey"
    ) %>%
    dplyr::mutate(
      dplyr::across(
        c(dplyr::starts_with("fraction_of"), "coarse_fragment_pct"),
        ~ . / 100
      ),
      horizon_thickness_cm = .data$soil_depth_bottom - .data$soil_depth,
      soil_organic_carbon_stock = PEcAn.data.land::soc2ocs(
        soc_percent = PEcAn.data.land::om2soc(.data$organic_matter_pct),
        bulk_density = .data$bulk_density,
        thickness = .data$horizon_thickness_cm,
        coarse_fraction = .data$coarse_fragment_pct
      )
    ) %>%
    dplyr::filter(stats::complete.cases(.))

  if (nrow(soilprop.new) == 0) {
    PEcAn.logger::logger.error("No valid soil properties after filtering")
    return(NULL)
  }

  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  soil.data.gssurgo <- list(
    fraction_of_sand_in_soil = soilprop.new$fraction_of_sand_in_soil,
    fraction_of_silt_in_soil = soilprop.new$fraction_of_silt_in_soil,
    fraction_of_clay_in_soil = soilprop.new$fraction_of_clay_in_soil,
    soil_depth = soilprop.new$soil_depth,
    soil_organic_carbon_stock = soilprop.new$soil_organic_carbon_stock
  )

  all.soil.ens <- c(all.soil.ens, list(soil.data.gssurgo))

  # Generate modeled ensembles
  tryCatch({
    # Dirichlet modeling per component
    simulated.soil.props <- soilprop.new %>%
      split(list(.$cokey, .$depth_layer)) %>%
      purrr::map_df(function(component_group) {
        tryCatch({
          texture_data <- component_group[, c("fraction_of_sand_in_soil",
                                              "fraction_of_silt_in_soil",
                                              "fraction_of_clay_in_soil")] %>%
            as.matrix()

          if (nrow(texture_data) == 0) return(NULL)

          dir.model <- sirt::dirichlet.mle(texture_data)
          alpha <- matrix(dir.model$alpha, nrow = size, ncol = length(dir.model$alpha), byrow = TRUE)
          simulated.soil <- sirt::dirichlet.simul(alpha)

          # SOC modeling
          soc_mean <- component_group$soil_organic_carbon_stock
          soc_sd <- stats::sd(soc_mean, na.rm = TRUE)
          n_depths <- length(soc_mean)

          if (n_depths == 1 || is.na(soc_sd) || soc_sd == 0) {
            simulated_soc <- rep(NA_real_, size)
          } else {
            shape <- (mean(soc_mean, na.rm = TRUE)^2) / (soc_sd^2)
            rate <- mean(soc_mean, na.rm = TRUE) / (soc_sd^2)
            simulated_soc <- stats::rgamma(size, shape = shape, rate = rate)
          }

          result_df <- data.frame(
            fraction_of_sand_in_soil = simulated.soil[, 1],
            fraction_of_silt_in_soil = simulated.soil[, 2],
            fraction_of_clay_in_soil = simulated.soil[, 3],
            soil_depth = component_group$soil_depth,
            cokey = unique(component_group$cokey),
            soil_organic_carbon_stock = simulated_soc
          )

          return(result_df)
        },
        error = function(e) {
          PEcAn.logger::logger.warn(conditionMessage(e))
          return(NULL)
        })
      })

    # generate weighted profiles
    soil.profiles <- simulated.soil.props %>%
      split(.$cokey) %>%
      purrr::map(function(soiltype.sim) {
        wt <- component_weights %>%
          dplyr::filter(.data$cokey == unique(soiltype.sim$cokey))
        sizein <- wt$comp_wt * size
        1:ceiling(sizein) %>%
          purrr::map(function(x) {
            soiltype.sim %>%
              split(.$soil_depth) %>%
              purrr::map_dfr(~.x[x, ])
          })
      }) %>%
      purrr::flatten()

    # convert profiles to ensemble arrays
    all.soil.ens <- soil.profiles %>%
      purrr::map(function(SEns) {
        SEns <- SEns[, names(SEns) != "cokey"]
        names(SEns) %>%
          purrr::map(function(var) {
            as.numeric(unlist(SEns[, var]))
          }) %>%
          stats::setNames(names(SEns))
      }) %>%
      c(all.soil.ens, .)
  },
  error = function(e) {
    PEcAn.logger::logger.warn(conditionMessage(e))
  })

  # generate NetCDF files
  out.ense <- (1:length(all.soil.ens)) %>%
    purrr::map(function(i) {
      tryCatch({
        prefix <- paste0("gSSURGO_soil_", i)
        new.file <- file.path(outdir, paste0(prefix, ".nc"))
        suppressWarnings({
          PEcAn.data.land::soil2netcdf(all.soil.ens[[i]], new.file)
        })
        new.file
      },
      error = function(e) {
        PEcAn.logger::logger.warn(conditionMessage(e))
        return(NULL)
      })
    }) %>%
    # remove nulls
    purrr::discard(is.null) %>%
    stats::setNames(rep("path", length(.)))
  return(out.ense)
}


#' Extract soil data from the gridpoint closest to a location
#'
#' @param in.file path to netcdf file containing soil data
#' @param outdir directory in which to write netcdf file of extracted data.
#'  Output filename will be the same as input filename.
#' @param lat,lon location in decimal degrees.
#'  Data will be extracted from the point in `in.file` that is nearest this
#'
#' @return path to netCDF file containing extracted data
#' @export
#'
#' @examples
#' \dontrun{
#' in.file <- "~/paleon/env_paleon/soil/paleon_soil.nc"
#' outdir  <- "~/paleon/envTest"
#' lat     <- 40
#' lon     <- -80
#' PEcAn.data.land::extract_soil_nc(in.file,outdir,lat,lon)
#' }
extract_soil_nc <- function(in.file,outdir,lat,lon){
  
  ## open soils
  nc <- ncdf4::nc_open(in.file)
  
  ## extract lat/lon
  dims <- names(nc$dim)
  lat.dim <- dims[grep("^lat",dims)]
  lon.dim <- dims[grep("^lon",dims)]
  soil.lat <- ncdf4::ncvar_get(nc, lat.dim)
  soil.lon <- ncdf4::ncvar_get(nc, lon.dim)
  
  ## check in range
  dlat <- abs(stats::median(diff(soil.lat)))
  dlon <- abs(stats::median(diff(soil.lon)))
  if(lat < (min(soil.lat)-dlat) | lat > (max(soil.lat)+dlat)){
    PEcAn.logger::logger.error("site lat out of bounds",lat,range(soil.lat))
  }
  if(lon < (min(soil.lon)-dlon) | lon > (max(soil.lon)+dlon)){
    PEcAn.logger::logger.error("site lon out of bounds",lon,range(soil.lon))
  }
  if(dims[1] == lat.dim){
    soil.row <- which.min(abs(lat-soil.lat))
    soil.col <- which.min(abs(lon-soil.lon))
  } else if(dims[1] == lon.dim){
    soil.col <- which.min(abs(lat-soil.lat))
    soil.row <- which.min(abs(lon-soil.lon))
  } else {
    PEcAn.logger::logger.error("could not determine lat/lon dimension order:: ",dims)
  }
  
  ## extract raw soil data
  soil.data <- list()
  soil.vars <- names(nc$var)
  for(i in seq_along(soil.vars)){
    if(length(dims) == 2){
      soil.data[[soil.vars[i]]] <- ncdf4::ncvar_get(nc,soil.vars[i])[soil.row,soil.col]
    } else {
      ## assuming there's a 3rd dim of soil depth profile
      soil.data[[soil.vars[i]]] <- ncdf4::ncvar_get(nc,soil.vars[i])[soil.row,soil.col,]
    }
  }
  ncdf4::nc_close(nc)
  
  ## PalEON / MSTMIP / UNASM hack
  # t_ variables are topsoil layer (0– 30 cm) and
  # s_ variables are subsoil layer (30–100 cm)
  depth <- ncdf4::ncdim_def(name = "depth", units = "meters", vals = c(0.3,1.0), create_dimvar = TRUE)  
  dvars <- soil.vars[grep("t_",soil.vars,fixed=TRUE)]
  for(i in seq_along(dvars)){
    svar <- sub("t_","s_",dvars[i])
    soil.data[[dvars[i]]] <- c(soil.data[[dvars[i]]],soil.data[[svar]]) ## combine different depths
    soil.data[[svar]] <- NULL  ## drop old variable
    names(soil.data)[which(names(soil.data) == dvars[i])] <- sub("t_","",dvars[i]) ## rename original
  }
  
  
  ## name/unit conversions 
  soil.data$sand   <- soil.data$sand/100
  soil.data$silt   <- soil.data$silt/100
  soil.data$clay   <- soil.data$clay/100
  soil.data$oc     <- soil.data$oc/100
  soil.data$gravel <- soil.data$gravel/100
  soil.data$ref_bulk <- PEcAn.utils::ud_convert(soil.data$ref_bulk,"g cm-3","kg m-3")
  names(soil.data)[which(names(soil.data) == "clay")] <- "fraction_of_clay_in_soil"
  names(soil.data)[which(names(soil.data) == "sand")] <- "fraction_of_sand_in_soil"
  names(soil.data)[which(names(soil.data) == "silt")] <- "fraction_of_silt_in_soil"
  names(soil.data)[which(names(soil.data) == "gravel")] <- "fraction_of_gravel_in_soil"
  names(soil.data)[which(names(soil.data) == "ref_bulk")] <- "soil_bulk_density"
  names(soil.data)[which(names(soil.data) == "ph")]   <- "soil_ph"
  names(soil.data)[which(names(soil.data) == "cec")]  <- "soil_cec" ## units = meq/100g
  names(soil.data)[which(names(soil.data) == "oc")]   <- "soilC"  ## this is currently the BETY name, would like to change and make units SI
  
  ## calc new filename
  prefix <- tools::file_path_sans_ext(basename(in.file))
  new.file <- file.path(outdir,paste0(prefix,".nc"))
  
  ## Calculate soil parameters and export to netcdf
  PEcAn.data.land::soil2netcdf(soil.data,new.file)
  
  return(new.file)
  
}


#' Get standard units for a soil variable
#'
#' Given SSURGO names for soil properties, looks up their standard units.
#' Note that names must match exactly.
#'
#' Supported variables are:
#'  * `soil_depth`
#'  * `soil_cec`
#'  * `fraction_of_clay_in_soil`
#'  * `fraction_of_sand_in_soil`
#'  * `fraction_of_silt_in_soil`
#'  * `fraction_of_gravel_in_soil`
#'  * `volume_fraction_of_water_in_soil_at_saturation`
#'  * `volume_fraction_of_water_in_soil_at_field_capacity`
#'  * `volume_fraction_of_condensed_water_in_dry_soil`
#'  * `volume_fraction_of_condensed_water_in_soil_at_wilting_point`
#'  * `soilC`
#'  * `soil_ph`
#'  * `soil_bulk_density`
#'  * `soil_type`
#'  * `soil_hydraulic_b`
#'  * `soil_water_potential_at_saturation`
#'  * `soil_hydraulic_conductivity_at_saturation`
#'  * `thcond0`
#'  * `thcond1`
#'  * `thcond2`
#'  * `thcond3`
#'  * `soil_thermal_conductivity`
#'  * `soil_thermal_conductivity_at_saturation`
#'  * `soil_thermal_capacity`
#'  * `soil_albedo`
#'  * `slpotwp`
#'  * `slpotcp`
#'  * `slcpd`
#'  * `slden`
#'  * `soil_organic_carbon_stock`
#'
#' @param varname character vector. See details
#'
#' @return character matrix with columns `var` and `unit`
#' @md
#' @export
#'
#' @examples
#' soil.units("soil_albedo")
soil.units <- function(varname = NA){
  variables <- as.data.frame(matrix(c("soil_depth","m",
                                      "soil_cec","meq/100g",
                                      "fraction_of_clay_in_soil","1",
                                      "fraction_of_sand_in_soil","1",
                                      "fraction_of_silt_in_soil","1",
                                      "fraction_of_gravel_in_soil","1",
                                      "volume_fraction_of_water_in_soil_at_saturation","m3 m-3",
                                      "volume_fraction_of_water_in_soil_at_field_capacity","m3 m-3",
                                      "volume_fraction_of_condensed_water_in_dry_soil","m3 m-3",
                                      "volume_fraction_of_condensed_water_in_soil_at_wilting_point","m3 m-3",
                                      "soilC","percent",
                                      "soil_ph","1",
                                      "soil_bulk_density","kg m-3",
                                      "soil_type","string",
                                      "soil_hydraulic_b","1",
                                      "soil_water_potential_at_saturation","m",
                                      "soil_hydraulic_conductivity_at_saturation","m s-1",
                                      "thcond0","W m-1 K-1",
                                      "thcond1","W m-1 K-1",
                                      "thcond2","1",
                                      "thcond3","1",
                                      "soil_thermal_conductivity","W m-1 K-1", 
                                      "soil_thermal_conductivity_at_saturation","W m-1 K-1", 
                                      "soil_thermal_capacity","J kg-1 K-1",
                                      "soil_albedo","1",
                                      "slpotwp","m",
                                      "slpotcp","m",
                                      "slcpd","J m-3 K-1",
                                      "slden","kg m-3",
                                      "soil_organic_carbon_stock","kg m-2"
    ),
  ncol=2,byrow = TRUE))
  colnames(variables) <- c('var','unit')
  
  unit = which(variables$var == varname)
  
  if(length(unit) == 0){
    if(is.na(varname)){
      return(variables)
    } else {
      return(NA)
    }
  }else{
    unit = as.character(variables$unit[unit])
    return(unit)
  }
  
}
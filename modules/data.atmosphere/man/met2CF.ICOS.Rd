% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/met2CF.ICOS.R
\name{met2CF.ICOS}
\alias{met2CF.ICOS}
\title{Convert variables ICOS variables to CF format.}
\usage{
met2CF.ICOS(
  in.path,
  in.prefix,
  outfolder,
  start_date,
  end_date,
  format,
  overwrite = FALSE,
  ...
)
}
\arguments{
\item{in.path}{path to the input ICOS product CSV file}

\item{in.prefix}{name of the input file}

\item{outfolder}{path to the directory where the output file is stored. If specified directory does not exists, it is created.}

\item{start_date}{start date of the input file}

\item{end_date}{end date of the input file}

\item{format}{format is data frame or list with elements as described below
  REQUIRED:
  format$header = number of lines of header
  format$vars is a data.frame with lists of information for each variable to read, at least airT is required
    format$vars$input_name = Name in CSV file
    format$vars$input_units = Units in CSV file
    format$vars$bety_name = Name in BETY
  OPTIONAL:
  format$lat = latitude of site
  format$lon = longitude of site
  format$na.strings = list of missing values to convert to NA, such as -9999
  format$skip = lines to skip excluding header
  format$vars$column_number = Column number in CSV file (optional, will use header name first)
Columns with NA for bety variable name are dropped.}

\item{overwrite}{overwrite should existing files be overwritten. Default False.}

\item{...}{used when extra arguments are present.}
}
\value{
information about the output file
}
\description{
Variables present in the output netCDF file:
air_temperature, air_temperature, relative_humidity,
specific_humidity, water_vapor_saturation_deficit,
surface_downwelling_longwave_flux_in_air,
surface_downwelling_shortwave_flux_in_air,
surface_downwelling_photosynthetic_photon_flux_in_air, precipitation_flux,
eastward_wind, northward_wind
}

############ Retrives soil data from gssurgo
#' This function queries the gSSURGO database for a series of map unit keys
#'
#' @param mukeys map unit key from gssurgo
#' @param fields a character vector of the fields to be extracted. See details and the default argument to find out how to define fields.
#'
#' @return a dataframe with soil properties. units can be looked up from database documentation
#' @export
#'
#' @details 
#' Full documention of available tables and their relationships can be found here \url{www.sdmdataaccess.nrcs.usda.gov/QueryHelp.aspx}
#' There have been occasions where NRCS made some minor changes to the structure of the API which this code is where those changes need
#' to be implemneted here.
#' Fields need to be defined with their associate tables. For example, sandtotal is a field in chorizon table which needs to be defined as chorizon.sandotal_(r/l/h), where 
#' r stands for the representative value, l stand for low and h stands for high. At the momeent fields from mapunit, component, muaggatt, and chorizon tables can be extracted. 
#'
#' @examples
#' \dontrun{
#'  PEcAn.data.land::gSSURGO.Query(
#'    fields = c(
#'      "chorizon.cec7_r", "chorizon.sandtotal_r",
#'      "chorizon.silttotal_r","chorizon.claytotal_r",
#'      "chorizon.om_r","chorizon.hzdept_r","chorizon.frag3to10_r",
#'      "chorizon.dbovendry_r","chorizon.ph1to1h2o_r",
#'      "chorizon.cokey","chorizon.chkey"))
#' }
gSSURGO.Query <- function(mukeys=2747727,
                        fields=c("chorizon.sandtotal_r",
                                 "chorizon.silttotal_r",
                                 "chorizon.claytotal_r")){
  #browser()
  #               ,
  ######### Reteiv soil
  headerFields <-
    c(Accept = "text/xml",
      Accept = "multipart/*",
      'Content-Type' = "text/xml; charset=utf-8",
      SOAPAction = "http://SDMDataAccess.nrcs.usda.gov/Tabular/SDMTabularService.asmx/RunQuery")
  
  body <- paste('<?xml version="1.0" encoding="utf-8"?>
               <soap:Envelope xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
               <soap:Body>
               <RunQuery xmlns="http://SDMDataAccess.nrcs.usda.gov/Tabular/SDMTabularService.asmx">
               <Query>
               SELECT mapunit.mukey, component.cokey, component.mukey, component.comppct_r, ',paste(fields, collapse = ", "),',
               muaggatt.aws050wta from mapunit
               join muaggatt on mapunit.mukey=muaggatt.mukey
               join component on mapunit.mukey=component.mukey
               join chorizon on component.cokey=chorizon.cokey
               where mapunit.mukey in (', paste(mukeys,collapse = ", "),');
               </Query>
               </RunQuery>
               </soap:Body>
               </soap:Envelope>')
  reader <- RCurl::basicTextGatherer()
  out <- RCurl::curlPerform(url = "https://SDMDataAccess.nrcs.usda.gov/Tabular/SDMTabularService.asmx",
                          httpheader = headerFields,  postfields = body,
                          writefunction = reader$update
  )
  suppressWarnings(
    suppressMessages({
      xml_doc <- XML::xmlTreeParse(reader$value())
      xmltop  <- XML::xmlRoot(xml_doc)
      tablesxml <- (xmltop[[1]]["RunQueryResponse"][[1]]["RunQueryResult"][[1]]["diffgram"][[1]]["NewDataSet"][[1]])
    })
  )
  
  #parsing the table  
  tryCatch({
    suppressMessages(
      suppressWarnings({
        tables <- XML::getNodeSet(tablesxml,"//Table")
        
        ##### All datatables below newdataset
        # This method leaves out the variables are all NAs  - so we can't have a fixed naming scheme for this df
        dfs <- tables %>%
          purrr::map_dfr(function(child){
            #converting the xml obj to list
            allfields <- XML::xmlToList(child)
            remov <- names(allfields) %in% c(".attrs")
            #browser()
            names(allfields)[!remov] %>%
              purrr::map_dfc(function(nfield){
                #browser()
                outv <- allfields[[nfield]] %>% unlist() %>% as.numeric
                ifelse(length(outv) > 0, outv, NA)
              })%>%
              as.data.frame() %>%
              `colnames<-`(names(allfields)[!remov])
          })%>%
          dplyr::select(comppct_r:mukey) %>%
          dplyr::select(-aws050wta)
        
      })
    )
    
    
    return(dfs)
  },
  error=function(cond) {
    print(cond)
    return(NULL)
  })
  
}



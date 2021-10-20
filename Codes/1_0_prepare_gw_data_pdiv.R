#' ---
#' title: "Base Data Construction"
#' author: "Taro Mieno"
#' output:
#'   html_document:
#'     number_sections: yes
#'     theme: flatly
#'     highlight: zenburn
#'     toc_float: yes
#'     toc: yes
#'     toc_depth: 3
#' geometry: margin=1in
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(
  echo = TRUE,
  cache = FALSE,
  comment = NA,
  message = FALSE,
  warning = FALSE,
  tidy = FALSE,
  cache.lazy = FALSE,
  #--- figure ---#
  dpi = 400,
  fig.width = 7.5,
  fig.height = 5,
  out.width = "750px",
  out.height = "500px"
)

#+ library, include = FALSE, cache = FALSE
library(fields)
library(tidyverse)
library(data.table)
library(sf)
library(here)

# /*=================================================*/
#' # Objectives
# /*=================================================*/
#' + combine 1991-2014 data and 2015-2016 data together
#' + merge with weather and water rights data

# /*=================================================*/
#' # Combine 1991-2014 data and 2015-2016 data together
# /*=================================================*/

# /*----------------------------------*/
#' ## 1991-2014 (GW-data)
# /*----------------------------------*/
data_91_14 <-
  here("Shared/Data/Kansas/PDIV_YearLevel.csv") %>% 
  fread() %>%
  setnames(names(.), tolower(names(.))) %>%
  # get rid of surface water use observations
  .[source_g == 1, ] 

#=== groundwater use data ===#
gw_data_91_14 <- 
  data_91_14[, .(wr_id, pdiv_id, wuapers_id, wua_year, acres_irr, af_used)] 

# data_91_14[is.na(predev_satthick), ]

#/*----------------------------------*/
#' ## 1991-2014 (time-invariant data)
#/*----------------------------------*/
#' Some of the time-invariant variables are extracted
#' and merged back later because newer datasets do not hold
#' the information 

#=== time-invariant variables to keep ===#
ti_vars <- c(
  "fpdiv_key", "wrf_key", "wr_num", "trs", "gwmd_num", "priority_days", "latitude", "longitude",
  "awc", "bulkdensity", "claytotal", "ec", "ksat", "kffact", "kwfact", "lep", "om", "sandtotal", "slope", "silttotal"
)

#=== get all the time-invariant variables ===#
ti_data <- data_91_14[, c("pdiv_id", ti_vars), with = FALSE] %>%
  unique(by = "pdiv_id")

#/*----------------------------------*/
#' ## 2015-2016 (GW-data)
#/*----------------------------------*/
gw_data_15_16 <- 
  here("Shared/Data/Kansas/water_use2015-2016.csv") %>% 
  fread() %>% 
  .[, .(wr_id, pdiv_id, wuapers_id, wua_year, acres_irr, af_used)]


#/*----------------------------------*/
#' ## 2017-2019 (GW-data)
#/*----------------------------------*/
gw_data_17_19 <- 
  here("Shared/Data/Kansas/water_use_data_2017-2019.dta") %>% 
  haven::read_dta() %>% 
  data.table() %>% 
  .[, .(wr_id, pdiv_id, wuapers_id, wua_year, acres_irr, af_used)]

#/*----------------------------------*/
#' ## Combine the three datasets
#/*----------------------------------*/
data_combined <- 
  rbindlist(
    list(gw_data_91_14, gw_data_15_16, gw_data_17_19), 
    fill = TRUE
  ) 

#/*=================================================*/
#' # Merge with the time-invariant data and save
#/*=================================================*/
#/*----------------------------------*/
#' ## Merge with the time-invariant data
#/*----------------------------------*/
data_ti <- ti_data[data_combined, on = "pdiv_id"] %>%
  setnames(c("pdiv_id", "wua_year"), c("site", "year"))

#/*----------------------------------*/
#' ## Merge with water rights data
#/*----------------------------------*/
#--- read the water rights data ---#
wr_data <- 
  here("Shared/Data/RawData/wr_id_authorized_quantities.csv") %>% 
  fread() %>%
  setnames(names(.),tolower(names(.)))

#--- merge with water rights data ---#
data <- wr_data[data_ti, on = 'wr_id'] %>% 
  .[af_used != 0, ] %>% 
  .[!is.na(af_used), ] %>% 
  .[!is.na(auth_irr), ]

#/*=================================================*/
#' # Save
#/*=================================================*/
saveRDS(data, here("Shared/Data/ProcessedData/gw_ks_2019.rds"))



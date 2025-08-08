library(tidyverse)
library(skimr)
library(leaflet)
library(sf)
library(magrittr)

setwd("C:/Users/hg000051/Downloads/")

# Function which reads tsv
  readtsv <- function(filepath){
    
    read.table(filepath, sep = '\t', header =T)
  }

# Load in data files from CANARY. Are these all the files I need, not sure
  Tgt_Screening_GCMS <- readtsv("GCMS Target Screening detections raw data 2025-08-08.tsv")
  EA_NLS_Target_GCMS <- readtsv("GCMS EA NLS Target Database 2025-08-08.tsv")
  GCMS_Open_Tgt <- readtsv("GCMS Target and Non-Targeted Screening for Open Data (detes 3106 7299 4994) 2025-08-08.tsv")
  
  LCMS_Open_Tgt <- readtsv("LCMS Target and Non-Targeted Screening for Open Data (detes 3106 7299 4994) 2025-08-08.tsv")  
  Tgt_Screening_LCMS <- readtsv("LCMS Target Screening detections raw data 2025-08-08.tsv")
  EA_NLS_Target_LCMS <- readtsv("LCMS EA NLS Target Database 2025-08-08.tsv")

#RSN sites
  RSN_Sites <- readxl::read_excel("RSN_Site_Locations_PBI_Download.xlsx")

# Explore
  GCMS_Open_Tgt %>% skim()
  
  LCMS_Open_Tgt %>% skim()

# Transform
  names(LCMS_Open_Tgt)==names(GCMS_Open_Tgt)
  
  Open_Tgt <- rbind(LCMS_Open_Tgt, GCMS_Open_Tgt)
  
  Open_Tgt_sf <- Open_Tgt %>% st_as_sf(coords = c("Longitude","Latitude"), crs=4326)
  RSN_Sites_sf <- RSN_Sites %>% st_as_sf(coords = c("GRTS longitude (WGS84)","GRTS latitude (WGS84)"), crs=4326)
  
# Confirm if site ids match
  
  Open_Tgt_sf %>% 
    filter(str_detect(Sample_Site_ID, "RSN"))
  
# Visualise
  
  Open_Tgt %>% 
    filter(str_detect(Sample_Site_ID, "RSN")) %>% 
    leaflet() %>% 
    addProviderTiles(providers$Esri) %>% 
    addCircleMarkers(radius = 3, color = "seagreen2")
  
  
  
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

# LC/GC-MS library
  LCGC_Lib <- readxl::read_excel("Combined GCMS  LCMS Database MAY 2025.xlsx", skip=2) %>% 
      rename("Compound_Name" ="...1" ,
               "CAS#" = "...2",
             "Description" = "...3",
             "GCMS_LOD" = "GC-MS",
             "LCMS_LOD" = "LC-MS",
             "Added_LC" = "...6",
             "Added_GC" = "...7")
  
  
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
  
 RSN_MS <-  Open_Tgt_sf %>% 
    filter(str_detect(Sample_Site_ID, "RSN"))

# Identify which library compounds we didn't find on the RSN network
 # Do we need to screen for these?
  
  Comps_exc <-  LCGC_Lib %>% 
          filter(!Compound_Name %in%  RSN_MS$Compound_Name)
  
  Comps_exc %>% distinct()  # list of compounds on the library but not found in the environment.
  
# Identify how many compounds are found at RSN sites by LCMS methods vs GCMS.
  
  RSN_MS %>% 
    group_by(Sample_Site_ID) %>% 
    distinct(Sample_Site_ID, Compound_Name, method) %>% 
    mutate(Detected = TRUE) %>% 
  pivot_wider(names_from = method, values_from = Detected, values_fill = FALSE)
  
# Identify if there are any national trends
  
  
 
  
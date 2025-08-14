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
# Extract month and year for temporal variance.

  RSN_MS <-  Open_Tgt_sf %>% 
    filter(str_detect(Sample_Site_ID, "RSN"))
  
 # Observe concentration normality- large right tail
  ggplot(RSN_MS, aes(x = Concentration)) +
         geom_histogram(aes(y = ..density..), bins = 50) +
         geom_density(color = "seagreen2", size = 1)
  
  ggplot(RSN_MS, aes(x = logConc)) +
    geom_histogram(aes(y = ..density..), bins = 50) +
    geom_density(color = "seagreen2", size = 1)
  

 # Log transform concentration and capture month to test temporal component
 
  RSN_MS %<>% 
    mutate(
      month = as.factor(month(Sample_datetime)),
      season = case_when(
               month %in% c(3,4,5) ~ "Spring",
               month %in% c(6,7,8) ~ "Summer",
               month %in% c(9,10,11) ~ "Autumn",
               month %in% c(12,1,2) ~ "Winter"
      ),
      logConc = log(Concentration)
    )
 
 # Test to see if temporal variability is statistically significant and see how it compares to spatial variability 

 
  lme1 <- lmer(logConc ~ 1 + (1 | Sample_Site_ID), data = RSN_MS)
  lme1
 
  # Add seasonality as an effect instead of just spatial. 
 lme2 <- lme4::lmer(logConc ~ 1 + (1| Sample_Site_ID) + (1 + season), data = RSN_MS)
 lme2
 # spring and winter have notably lower concs -0.30 ish whilst summer is 0.03 different to autumn. 
 
 anova(lme1,lme2)
 # seasonality effect is statistically significant. AIC difference >10 so lme1 better fit.
 
 # Try seasonality as random effect 
 
 # Add seasonality as a random effect instead of just spatial. 
 lme2r <- lme4::lmer(logConc ~ 1 + (1| Sample_Site_ID) + (1 | season), data = RSN_MS)
 lme2r
 
 # residuals explains more than spatial or temporal - what we can't account for. 
 # Spaital explains ~double the amount of variance than temporal. 
 
 lme3 <- lme4::lmer(logConc ~ 1 + (1| Sample_Site_ID) + (1 | month), data = RSN_MS)
 lme3
 
 # Resid 1.693
 # concs rise in summer/ autumn months.

 anova(lme2,lme3)
 # 
 
 # Should probably group by determinand group, like pesticide, PFAS etc.
 
 
 
 
 # ICC after 
 
 vc <- as.data.frame(VarCorr(lme2r))
 total_var <- sum(vc$vcov)
  
    icc_site  <- vc$vcov[vc$grp == "Sample_Site_ID"] / total_var
 icc_month <- vc$vcov[vc$grp == "month"] / total_var
  icc_resid <- vc$vcov[vc$grp == "Residual"] / total_var
  
    icc_table <- data.frame(
          Component = c("Spatial (Site)", "Temporal (Month)", "Residual"),
          Proportion = c(icc_site, icc_month, icc_resid)
      )

 
 #9% variance is spatial
# 2% is temporal 
# 90% is within site noise
 
 # spatial variance > temporal variance. More gain from increasing number of sites than frequency at each site. 
 # Tes
 
 
    vc <- as.data.frame(VarCorr(lme2r))
    total_var <- sum(vc$vcov)
    
    # Calculate proportion of variance (ICC)
    icc_site   <- vc$vcov[vc$grp == "Sample_Site_ID"] / total_var
    icc_season <- vc$vcov[vc$grp == "season"] / total_var
    icc_resid  <- vc$vcov[vc$grp == "Residual"] / total_var
    
    # Create table for reporting
    icc_table <- data.frame(
      Component  = c("Spatial (Site)", "Temporal (Season)", "Residual"),
      Proportion = c(icc_site, icc_season, icc_resid)
    )
    
    # View ICC table
    icc_table
  # More site variance explained than seasonal still, checks out.
  # temp variation is halved when using seasons to compare temporal variability. Meaning variability between months inside a season gets lost. 
 # seasonal monitoring would miss within month variability. 
 
# Identify which library compounds we didn't find on the RSN network
 # Do we need to screen for these?
  
  Comps_exc <-  LCGC_Lib %>% 
          filter(Compound_Name %in%  RSN_MS$Compound_Name)
  
  Comps_exc %>% distinct()  # list of compounds on the library but not found in the environment.
  
# Identify how many compounds are found at RSN sites by LCMS methods vs GCMS.
  
  RSN_MS %>% 
    group_by(Sample_Site_ID) %>% 
    distinct(Sample_Site_ID, Compound_Name, method) %>% 
    mutate(Detected = TRUE) %>% 
  pivot_wider(names_from = method, values_from = Detected, values_fill = FALSE)
  
# 1 

  # Test for temporal vs spatial variability. 
  

  
 
  
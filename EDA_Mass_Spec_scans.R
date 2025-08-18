library(tidyverse)
library(skimr)
library(leaflet)
library(sf)
library(magrittr)
library(lme4)

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
               month %in% c(12,1,2) ~ "Winter"),
      logConc = log(Concentration)
    ) %>% 
    filter(!is.na(logConc))
 
# Visualise   
  # Concentration will vary between sites.
  
  RSN_MS$season <- factor(RSN_MS$season, levels = c("Winter", "Spring", "Summer", "Autumn"))
  
  ggplot(RSN_MS, aes(x = season, y =logConc, fill=season))+
    geom_violin() + geom_boxplot(width=0.3)+labs(title = "Concentraion by Season")

  pal <- MetBrewer::met.brewer("Peru1", type="continuous")
  
  leaflet::leaflet(RSN_MS) %>% 
    addProviderTiles(providers$Esri) %>% 
    addCircleMarkers(fillColor = pal, 
                     color = NA,
                     radius = 3,
                     fillOpacity = 1) 
    
  leaflet::leaflet(RSN_MS) %>% 
    addProviderTiles(providers$Esri) %>% 
    leaflet.extras::addHeatmap(intensity =~logConc,
                               blur=30, radius=20)
 
  
  
  
  lme <- 
    lmer(logConc ~ month + (1 | Sample_Site_ID), data = RSN_MS)
  lme
  
  coef(lme)$Sample_Site_ID
  # Fix the 
  lmef <- lmer(logConc ~ month + (month | Sample_Site_ID), data=RSN_MS)
  
  # Fix season, but for each site allow seasonality to vary.
  summary(lmef)
  coef(lmef)$Sample_Site_ID
   # slope coefficients are the same per site for each month.
 # Test to see if temporal variability is statistically significant and see how it compares to spatial variability 
 # Site is our random effect
 
  
  
  
  lme1 <- 
    lmer(logConc ~ (1 | Sample_Site_ID), data = RSN_MS)
  lme1
 
  # Add seasonality as an random effect instead of just spatial. 
 lme2 <- lme4::lmer(logConc ~ (1| Sample_Site_ID) + (1 + season), data = RSN_MS)
 lme2
 # season and sites are both random. How much of variation in conc is due to season or temp
 # average seasonal difference across sites but preventing them from varying by site.
 # Assumes all sites respond to seasons in the same way
 
 lme2 <- lme4::lmer(logConc ~ (1| Sample_Site_ID) + (1 | season), data = RSN_MS)
 lme2
 # (1| between site variation), (1| between season variation) residuals are within-site or within-season unexplained variation
 
 summary(lme2)
 
 # spring and winter have notably lower concs -0.30 ish whilst summer is 0.03 different to autumn. 
 
 anova(lme1,lme2)
 # seasonality effect is statistically significant. AIC difference >10 so lme1 better fit.
 
 # Try seasonality as random effect 
 
 # Add seasonality as a random effect instead of just spatial. 
 lme2r <- lme4::lmer(logConc ~ (1| Sample_Site_ID) + (1 | season), data = RSN_MS)
 lme2r  # can remove the 1 here
 
 # residuals explains more than spatial or temporal - what we can't account for. 
 # Spaital explains ~double the amount of variance than temporal. 
 # lets try month instead of season to bring the residuals down and explain the intra-season varation better.
 # explains marginally more. So seasons is good enough, our residuals still stay high when 
 # we account for variance between months.
 
 lme3 <- lme4::lmer(logConc ~ (1| Sample_Site_ID) + (1 | month), data = RSN_MS)
 lme3
 
 
 # Resid 1.693
 # concs rise in summer/ autumn months.
 
 # Above only allowed random intercept but fixes slopes, lets allow random slopes to better
 # capture site-specific season response.
 
 lme4 <- lme4::lmer(logConc ~ season + (season|Sample_Site_ID), data= RSN_MS)
 lme4
 # doesn't bring residuals down.
 
 
 lme5 <- lme4::lmer(logConc ~ season * year + (1 | Sample_Site_ID), data= RSN_MS)
 # fix year for the 4 years we have as only have a few. 
 
 # year doesn't do much - there's more variance between sites and seasons than seasons and years.
 
 anova(lme4,lme5)
 
 # lets do between catchments.
 
 
 
 lme6 <- lme4::lmer(logConc ~ season + (1 | Sample_Site_ID) + (1| ARE_DESC), data= RSN_MS)
 
 anova(lme4, lme6)
 
 lme7 <- lme4::lmer(logConc ~ season + (1| Sample_Site_ID), data=RSN_MS)
 
 anova(lme6, lme7)
 # 
 
 # Should probably group by determinand group, like pesticide, PFAS etc.
 
 lme8 <- lme4::lmer(logConc ~ (1| Sample_Site_ID) + (1| season), data=RSN_MS)
 
 
 plot(fitted(lme8), resid(lme8))
 abline(h = 0, col="red")
 # some heteroscadacity 
 
 qqnorm(resid(lme8))
 qqline(resid(lme8), col="red")
 
 
 
 
 
 lme9 <- lme4::lmer(logConc ~ (1| Sample_Site_ID) + (1| season) + (1|OPCAT_NAME), data=RSN_MS)
 # ICC after 
 
 # Use ICC to look at intragroup variability and compare against intergroup variability.
 
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
 
    
    #Q - IDENTIFY IF BIMONTHLY OR SEASONALLY IS STATITICALLY SIGNIFICANT.
    # Q - temporal variability isn't 0 so has some impact, could sample some sites montlhy and others seasonally.
    #       filter out monitoring sites to be bi-monthly to simulate this monitoring, what effect does it have?
    
    
    
    
    
    
    
    
    
    
    
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
  

  
 
  
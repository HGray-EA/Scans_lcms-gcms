library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(dplyr)
library(brickster)


# Load in data and transform 

# Function which reads tsv
readtsv <- function(filepath){
  
  read.table(paste0("/dbfs/mnt/lab/unrestricted/harry.gray@environment-agency.gov.uk/SEDA_Analysis/",filepath), sep = '\t', header =T)
}
GCMS_Open_Tgt <- readtsv("GCMS_Target_and_Non_Targeted_Screening_for_Open_Data__detes_3106_7299_4994__2025_08_08.tsv")

LCMS_Open_Tgt <- readtsv("LCMS_Target_and_Non_Targeted_Screening_for_Open_Data__detes_3106_7299_4994__2025_08_08.tsv")  

RSN_Sites <- readxl::read_excel("/dbfs/mnt/lab/unrestricted/harry.gray@environment-agency.gov.uk/SEDA_Analysis/Site_locations__showing_current_status_.xlsx", skip =2)

# Transform 
# names(LCMS_Open_Tgt)==names(GCMS_Open_Tgt) check colnames are the same before rbind

# Merge LCMS & GCMS data into one dataframe
Open_Tgt <- rbind(LCMS_Open_Tgt, GCMS_Open_Tgt)

# Convert to sf for leaflet map later
Open_Tgt_sf <- Open_Tgt %>% st_as_sf(coords = c("Longitude","Latitude"), crs=4326)
RSN_Sites_sf <- RSN_Sites %>% st_as_sf(coords = c("GRTS longitude (WGS84)","GRTS latitude (WGS84)"), crs=4326)

# Confirm if site ids match
# Extract month and year for temporal variance and log concentration for normality.

# Select scans sites from RSN network only  
RSN_MS <-  Open_Tgt_sf %>% 
  filter(str_detect(Sample_Site_ID, "RSN"))

RSN_MS %<>% 
  mutate(
    month = as.factor(month(Sample_datetime)),
    season = case_when(
      month %in% c(3,4,5) ~ "Spring",
      month %in% c(6,7,8) ~ "Summer",
      month %in% c(9,10,11) ~ "Autumn",
      month %in% c(12,1,2) ~ "Winter"),
    logConc = log(Concentration)
  )

# Order for plotting 
RSN_MS$season <- factor(RSN_MS$season, levels = c("Winter", "Spring", "Summer", "Autumn"))

#-----------------------------------------------------------------------------#

ui <- page_sidebar(
  title = "Scans Compound Seasonality",
  sidebar = sidebar(
        width = 300,
    selectInput(
      "compound",
      "Select Compound:",
      choices = unique(RSN_MS$Compound_Name),
      selected = unique(RSN_MS$Compound_Name)[1]
    )
  ),
  card(
    card_header("Concentration by Season"),
    plotlyOutput("concentration_plot", height = "500px")
  )
)

#-----------------------------------------------------------------------------#

server <- function(input, output, session) {
# Reactive data filtering based on selected compound
  filtered_data <- reactive({
    RSN_MS %>%
      filter(Compound_Name == input$compound)
  })
  
  output$concentration_plot <- renderPlotly({
   
# Create the ggplot
  plt <-   ggplot(filtered_data(), aes(season, logConc, fill = season)) +
      ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, fill = "skyblue") +
      geom_boxplot(width = .15, outlier.shape = NA) +
      geom_jitter(width = .1, alpha = 0.5)+
      labs(y = "log(Concentration [ug/L])", x= "Season")
    
    
# Convert to plotly
    ggplotly(plt)
  })
}

#-----------------------------------------------------------------------------#
shinyApp(ui = ui, server = server)

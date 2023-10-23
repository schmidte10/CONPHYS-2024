#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#--- import libraries ---#
library(readr)
library(tidyverse)
library(chron)
library(lubridate)
library(shiny)
library(ggpubr)
library(shinythemes)
library(DT)
#--- import data ---#
cs <- read_delim("CS_LocalAdapt6.txt", 
                 delim = "\t", escape_double = FALSE, 
                 col_types = cols(...21 = col_skip(), 
                                  ...22 = col_skip()), trim_ws = TRUE) %>% 
  mutate(creation_time = as.POSIXct(`Creation time`, format = "%d/%m/%Y %H:%M:%S"))

cs2 <- cs %>%
  unite("UNIQUE_SAMPLE_ID", c(FISH_ID,TEMPERATURE,`Sample index`), sep="_", remove = FALSE) %>% 
  separate(creation_time, into=c('DATE','TIME'), sep = " ", remove = FALSE) %>% 
  arrange(`Sample ID 1`, DATE, TIME) 

cs3 <- cs2 %>% 
  mutate(DATE = as.Date(creation_time), 
         TIME = as.POSIXct(creation_time, "%H:%M:%S")) %>% 
  mutate(TIME = chron(times=cs2$TIME)) %>% 
  arrange(TIME) %>%
  group_by(UNIQUE_SAMPLE_ID, `Sample ID 1`) %>% 
  mutate(TIME_DIFF = TIME - first(TIME)) %>% 
  filter(TIME != first(TIME)) %>%
  ungroup() %>% 
  mutate(TIME_DIFF_SECS = period_to_seconds(hms(TIME_DIFF))) %>% 
  mutate(MINUTES = TIME_DIFF_SECS/60) %>% 
  mutate(MINUTES = round(MINUTES, digits = 2)) %>% 
  dplyr::rename(CUVETTE = `Sample ID 1`) %>% 
  mutate(REGION = substr(FISH_ID, 1, 1 ), 
         POPULATION = substr(FISH_ID, 2, 4), 
         SAMPLE_NO = substr(FISH_ID, 5, 7)) %>% 
  mutate(REGION = case_when( REGION =="L"~ "Leading", 
                             REGION == "C" ~ "Core", 
                             TRUE ~ "na")) 

cs3$CUVETTE <- paste("Cuvette", cs3$CUVETTE, sep="_") 
cs3 <- cs3 %>% 
    mutate(CUVETTE = as.factor(CUVETTE)) 

cs3 <- cs3 %>%
  group_by(POPULATION) %>%
  mutate(SAMPLE_NO = dense_rank(FISH_ID)) %>% 
    ungroup() 


# Define UI for application that draws a histogram
ui <- fluidPage(
  #tags$head(
    #tags$style(type="text/css", "select { max-width: 240px; }"),
    #tags$style(type="text/css", ".span4 { max-width: 290px; }"),
    #tags$style(type="text/css", ".well { max-width: 280px; }")
  #),
    # Application title
  navbarPage("Enzyme Quailty Checks", theme = shinytheme("yeti"),
             tabPanel("CS", fluid = TRUE, icon = icon("fish"), 
                      
  

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            div(style="display:inline-block;",radioButtons(
                "TEMPERATURE",
                        "Temperature C",
                        c("20" = "20",
                          "30" = "30",
                          "40" = "40",
                          "50" = "50"))),  
            
            selectInput("POPULATION", "POPULATION:", 
                        c("Cockermouth Island" = "CKM",
                          "Keswick Island"= "KES",
                          "Chauvel Reef"= "CHA",
                          "Sudbury Reef"= "SUD",
                          "Tongue Reef"= "TON",
                          "Vlassof Cay"= "VLA")),
            selectInput("SAMPLE_NO", "Replicate:", 
                        c("1"="1", 
                          "2"="2",
                          "3"="3",
                          "4"="4",
                          "5"="5",
                          "6"="6",
                          "7"="7",
                          "8"="8",
                          "9"="9",
                          "10"="10"))
            ), 
        
          
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("plot"), 
           DT::dataTableOutput("mytable")) 
        )
    ))
  )
    


# Define server logic required to draw a histogram
server <- function(input, output, session) { 
  
  CS_LocalAdapt3_Finder <- reactive({ 
    req(input$TEMPERATURE)  
    req(input$POPULATION)
    filter(cs3, TEMPERATURE %in% input$TEMPERATURE,
           POPULATION %in% input$POPULATION,
           SAMPLE_NO %in% input$SAMPLE_NO)}) 

  # suppress warnings  
  storeWarn<- getOption("warn")
  options(warn = -1)
  
output$plot <- renderPlot({ 
  input$TEMPERATURE
  input$POPULATION 
  input$SAMPLE_NO
  isolate({   
    ggplot(CS_LocalAdapt3_Finder(),aes(MINUTES, Result)) + 
      geom_point() +
      facet_wrap(~CUVETTE) +
      geom_smooth(method = "lm") +
      theme_bw() + 
      ylim(-0.3,1.4)+
      ggtitle(paste(CS_LocalAdapt3_Finder()[1,8])) + 
      stat_regline_equation(label.y = 0.7) + 
      stat_cor(label.y = 0.6)})
}, height = 420, width = 700)  

output$mytable = DT::renderDataTable({ 
  input$TEMPERATURE 
  input$POPULATION 
  input$SAMPLE_NO
  
  CS_activity <- CS_LocalAdapt3_Finder() %>% 
    group_by(UNIQUE_SAMPLE_ID, CUVETTE) %>% 
    do({
    mod = lm(Result ~ MINUTES, data = .)
      data.frame(Intercept = coef(mod)[1],
                 Slope = coef(mod)[2], 
                 r2 = summary(mod)$r.squared)
    }) %>%
    ungroup() %>%
    datatable() %>%
    formatStyle('CUVETTE', target = "row",
                backgroundColor = styleEqual(c("Cuvette_1","Cuvette_2","Cuvette_3"), c('lightblue','lightblue','lightblue'))) %>% 
    formatStyle('CUVETTE', target = "row",
                backgroundColor = styleEqual(c("Cuvette_5"), c('springgreen')))
  
}) 

}


  

# Run the application 
shinyApp(ui = ui, server = server)

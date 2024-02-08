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
LDH_LocalAdapt <- read_delim("LDH_LocalAdapt.txt", 
                             delim = "\t", escape_double = FALSE, 
                             col_types = cols(`Creation time` = col_datetime(format = "%d/%m/%Y %H:%M")), 
                             trim_ws = TRUE)



#--- separate sample column into useful data ---# 
LDH_LocalAdapt2 <-  
    LDH_LocalAdapt %>% 
    unite("UNIQUE_SAMPLE_ID", c(FISH_ID,TEMPERATURE), sep = "_", remove = FALSE) %>% 
    mutate(Region = str_sub(FISH_ID, 1, 1),
         Population = str_sub(FISH_ID, 2, 4),
         ID = str_sub(FISH_ID, 5)) %>% 
    dplyr::rename(CUVETTE = `Sample index`, 
           DATETIME = `Creation time`) %>% 
    separate(DATETIME, into=c('DATE','TIME'), sep = " ", remove = FALSE) %>%  
    arrange(CUVETTE, TIME)

LDH_LocalAdapt3 <- LDH_LocalAdapt2 %>% 
    mutate(TIME = chron(times = LDH_LocalAdapt2$TIME)) %>%
    dplyr::group_by(UNIQUE_SAMPLE_ID, CUVETTE) %>% 
    mutate(TIME_DIFF = TIME - first(TIME)) %>%  
    ungroup() %>%
    mutate(TIME_DIFF_SECS = period_to_seconds(hms(TIME_DIFF))) %>% 
    mutate(MINUTES = TIME_DIFF_SECS/60) %>% 
    mutate(MINUTES = round(MINUTES, digits = 2)) %>% 
    dplyr::rename(POPULATION = Population)

LDH_LocalAdapt3$CUVETTE <- paste("Cuvette", LDH_LocalAdapt3$CUVETTE, sep="_") 
LDH_LocalAdapt3 <- LDH_LocalAdapt3 %>% 
    mutate(CUVETTE = as.factor(CUVETTE)) 

LDH_LocalAdapt3 <- LDH_LocalAdapt3 %>%
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
             tabPanel("LDH", fluid = TRUE, icon = icon("fish"), 
                      
  

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
  
  LDH_LocalAdapt3_Finder <- reactive({ 
    req(input$TEMPERATURE)  
    req(input$POPULATION)
    filter(LDH_LocalAdapt3, TEMPERATURE %in% input$TEMPERATURE,
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
    ggplot(LDH_LocalAdapt3_Finder(),aes(MINUTES, Result)) + 
      geom_point() +
      facet_wrap(~CUVETTE) +
      geom_smooth(method = "lm") +
      theme_bw() + 
      ylim(-0.3,1.4)+
      ggtitle(paste(LDH_LocalAdapt3_Finder()[1,8])) + 
      stat_regline_equation(label.y = 0.7) + 
      stat_cor(label.y = 0.6)})
}, height = 420, width = 700)  

output$mytable = DT::renderDataTable({ 
  input$TEMPERATURE 
  input$POPULATION 
  input$SAMPLE_NO
  
  LDH_activity <- LDH_LocalAdapt3_Finder() %>% 
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

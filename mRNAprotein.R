# Packages ----
library(shiny)  # Required to run any Shiny app
library(ggplot2)  # For creating pretty plots
library(dplyr)  # For filtering and manipulating data
library(tidyr)
library(tidyverse)
library(data.table)
library(mice)
library(corrplot)
library("FactoMineR")
library("factoextra")

source("shiny-app-functions.R")


# Define UI ----------
ui <- fluidPage(
  titlePanel("mRNA and Protein Expression"),
  
  navlistPanel(
    "By Tissue Type",
    tabPanel("Plots", 
             sidebarLayout(position="left", #default = left
                           sidebarPanel(
                             selectInput(inputId = "tissue",  # Give the input a name "tissue"
                                         label = "1. Select tissue type",  # Give the input a label to be displayed in the app
                                         choices = unique(colnames(mergedData[2:6])), selected="uterus"),
                           ),
                           mainPanel(
                             textOutput("plotTitle"),
                             plotOutput("tissuePlot"))
                          )
             ),
    tabPanel("Correlations", 
             fluidRow(
               column(width = 6,
                      plotOutput("corrPlot"))
               )),
    tabPanel("PCA Analysis",
             fluidRow(
               column(width = 6,
                      plotOutput("screePlot")
               ),
               column(width = 6,
                      plotOutput("cosPlot")),
               column(width = 6,
                      plotOutput("corrCircle"))
    )),
    
    "By Gene",
    tabPanel("Correlations across genes", plotOutput("boxPlot")),
    tabPanel("Correlations per genes")
  ),
)


# Define server logic ---------           
server <- function(input, output, session) {
  
  # For tissue type plots
  output$plotTitle <- renderText(input$tissue)
  output$tissuePlot <- renderPlot({  #NOTE: used if else statements display correct plot
    tissue_plot(input$tissue)
  })
  
  #For tissue type correlations
  output$corrPlot <- renderPlot({
    corrplot(spear_corr, method = "square")
  }) 
  output$boxPlot <- renderPlot({
    plotFigureA(orig_data_across_corrs,
                orig_predict_across_corrs,
                free_predict_across_corrs,
                rand_predict_across_corrs)
  }) 
  
  #For tissue type PCA
  output$screePlot <- renderPlot({
    scree_plot()
  })
  output$cosPlot <- renderPlot({
    cos_plot()
  })
  output$corrCircle <- renderPlot({
    corr_circle()
  })
  
}

# Run the app -----------
shinyApp(ui = ui, server = server)

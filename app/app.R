#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Final Project BF591 Pooja Savla

library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(magrittr)
library(tidyverse)
library(glue)

# UI
ui <- fluidPage(
  tags$head(
    tags$style(
      HTML("
        .nav-tabs > li > a {
          background-color: #ccc;
          width: 2in;
          margin-right: 10px;
          text-align: center;
          color: black;
        }
        .nav-tabs {
          display: flex;
          justify-content: flex-start;
        }
      ")
    )
  ),
  h1("My RShiny Project"),
  tabsetPanel(
    tabPanel("Tab 1",
             sidebarLayout(
               sidebarPanel(
                 h3("Sidebar Panel for Tab 1"),
                 # Add your sidebar content here
               ),
               mainPanel(
                 h2("Main Panel for Tab 1"),
                 tabsetPanel(
                   tabPanel("Nested Tab 1",
                            h3("Content for Nested Tab 1")
                            # Add your content for Nested Tab 1 here
                   ),
                   tabPanel("Nested Tab 2",
                            h3("Content for Nested Tab 2")
                            # Add your content for Nested Tab 2 here
                   ),
                   tabPanel("Nested Tab 3",
                            h3("Content for Nested Tab 3")
                            # Add your content for Nested Tab 3 here
                   )
                 )
               )
             )
    ),
    tabPanel("Tab 2",
             sidebarLayout(
               sidebarPanel(
                 h3("Sidebar Panel for Tab 2"),
                 # Add your sidebar content here
               ),
               mainPanel(
                 h2("Main Panel for Tab 2"),
                 p("Content for Tab 2 goes here.")
                 # Add your main panel content here
               )
             )
    ),
    tabPanel("Tab 3",
             sidebarLayout(
               sidebarPanel(
                 h3("Sidebar Panel for Tab 3"),
                 # Add your sidebar content here
               ),
               mainPanel(
                 h2("Main Panel for Tab 3"),
                 p("Content for Tab 3 goes here.")
                 # Add your main panel content here
               )
             )
    ),
    tabPanel("Tab 4",
             sidebarLayout(
               sidebarPanel(
                 h3("Sidebar Panel for Tab 4"),
                 # Add your sidebar content here
               ),
               mainPanel(
                 h2("Main Panel for Tab 4"),
                 p("Content for Tab 4 goes here.")
                 # Add your main panel content here
               )
             )
    )
  )
)

# Server
server <- function(input, output) {
  # Server logic goes here (if needed)
  load_data <- reactive({
    req(input$file)
    data <- read.csv(input$file$datapath)
    colnames(data)[1] <- "gene"
    return(data)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
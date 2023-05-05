## Author: Pooja Savla
## prj4998@bu.edu
## BU BF591
## Final Project



# install.packages('shinyHeatmaply')
# install.packages("devtools")
# devtools::install_github("thomasp85/patchwork")
# install.packages("gplots")
# install.packages("shinythemes")
library(shiny)
library(DT)
library(ggbeeswarm)
library(plotly)
library(gplots)
library(qwraps2)
library(heatmaply)
library(shinythemes)
library(patchwork)
library(tidyverse)
library(bslib)
library(ggplot2)
library(plotly)
library(colourpicker)

# UI
ui <- fluidPage(
  theme = shinytheme("cosmo"),
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
          overflow-x: auto;
        }
      ")
    )
  ),
  h1("My RShiny Project - BF591 Final Project - Pooja Paresh Savla"),
  h6("mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington's Disease and neurologically normal individuals"),
  tabsetPanel(
    # ---------------------------------------------------------Samples Tab------------------------------------------------------------
    tabPanel(
      "Samples",
      sidebarLayout(
        sidebarPanel(
          fileInput(
            "sample_metadata",
            label = "Load sample metadata",
            placeholder = "sample_metadata.csv"
          ),
          div(
            actionButton(
              inputId = "Submit", 
              label = "Submit", 
              style = "color: white; background-color:  #006994;"
            ), 
            style = "text-align: center;"
          )
        ),
        mainPanel(
          h2(),
          tabsetPanel(
            tabPanel("Summary", 
                     tableOutput("samples_summary")),
            tabPanel("Table", 
                     dataTableOutput("samples_table")),
            tabPanel("Plots", 
                     plotOutput("samples_plots"))
          )
        )
      )
    ),
    # ---------------------------------------------------------Counts Tab------------------------------------------------------------
    tabPanel(
      "Counts",
      sidebarLayout(
        sidebarPanel(
          fileInput("count_file", 
                    label = "Load sample normalized counts",
                    placeholder = "norm_counts.csv"),
          sliderInput("variance_slider", 
                      "Select percentile of variance", 
                      min = 0, 
                      max = 100, 
                      value = 50, 
                      step = 1),
          sliderInput("non_zero_slider", 
                      "Select the number of the samples that are non-zero:", 
                      min = 0, 
                      max = 100, 
                      value = 50, 
                      step = 1),
          div(
            actionButton(
              inputId = "Submit", 
              label = "Submit", 
              style = "color: white; background-color:  #006994;"
            ), 
            style = "text-align: center;"
          )
        ),
        mainPanel(
          h2("Main Panel for Tab 2"),
          tabsetPanel(
            tabPanel("Filter Summary", 
                     tableOutput("counts_summary")),
            tabPanel("Diagnostic Plots", 
                     plotOutput("scatter_plots")),
            tabPanel("Heatmap", 
                     plotOutput("counts_heatmap")),
            tabPanel("PCA",  
                     sidebarLayout(
                       sidebarPanel(
                         radioButtons("PC1_choice", 
                                      "Choose the PC for the x-axis", 
                                      choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
                                      selected = "PC1"),
                         radioButtons("PC2_choice", 
                                      "Choose the PC for the y-axis", 
                                      choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
                                      selected = "PC2"),
                         submitButton("Plot", width='100%')
                       ),
                       mainPanel(
                         plotOutput("counts_pca_plot")
                       )
                     )

          )
        )
      )
    )
    ),
# ---------------------------------------------------------differential expression results Tab------------------------------------------------------------
    tabPanel(
      "DE",
      sidebarLayout(
        sidebarPanel(
          fileInput("file", "Load differential expression results", accept = ".csv", placeholder = "deseq_res.csv"),
          HTML("A volcano plot can be generated with <b>\"log2 fold-change\"</b> on the x-axis and <b>\"p-adjusted\"</b> on the y-axis."),
          radioButtons("x_axis", "Choose the column for the x-axis",
                       choices = c("baseMean" = "baseMean","log2FoldChange" = "log2FoldChange", 
                                   "lfcSE"="lfcSE","stat"="stat","pvalue"="pvalue","padj"="padj"),
                       selected = "log2FoldChange"),
          radioButtons("y_axis", "Choose the column for the y-axis",
                       choices = c("baseMean" = "baseMean","log2FoldChange" = "log2FoldChange", 
                                   "lfcSE"="lfcSE","stat"="stat","pvalue"="pvalue","padj"="padj"),selected = "padj"),
          colourInput("base", "Base point color", value = "#22577A"),
          colourInput("highlight", "Highlight point color", value = "#FFCF56"),
          sliderInput(inputId = "slider", min = -300, max = 0,
                      label = "Select the magnitude of the p adjusted coloring:", value = -150, step = 1),
          div(actionButton(inputId = "Plot", label = "Plot", style = "color: green; background-color: #78c2ad;"), style = "text-align: center;")
        ),
        
        mainPanel(
               tabsetPanel(
                 type = "tabs",
                 tabPanel("Plot", plotOutput("volcano")),
                 tabPanel("Table", tableOutput("table"))
               )
             )
           )
         ), 
    tabPanel(
           "GSEA",
           sidebarLayout(
             sidebarPanel(
               h3("Sidebar Panel for Tab 4"),
               fileInput("file", "Load FGSEA sample file", accept = ".csv", placeholder = "fgsea.csv")
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
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024^2)  # increase file upload limit

# ---------------------------------------------------------Samples Tab------------------------------------------------------------  
  samples_load_data <- reactive({
    inFile <- input$sample_metadata
    
    if (is.null(inFile)) {
      return(NULL)
    }
    df <- read.csv(file = inFile$datapath, header = TRUE, stringsAsFactors = TRUE)
    
    return(df)
  })
  
  draw_sum_table <- function(dataf) {
    column_names <- colnames(dataf)
    col_type <- sapply(dataf, class)
    mean_sd_col <- sapply(dataf, function(col) {
      if (is.factor(col)) {
        paste(unique(col), collapse = ", ")
      } else {
        mean <- round(mean(col), 2)
        sd <- round(sd(col), 2)
        paste(mean, "(+-", sd, ")")
      }
    })
    
    df <- data.frame("Column Name" = column_names,
                     "Type" = col_type,
                     "Mean (sd) or Distinct Values" = mean_sd_col)
    
    return(df)
  }
  
  observeEvent(input$Submit, {
    req(input$sample_metadata)
    table <- samples_load_data()
    output$samples_summary <- renderTable({
      draw_sum_table(dataf = table)
    }, striped = TRUE)
  })
  
  output$samples_table <- renderDataTable({
    samples_load_data()
  })

  output$samples_plots <- renderPlot({
    samples_plot(samples_load_data())
  })
  
  #---------------------------------------------------------Counts Tab------------------------------------------------------------
  counts_load_data <- reactive({
    inFile <- input$count_file
    if (is.null(inFile)) {
      return(NULL)
    }
    df <- read.csv(file = inFile$datapath, header = TRUE, stringsAsFactors = TRUE)
    return(df)
  })
  
  draw_count_sum_table <- function(count_data, var_val, zero_val) {
    # number of samples
    num_samples <- NCOL(count_data)
    # total number of genes
    total_genes <- length(count_data$gene)
    
    # Filter data based on variance percentile and nonzeros
    nonzero_genes <- rowSums(count_data[-1] == 0) >= zero_val
    nonzero_counts <- count_data[nonzero_genes, ]
    zero_counts <- filter(count_data, !gene %in% nonzero_counts$gene)
    nonzero_counts$variance <- apply(nonzero_counts[, -c(1)], 1, var)
    
    # Calculation percentile
    percent <- quantile(nonzero_counts$variance, prob = var_val / 100)
    
    var_counts <- dplyr::filter(nonzero_counts, variance >= percent) %>%
      dplyr::select(-variance)
    
    # number and % of genes passing current filter
    num_pass <- length(var_counts$gene)
    n1 <- round((num_pass / total_genes) * 100, 2)
    pass <- c(num_pass, n1)
    num_pass <- toString(pass, width = 30)
    
    # number and % of genes not passing current filter
    np <- length(zero_counts$gene) + (length(nonzero_counts$gene) - length(var_counts$gene))
    n2 <- round(((np / total_genes) * 100), 2)
    not_pass <- c(np, n2)
    
    num_not_pass <- toString(not_pass, width = 30)
    
    df <- data.frame(num_samples, total_genes, num_pass, num_not_pass)
    colnames(df) <- c(
      "Number of Samples", "Total Number of Genes",
      "Number , % of genes passing the filter",
      "Number , % of genes not passing the filter"
    )
    return(df)
  }
  
  output$counts_summary <- renderTable(
    {
      req(input$count_file)
      table <- counts_load_data()
      if (!is.null(table)) {
        draw_count_sum_table(
          table,
          input$variance_slider,
          input$non_zero_slider
        )
      }
    },
    spacing = "m",
    bordered = T
  )

  
}


shinyApp(ui, server)
      
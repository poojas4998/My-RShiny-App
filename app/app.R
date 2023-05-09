## Author: Pooja Savla
## prj4998@bu.edu
## BU BF591
## Final Project



# install.packages('shinyHeatmaply')
# install.packages("devtools")
# devtools::install_github("thomasp85/patchwork")
# install.packages("gplots")
# install.packages("shinythemes")
# install.packages("shinythemes")
library(shiny)
library(DT)
library(plotly)
library(gplots)
library(qwraps2)
library(heatmaply)
library(shinythemes)
library(patchwork)
library(tidyverse)
library(bslib)
library(ggplot2)
library(colourpicker)
library(grid)
library(pheatmap)
library(gridExtra)
library(glue)
library(pheatmap)
library(ggbeeswarm)
library(GSEABase)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(magrittr)

ui <- fluidPage(
  # Use the "cosmo" theme for the app
  theme = shinytheme("cosmo"),
  # Define custom styles for the navigation tabs
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
  # Define the app title and subtitle
  h1("My RShiny Project - Pooja Paresh Savla"),
  h6("mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington's Disease and neurologically normal individuals"),
  # Define the tabset panel
  tabsetPanel(
    # Define the "Samples" tab
    tabPanel(
      "Samples",
      # Define the sidebar layout for the "Samples" tab
      sidebarLayout(
        # Define the sidebar panel for the "Samples" tab
        sidebarPanel(
          # Define a file input for the sample metadata
          fileInput(
            "sample_metadata",
            label = "Load sample metadata",
            placeholder = "sample_metadata.csv"
          ),
          # Define a submit button to submit the metadata file
          div(
            actionButton(
              inputId = "submit_1", 
              label = "Submit", 
              style = "color: white; background-color:  #006994;"
            ), 
            style = "text-align: center;"
          )
        ),
        # Define the main panel for the "Samples" tab
        mainPanel(
          # Define the summary table output for the "Summary" tab
          tabsetPanel(
            tabPanel("Summary", 
                     tableOutput("samples_summary")),
            # Define the table output for the "Table" tab
            tabPanel("Table", 
                     dataTableOutput("samples_table")),
            # Define the histogram plot output for the "Histogram Plots" tab
            tabPanel("Histogram Plots", 
                     fluidRow(
                       br(),
                       radioButtons("variable_choice", 
                                    "Choose the variable for the x-axis", 
                                    choices = c("Age_of_death", "PMI", "RIN"),
                                    selected = "Age_of_death",
                                    inline = TRUE),
                       submitButton("Plot")
                     ),
                     fluidRow(
                       mainPanel(
                         plotOutput("samples_histogram")
                       )
                     )
            )
          )
        )
      )
    ),
    # ---------------------------------------------------------Counts Tab------------------------------------------------------------
    tabPanel(
      "Counts",
      sidebarLayout(
        sidebarPanel(
          fileInput("counts_csvfile",
                    label = "Upload counts matrix file.", 
                    accept = ".csv", 
                    placeholder = "counts.csv"),
          sliderInput("counts_slider_var",
                      min = 0, 
                      max = 100,
                      label = "Select a minimum percentile of variance per gene:", 
                      value = 65, 
                      step = 1),
          sliderInput("counts_slider_num",
                      min = 0, 
                      max = 70,
                      label = "Select a minimum number of non-zero samples per gene:", 
                      value = 60, 
                      step = 1),
          submitButton("Plot", 
                       width = '100%')
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Filter Summary", 
                     tableOutput("counts_table")),
            tabPanel("Diagnostic Plots", 
                     plotOutput("counts_filter_plots")),
            tabPanel("Heatmap", 
                     plotOutput("counts_heatmap")),
            tabPanel("PCA",
                     sidebarLayout(
                       sidebarPanel(
                         selectInput("PC1_choice",
                                     "Choose the PC for the x-axis",
                                     choices = paste0("PC", 1:69),
                                     selected = "PC1"),
                         selectInput("PC2_choice",
                                     "Choose the PC for the y-axis",
                                     choices = paste0("PC", 1:69),
                                     selected = "PC2"),
                         submitButton("Plot", width='100%')
                       ),
                       mainPanel(
                         plotOutput("counts_PCA")
                       )
                     )
            )
            
          )
        )
      )
    ),
    tabPanel("Differential Expression",
             p("The differential expression tab allows the user to explore the results of differential expression on uploading a differential expression analysis results."),
             sidebarLayout(
               sidebarPanel(
                 fileInput("DE_Anaylsis", "Choose a Differential Expression Results CSV File", accept = ".csv",placeholder = 'Differential Expression Data File'),
                 submitButton(text = 'Submit'), width = 3
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Table",
                            p("The table tab provides sortable table of the differential expression analysis results."),
                            DT::dataTableOutput("table",width = "80%")),
                   tabPanel("Volcano Plot", 
                            p("The volcano tab provides a volcano plot based on the input paramters, namely, X-axis and Y-axis"),
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("x_axis", "X-Axis",
                                             choices = c("log2FoldChange" = "log2FoldChange","baseMean" = "baseMean", "lfcSE"="lfcSE","stat"="stat","pvalue"="pvalue","padj"="padj"),
                                             selected = "log2FoldChange"),
                                radioButtons("y_axis", "Y-Axis",
                                             choices = c("log2FoldChange" = "log2FoldChange","baseMean" = "baseMean", "lfcSE"="lfcSE","stat"="stat","pvalue"="pvalue","padj"="padj"),
                                             selected = "padj"),
                                sliderInput(inputId = "padjusted", min = -35, max = 0,label = "Select the magnitude of the p adjusted coloring:", value = -15, step = 1),
                                colourInput("de_base", "Base point color","pink"),
                                colourInput("de_highlight", "Highlight point color","darkgreen"),
                                submitButton("Submit"),width=3
                              ),
                              mainPanel(
                                plotOutput("volcano"))))
                 )
               )
             )
    ), # tab panel bracket
    tabPanel("GSEA",
             p("Use this tab to visualize gene set enrichment analysis results of the data."),
             sidebarLayout(
               sidebarPanel(
                 fileInput("GSEA_csvfile",
                           label = "Upload GSEA results table.", 
                           accept = ".csv", 
                           placeholder = "gsea_results.csv")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("NES Bar Plot",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("GSEA_Bar",
                                            min = 1, 
                                            max = 50,
                                            label = "Select number of top pathways:", 
                                            value = 40, 
                                            step = 1),
                                submitButton("Plot", 
                                             width = '100%')
                              ),
                              mainPanel(
                                plotOutput("GSEA_barplot")
                              )
                            )),
                   tabPanel("Table",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("GSEA_table_pvalue",
                                            min = -20, 
                                            max = 0,
                                            label = "Select adjusted p-value to filter by:", 
                                            value = -1, 
                                            step = 1),
                                radioButtons("GSEA_pathways_choice",
                                             "Choose the types of pathways", 
                                             choices = c("All", "Positive", "Negative"),
                                             selected = "All"),
                                submitButton("Submit", width="100%")
                              ),
                              mainPanel(
                                sidebarLayout(
                                  sidebarPanel(
                                    downloadButton('download_NES_table', "Download Table")
                                  ),
                                  mainPanel(
                                    dataTableOutput("GSEA_results_table")
                                  )
                                ))
                            )),
                   tabPanel("NES Scatter Plot",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("GSEA_scatter_pvalue",
                                            min = -20, 
                                            max = 0,
                                            label = "Select adjusted p-value to filter by:", 
                                            value = -1, 
                                            step = 1),
                                submitButton("Plot", width="100%")
                              ),
                              mainPanel(
                                plotOutput("GSEA_scatterplot")
                              )
                            ))
                 )
               )
             )
    ) # tab panel

    
  ) # tabset panel
) # fluidpage bracket




#input file size setting


server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  
  samples_load_data <- reactive({
    inFile <- input$sample_metadata
    
    
    if (is.null(inFile)) {
      return(NULL)
    }
    df <- read.csv(file = inFile$datapath, header = TRUE, stringsAsFactors = TRUE)
    
    return(df)
  })
  
  
  load_meta <- reactive({
    if (!is.null(input$sample_metadata)){
      meta <- read_csv(input$sample_metadata$datapath)
      return(meta)}
    else{return(NULL)}
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

  samples_histogram_graph <- function(meta_tib, variable_choice = Age_of_death) {
    if (!is.null(input$sample_metadata)){
          histo <- ggplot(meta_tib, aes(x=!!sym(variable_choice)))+
            geom_histogram(bins = 10, color = "black", fill = "pink")+
            labs(title = paste("Histogram of",variable_choice))+
            xlab(paste(variable_choice))+
            ylab('Count')+
            theme_bw()
          return(histo)}
        else{return(NULL)}
  }

  
  output$samples_summary <- renderTable(
    {
      req(input$sample_metadata)
      table <- samples_load_data()
      return(draw_sum_table(dataf = table))
    },
    striped = T
  )
  
  output$samples_table <- renderDataTable({
    samples_load_data()
  })
  
  
  output$samples_histogram <- renderPlot({
    samples_histogram_graph(load_meta(), 
                            variable_choice = input$variable_choice)
  })

  
  #---------------------------------------------------------Counts Tab------------------------------------------------------------
  counts_load_data <- reactive({
    req(input$counts_csvfile)
    file=input$counts_csvfile
    if (is.null(file))
    {return(NULL)} 
    else
    {datafile=read_csv(file$datapath)} %>% 
      return()
  })
  
  counts_filter_summary <- function(dataf, slider_var, slider_num){
    stats_df <- dataf
    stats_df$variance <- apply(dataf[,-1], 1, var)
    stats_df <- arrange(stats_df, by=variance)
    stats_df$rank <- rank(stats_df$variance)
    stats_df$nonzero <- rowSums(stats_df!=0)
    filtered <- filter(stats_df, rank >= nrow(stats_df)*(slider_var/100) & nonzero >= slider_num)
    num_samples = ncol(stats_df)-4
    num_genes = nrow(stats_df)
    num_genes_filtered = nrow(filtered)
    perc_genes_filtered = (num_genes_filtered/num_genes)*100
    num_genes_out = num_genes - num_genes_filtered
    perc_genes_out = (num_genes_out/num_genes)*100
    final_stats <- tibble(num_samples=num_samples,
                          num_genes=num_genes,
                          num_genes_filtered=num_genes_filtered,
                          percent_filtered=perc_genes_filtered,
                          num_genes_filtered_out=num_genes_out,
                          percent_filtered_out=perc_genes_out)
    return(final_stats)
  }
  
  counts_filter_plots <- function(dataf, slider_var, slider_num){
    stats_df <- dataf
    stats_df$variance <- apply(dataf[,-1], 1, var)
    stats_df <- arrange(stats_df, by=variance)
    stats_df$rank <- rank(stats_df$variance)
    stats_df$nonzero <- rowSums(dataf[,-1]!=0)
    stats_df$median <- apply(dataf[,-1], 1, median)
    filtered_status <- stats_df %>% 
      mutate(filter_status=ifelse(rank >= nrow(stats_df)*(slider_var/100) & nonzero >= slider_num, 'pass filter', 'filtered out'))
    p1 <- ggplot(filtered_status, aes(x=variance, y=median, color=filter_status)) +
      geom_point() +
      scale_y_log10() +
      scale_x_log10() +
      labs(x='Percentile of Variance', y='Median Number of Counts', color='Filter Status', main='Median Counts vs. Percentile Variance') +
      theme_dark() +
      scale_color_brewer(palette="Greens")
    p2 <- ggplot(filtered_status, aes(x=nonzero, y=median, color=filter_status)) +
      geom_point() +
      scale_y_log10() +
      labs(x='Number of Non-Zero Samples', y='Median Number of Counts', color='Filter Status', main='Median Counts vs. Number of Non-Zero Samples') +
      theme_dark() +
      scale_color_brewer(palette="Oranges")
    grid.arrange(p1, p2, nrow = 1)
    return()
  }
  
  counts_heatmap <- function(dataf, slider_var, slider_num) {
    stats_df <- dataf
    stats_df$variance <- apply(dataf[,-1], 1, var)
    stats_df <- arrange(stats_df, by=variance)
    stats_df$rank <- rank(stats_df$variance)
    stats_df$nonzero <- rowSums(stats_df!=0)
    filtered <- filter(stats_df, rank >= nrow(stats_df)*(slider_var/100) & nonzero >= slider_num)
    mat <- as.matrix(filtered[,2:70])
    rownames(mat) <- filtered$gene
    p <- pheatmap(log10(mat+1),
                  scale="row",
                  color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(11),
                  cluster_rows=TRUE,
                  cluster_cols=TRUE,
                  show_rownames=FALSE,
                  show_colnames=FALSE,
                  legend=TRUE
    )
  }
  
  counts_PCA <- function(dataf, PC1_choice, PC2_choice) {
    expr_mat <- t(as.matrix(dataf[,-1]))
    pca <- prcomp(
      expr_mat,
      center=TRUE,
      scale=TRUE 
    )
    pca_var <- tibble(
      PC=factor(str_c("PC",1:69),str_c("PC",1:69)),
      Variance=pca$sdev**2,
      Explained_Variance=Variance/sum(Variance)*100,
      `Cumulative % Explained Variance`=cumsum(Explained_Variance)
    )
    as_tibble(pca$x) %>%
      ggplot(aes(x=!!sym(PC1_choice),y=!!sym(PC2_choice))) +
      geom_point() +
      labs(x=str_glue("{PC1_choice} Percent Variance Explained: {pca_var$Explained_Variance[pca_var$PC=={PC1_choice}]}"), 
           y=str_glue("{PC2_choice} Percent Variance Explained: {pca_var$Explained_Variance[pca_var$PC=={PC2_choice}]}")) %>% 
      return()
  }
  
  output$counts_table <- renderTable({counts_filter_summary(dataf = counts_load_data(), 
                                                            slider_var = input$counts_slider_var, 
                                                            slider_num=input$counts_slider_num)})
  
  output$counts_filter_plots <- renderPlot({counts_filter_plots(dataf = counts_load_data(),
                                                                slider_var = input$counts_slider_var, 
                                                                slider_num = input$counts_slider_num)})
  
  output$counts_heatmap <- renderPlot({counts_heatmap(dataf = counts_load_data(),
                                                      slider_var = input$counts_slider_var,
                                                      slider_num = input$counts_slider_num)})
  
  output$counts_PCA <- renderPlot({counts_PCA(dataf = counts_load_data(),
                                              PC1_choice = input$PC1_choice,
                                              PC2_choice = input$PC2_choice)})
  
  # ---------------------------------------------------------DE Tab------------------------------------------------------------

  DE_data <- reactive({
    req(input$DE_Anaylsis)
    data <- read.csv(input$DE_Anaylsis$datapath)
    colnames(data)[1] <- "Gene"
    return(data)
  })

  volcano_plot <-
    function(dataf, x_name, y_name, slider, color1, color2) {
      if(is.null(dataf))     
        return(NULL)
      dataf<-na.omit(dataf)
      out<- ggplot(dataf,aes(x= !!sym(x_name),y= -log10(!!sym(y_name)),color= !!sym(y_name)<(10^slider)))+ 
        theme_bw()+
        theme(legend.position="bottom")+
        ggtitle('Volcano plot')+
        scale_color_manual(name= paste0("padj < 1 x 10^",toString(slider)), values=c(color1, color2))+
        geom_point()
      return(out)
    }
  
  # Using renderPlot to display the plot
  observeEvent(list(input$x_axis, input$y_axis, input$padjusted, input$de_base, input$de_highlight), {
    output$volcano <- renderPlot( {
      data <- DE_data()
      if (is.null(data)) {
        return()
      }
      isolate({
        volcano_plot(data, input$x_axis, input$y_axis,input$padjusted, input$de_base, input$de_highlight)
      })
    }) 
  })
  
  # Using renderDataTable to display the table
  output$table <- DT::renderDataTable({ 
    data <- DE_data()
    if (is.null(data)) {
      return()
    }
    isolate({
      return(data)
    })
  })
  # ---------------------------------------------------------GSEA Tab------------------------------------------------------------  
  GSEA_load_data <- reactive({
    req(input$GSEA_csvfile)
    file=input$GSEA_csvfile
    if (is.null(file))
    {return(NULL)} 
    else
    {datafile=read_csv(file$datapath)} %>% 
      return()
  })
  
  GSEA_bar <- function(dataf, pathways_slider) {
    filtered <- dataf %>% arrange(padj) %>% slice_head(n=pathways_slider) %>% mutate(status= ifelse(NES > 0, 'positive', 'negative'))
    filtered$pathway <- gsub("\\_", " ", filtered$pathway)
    filtered <- filtered %>% arrange(status)
    pathways <- factor(filtered$pathway)
    filtered$pathway <- factor(filtered$pathway, levels=unique(pathways))
    status <- factor(filtered$status)
    filtered$status <- factor(filtered$status, levels=unique(status))
    plot <- filtered %>% ggplot(aes(x=stringr::str_wrap(pathway, 40), y=NES, fill=status)) +
      geom_col() +
      coord_flip() +
      theme_light() +
      theme(axis.text = element_text(size = 5)) +
      xlab("")
    return(plot)
  }
  
  GSEA_table <- function(dataf, pvalue_choice, type_pathways_choice) {
    data <- dataf %>% mutate(status = ifelse(NES > 0, 'Positive', 'Negative')) %>% filter(padj<=10**pvalue_choice)
    if (type_pathways_choice != 'All') {
      filtered <- data %>% filter(status == type_pathways_choice)
      return(filtered)
    }
    else {
      return(data)
    }
  }
  
  GSEA_scatter <- function(dataf, pvalue_choice) {
    filtered <- dataf %>% mutate(filter_status=ifelse(padj<10**pvalue_choice, "passed filter", "filtered out"))
    plot <- filtered %>% ggplot(aes(x=NES, y=-log10(padj), color=filter_status)) +
      geom_point()
    return(plot)
  }
  
  observeEvent(input$GSEA_Bar, {
    output$GSEA_barplot <- renderPlot({
      dataf <- GSEA_load_data()
      if (is.null(dataf)) {
        return()
      }
      GSEA_bar(dataf = dataf, pathways_slider = input$GSEA_Bar)
    })
  })
  
  observeEvent(c(input$GSEA_table_pvalue, input$GSEA_pathways_choice), {
    output$GSEA_results_table <- renderDataTable({
      dataf <- GSEA_load_data()
      if (is.null(dataf)) {
        return()
      }
      GSEA_table(dataf = dataf,
                 pvalue_choice = input$GSEA_table_pvalue,
                 type_pathways_choice = input$GSEA_pathways_choice)
    })
  })
  
  observeEvent(input$GSEA_scatter_pvalue, {
    output$GSEA_scatterplot <- renderPlot({
      dataf <- GSEA_load_data()
      if (is.null(dataf)) {
        return()
      }
      GSEA_scatter(dataf = dataf,
                   pvalue_choice = input$GSEA_scatter_pvalue)
    })
  })
  
  output$download_NES_table <- downloadHandler(
    filename = function(){"NES_Results.csv"},
    content = function(filename){
      dataf <- GSEA_load_data()
      if (is.null(dataf)) {
        return(NULL)
      }
      write.csv(GSEA_table(dataf = dataf,
                           pvalue_choice = input$GSEA_table_pvalue,
                           type_pathways_choice = input$GSEA_pathways_choice), filename)
    }
  )
  

}

shinyApp(ui = ui, server = server)
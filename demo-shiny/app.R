library(shiny)
library(plotly)
library("heatmaply")
library("shinyHeatmaply")

# demo data to use for plots
source("demo-data.R")
# install.packages(c("heatmaply", "shinyHeatmaply"))

# get user inputs from file
# https://stackoverflow.com/a/39058108/5359531
# https://stackoverflow.com/a/40099556/5359531
# https://shiny.rstudio.com/articles/persistent-data-storage.html
# https://shiny.rstudio.com/reference/shiny/latest/fileInput.html

# Plotly and Shiny
# https://plot.ly/ggplot2/
# https://plot.ly/r/shiny-tutorial/#plotly-graphs-in-shiny
# https://plot.ly/r/shiny-coupled-events/

# Shiny plotly heatmap
# https://www.r-statistics.com/2017/03/shinyheatmaply-a-shiny-app-for-creating-interactive-cluster-heatmaps/
# https://cran.r-project.org/web/packages/shinyHeatmaply/index.html
# https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html

# https://plot.ly/r/shinyapp-plotly-events/

# https://shiny.rstudio.com/articles/modules.html

# https://github.com/talgalili/heatmaply/issues/80
# ~~~~~ UI ~~~~~ #
ui <- shinyUI(fluidPage(
    
    # Shiny
    # plotlyOutput("plot"), # demo plotly mtcars
  
  #####Nick's additions#####
  fixedRow(
    column(6, plotlyOutput("my_pca", height = "600px")),
    column(6, plotlyOutput("my_3dpca", height = "600px"))),
  verbatimTextOutput("id_box"),
  # fixedRow(
  #   column(6, plotlyOutput("my_3dmds", height = "600px")),
  #   column(6, plotlyOutput("my_mds", height = "600px"))),
  # fixedRow(
  #   column(6, plotlyOutput("my_3dtsne", height = "600px")),
  #   column(6, plotlyOutput("my_tsne", height = "600px"))),
  ##########
  
  fixedRow(
        column(6, plotlyOutput("my_heatmap", height = "600px")),
        column(6, plotlyOutput("my_profile_plot", height = "600px"))),
    # plotlyOutput("my_heatmap"),
    verbatimTextOutput("heatmap_hover"),
    verbatimTextOutput("heatmap_selected"),
    verbatimTextOutput("selected_genes"),
    verbatimTextOutput("brush_info"),
    verbatimTextOutput("eventdata"),
    # plotlyOutput("my_profile_plot"),
    
    
    # User UI settings
    selectizeInput(inputId = 'select_gene_ids', label = NULL, choices = NULL, multiple=TRUE),
    
    textInput("control_label",
              "This controls some of the labels:",
              "LABEL TEXT"),
    numericInput("inNumber", "Number input:",
                 min = 1, max = 20, value = 5, step = 0.5),
    radioButtons("inRadio", "Radio buttons:",
                 c("label 1" = "option1",
                   "label 2" = "option2",
                   "label 3" = "option3")),
    
    actionButton("load_inputs", "Load inputs"),
    actionButton('save_inputs', 'Save inputs')
    
))



# ~~~~~ SERVER ~~~~~ #
server <-  shinyServer(function(input, output,session) {
    
    # user selected genes
    select_gene_ids <- character()
    
    # fill the drop down box
    updateSelectizeInput(session = session, inputId = 'select_gene_ids', label = NULL, choices = unique(as.character(wt1[["id"]])), server = TRUE)
    
    # get the user entries
    observeEvent( input$select_gene_ids, {
        if(is.null(input$select_gene_ids)){
            select_gene_ids <- NULL
        } else if (! is.null(input$select_gene_ids)){
            select_gene_ids <- input$select_gene_ids
        }
        }, ignoreNULL=FALSE)
    
    
    #####Nick's additions#####
    # PCA
    output$my_pca <- renderPlotly({
      my_pca %>% layout(dragmode = "select")
    })
    # 3D PCA
    output$my_3dpca <- renderPlotly({

        my_3dpca %>% layout(dragmode = "select")
        
        
    })
    
    # Coupled event-outputting subset ids
    output$id_box <- renderPrint({
      
      cat("Subset Protein IDs\n\nCopy and Paste proteins into the search box\n")
      cat("at https://string-db.org to retrieve network/functional information\n\n")
      
      # Get subset based on selection
      d <- event_data("plotly_selected",source="my_pca")
      
      # If NULL dont do anything
      if(is.null(d) == T) return(NULL)
      
      #Get protein identifiers from subset
      for(i in c(d$pointNumber)){
        name<-rownames(pr_mat$rotation)[i]
        if(grepl(";",name)){
          names<-unlist(strsplit(name,";"))
          for(i in names){
            cat(i,"\n")
          }
        }else{cat(name,"\n")}
      }
      
    })
    # # MDS
    # output$my_mds <- renderPlotly({
    #   my_mds
    # })
    # # 3D MDS
    # output$my_3dmds <- renderPlotly({
    #   my_3dmds %>% layout(dragmode = "select")
    # })
    # # tSNE
    # output$my_tsne <- renderPlotly({
    #   my_tsne
    # })
    # # 3D tSNE
    # output$my_3dtsne <- renderPlotly({
    #   my_3dtsne %>% layout(dragmode = "select")
    # })
    ##########
    
    # Heatmap
    output$my_heatmap <- renderPlotly({
        my_heatmap %>% layout(dragmode = "select")
    })
    
    # Profile Plot
    output$my_profile_plot <- renderPlotly({
        if( ! is.null(input$select_gene_ids) & length(input$select_gene_ids) > 0 & length(input$select_gene_ids) < 50){
            # plot based on input$select_gene_ids
            plot_profile(wt1[which(as.character(wt1[["id"]]) %in% input$select_gene_ids) ,], what = c("id", "expt_id"), color.by = "id", line.smooth = FALSE)
        } else{
            # default plot
            my_profile_plot
        }
    })
    output$heatmap_hover <- renderPrint({
        d <- event_data("plotly_hover","my_heatmap")
        if (is.null(d)) "Hover on a point!" else d
    })
    output$brush_info <- renderPrint({
        cat("input$plot_brush:\n")
        str(input$plot_brush)
    })
    output$heatmap_selected <- renderPrint({
        cat("event_data('plotly_selected'):\n")
        # d <- event_data("plotly_selected")
        # if (is.null(d)) "Select some points!" else d
        event_data("plotly_selected","my_heatmap")
    })
    output$eventdata <- renderPrint({
        cat("str(event_data()):\n")
        str(event_data())
    })
    output$selected_genes <- renderPrint({
        # input$select_gene_ids
        cat("input$select_gene_ids:\n")
        input$select_gene_ids
    })
    
    
    
    # Retrieve saved User UI settings
    observeEvent(input$load_inputs,{
        
        if(!file.exists('inputs.RDS')) {return(NULL)}
        
        savedInputs <- readRDS('inputs.RDS')
        
        inputIDs      <- names(savedInputs)
        inputvalues   <- unlist(savedInputs)
        for (i in 1:length(savedInputs)) {
            session$sendInputMessage(inputIDs[i],  list(value=inputvalues[[i]]) )
        }
    })
    
    observeEvent(input$save_inputs,{
        saveRDS( reactiveValuesToList(input) , file = 'inputs.RDS')
    })
})

shinyApp(ui = ui, server = server)

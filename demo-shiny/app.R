library("shiny")
library("plotly")
library("heatmaply")
library("shinyHeatmaply")

# demo data to use for plots
source("demo-data.R")
source("user_settings_io.R")
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

# https://shiny.rstudio.com/gallery/file-upload.html
# ~~~~~ UI ~~~~~ #
ui <- shinyUI(fluidPage(
  
  # Shiny
  h1("QC & FILTERING"),
  br(), br(),
  h2("% Contamination"),
  br(),
  h2("Peptide counts"),
  br(),
  fixedRow(
  column(6, plotOutput("qc_nPeptides", height = "600px")),
  column(6, plotOutput("qc_rawIntensity", height = "600px"))),
  br(),
  
  ## HEATMAPS
  h1("FIND CO-ELUTING PROTEINS"),
  br(), br(),
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
    rdsFileInput("user_settings_file"), # from the external user_settings_io.R script
    verbatimTextOutput("user_settings"),
    
    
    actionButton("load_inputs", "Load Inputs From Selected File"),

    textInput("control_label",
              "This controls some of the labels:",
              "LABEL TEXT"),
    numericInput("inNumber", "Number input:",
                 min = 1, max = 20, value = 5, step = 0.5),
    radioButtons("inRadio", "Radio buttons:",
                 c("label 1" = "option1",
                   "label 2" = "option2",
                   "label 3" = "option3")),
    
    
    actionButton('save_inputs', 'Save inputs')
    
))



# ~~~~~ SERVER ~~~~~ #
server <-  shinyServer(function(input, output,session) {
    
    # user selected genes
    select_gene_ids <- character()
    
    # fill the drop down box
    updateSelectizeInput(session = session, inputId = 'select_gene_ids', label = NULL, choices = unique(as.character(profile_plot_data[["id"]])), server = TRUE)
    
    # print the user's setting entries
    user_settings_file_print <- callModule(rdsFile_print, "user_settings_file")
    
    output$user_settings <- renderPrint({
        user_settings_file_print()
    })
        
    observeEvent( input$select_gene_ids, {
        if(is.null(input$select_gene_ids)){
            select_gene_ids <- NULL
        } else if (! is.null(input$select_gene_ids)){
            select_gene_ids <- input$select_gene_ids
        }
        }, ignoreNULL=FALSE)
    
    
    
    # QC plots
    output$qc_nPeptides <- renderPlot({p_nPep})
    output$qc_rawIntensity <- renderPlot({p_rawInt})
    
    # Heatmap
    heatmapInput <- reactive({
      key <- data_heatmap$id  ## key identifies brushed subjects
      
      gg1 <- ggplot(data_heatmap, aes(x= fraction, y = order, color = value, key = key)) + 
        geom_point(shape = 15, size = 8) + theme_bw() + 
        theme(legend.position = "none", panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.text.x=element_blank()) +
        scale_colour_gradientn(colours = heat.colors(10)) + ylab("proteins") 
      
      ggplotly(gg1, source = "heatmap") %>% layout(dragmode = "select")
    })
    
    output$my_heatmap <- renderPlotly({
      print(heatmapInput())
    })

    
    
    # Profile Plot
    output$my_profile_plot <- renderPlotly({
      brush <- event_data("plotly_selected", source = "heatmap")
      if( ! is.null(brush)){
        make_profile_plot(df = profile_plot_data, selected_ID = brush$key)
      } else{
        # default plot
        my_profile_plot
      }
    })
    
    
    output$heatmap_hover <- renderPrint({
        d <- event_data("plotly_hover")
        if (is.null(d)) "Hover on a point!" else d
    })
    output$brush_info <- renderPrint({
        cat("input$plot_brush:\n")
        str(input$plot_brush)
    })
    output$heatmap_selected <- renderPrint({
      event_data("plotly_selected", source = "heatmap")
      
    })
    output$eventdata <- renderPrint({
        cat("str(event_data()):\n")
        str(event_data())
    })
    output$selected_genes <- renderPrint({
        cat("input$select_gene_ids:\n")
        input$select_gene_ids
    })
    
    
    # Load the user's UI settings from a Load'ed RDS file
    user_settings_file_load <- callModule(rdsFile_load, "user_settings_file")
    observeEvent(input$load_inputs,{
        savedInputs <- user_settings_file_load()
        inputIDs      <- names(savedInputs)
        inputvalues   <- unlist(savedInputs)
        for (i in 1:length(savedInputs)) {
            session$sendInputMessage(inputIDs[i],  list(value=inputvalues[[i]]) )
        }
    })
    
    # Save the user's current UI settings to a RDS file; path currently hard-coded
    observeEvent(input$save_inputs,{
        saveRDS( reactiveValuesToList(input) , file = 'inputs.RDS')
    })

})

shinyApp(ui = ui, server = server)

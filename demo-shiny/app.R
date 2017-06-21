library(shiny)
library(plotly)
library("heatmaply")
library("shinyHeatmaply")

# demo data to use for plots
source("demo-data2.R")
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
    #updateSelectizeInput(session = session, inputId = 'select_gene_ids', label = NULL, choices = unique(as.character(wt1[["id"]])), server = TRUE)
    updateSelectizeInput(session = session, inputId = 'select_gene_ids', label = NULL, choices = unique(as.character(outcome.df[["subj"]])), server = TRUE)
    
    # get the user entries
    observeEvent( input$select_gene_ids, {
        if(is.null(input$select_gene_ids)){
            select_gene_ids <- NULL
        } else if (! is.null(input$select_gene_ids)){
            select_gene_ids <- input$select_gene_ids
        }
        }, ignoreNULL=FALSE)
    
    
    stuff <- reactive({
      key <- outcome.df$subj  ## key identifies brushed subjects
      
      gg1 <- ggplot(outcome.df, aes(x= grid, y = subj, color = value, key = key)) + 
        geom_point(shape = 15, size = 8) + theme_bw() + 
        theme(legend.position = "none", panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background=element_blank()) +
        scale_colour_gradientn(colours = heat.colors(10))
      
      
        
      gg1_ok <- ggplot(outcome.df, aes(x = grid, y = value, key = key)) + geom_point()
      
      ggplotly(gg1, source = "heatM") %>% layout(dragmode = "select")
    })
    
    
    # Heatmap
    output$my_heatmap <- renderPlotly({
        print(stuff())
        #my_heatmap  %>% layout(dragmode = "select")
    })
    

    
    stuff2 <- reactive({
      
      brush <- event_data("plotly_selected", source = "heatM")
      baseplot = ggplot(outcome.df, aes(grid, value, group = subj)) + theme_bw() + 
        geom_point() + geom_line()
      
      if (is.null(brush)) {
        brushed_subjs = NULL
        baseplot.gg =  baseplot
      }else{
        brushed_subjs = brush$key
        baseplot.gg = outcome.df %>% filter(subj %in% brushed_subjs) %>%
          ggplot(., aes(grid, value, group = subj)) + geom_line(color = "green")
      }
      ggplotly(baseplot.gg)
    })
    
    
    # Profile Plot
    output$my_profile_plot <- renderPlotly({
      #if( ! is.null(input$select_gene_ids) & length(input$select_gene_ids) > 0 & length(input$select_gene_ids) < 50){
      # plot based on input$select_gene_ids
      #   plot_profile(wt1[which(as.character(wt1[["id"]]) %in% input$select_gene_ids) ,], what = c("id", "expt_id"), color.by = "id", line.smooth = FALSE)
      #} else{
      # default plot
      #   my_profile_plot 
      #}
      #ggplot(outcome.df, aes(grid, value, group = subj)) + theme_bw() + 
      #geom_point() + geom_line()
      print(stuff2())
    })
    
    


    output$heatmap_selected <- renderPrint({
        #cat("event_data('plotly_selected'):\n")
        # d <- event_data("plotly_selected")
        # if (is.null(d)) "Select some points!" else d
        #event_data("plotly_selected")
        #brush <- event_data("plotly_selected", source = "heatM")
        #brush$key
       event_data("plotly_selected", source = "heatM")
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
    
    

    observeEvent(input$save_inputs,{
        saveRDS( reactiveValuesToList(input) , file = 'inputs.RDS')
    })
})

shinyApp(ui = ui, server = server)

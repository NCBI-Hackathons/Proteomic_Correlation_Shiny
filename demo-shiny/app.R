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

# ~~~~~ UI ~~~~~ # 
ui <- shinyUI(fluidPage(
    
    # Shiny
    # plotlyOutput("plot"), # demo plotly mtcars
    plotlyOutput("my_profile_plot"),
    verbatimTextOutput("event"),
    plotlyOutput("my_heatmap"),

    # User UI settings
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
    
    # Plotly
    # renderPlotly() also understands ggplot2 objects!
    # output$plot <- renderPlotly({
    #     plot_ly(mtcars, x = ~mpg, y = ~wt)
    # }) # demo plotly mtcars
    output$my_profile_plot <- renderPlotly({
        my_profile_plot
    })
    output$my_heatmap <- renderPlotly({
        my_heatmap
    })
    my_heatmap
    
    output$event <- renderPrint({
        d <- event_data("plotly_hover")
        if (is.null(d)) "Hover on a point!" else d
    })
    
    # User UI settings
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

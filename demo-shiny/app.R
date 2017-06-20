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
    verbatimTextOutput("eventdata"),
    # plotlyOutput("my_profile_plot"),

    # Some help text from the coupled events
    h2("Coupled events in plotly charts using Shiny"),
    h4("This Shiny app showcases coupled events using Plotly's ", tags$code("event_data()"), " function."),
    tags$ol(
        tags$li("The first chart showcases", tags$code("plotly_selected")),
        tags$li("The third chart showcases", tags$code("plotly_click"))
    ),


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


    # Profile Plot
    output$my_profile_plot <- renderPlotly({
        my_profile_plot
    })

    # Heatmap
    output$my_heatmap <- renderPlotly({
        # Get subset based on selection
        # event.data <- event_data("plotly_selected", source = "subset")
        # If NULL dont do anything
        # if(is.null(event.data) == T) return(NULL)

        # Get subset of original data from selection
        # malig.class <- subset(plot.df, Class == "malignant")[subset(event.data, curveNumber == 0)$pointNumber + 1,]
        # benign.class <- subset(plot.df, Class == "benign")[subset(event.data, curveNumber == 1)$pointNumber + 1,]

        # Plot object call
        # my_heatmap # wt1_mat_subset
        
        my_heatmap#  %>% layout(dragmode = "select")
        # iris_heatmap %>% layout(dragmode = "select")

        })
    

    output$heatmap_hover <- renderPrint({
        d <- event_data("plotly_hover")
        if (is.null(d)) "Hover on a point!" else d
    })
    output$heatmap_selected <- renderPrint({
        # d <- event_data("plotly_selected")
        # if (is.null(d)) "Select some points!" else d
        event_data("plotly_selected")
    })
    output$eventdata <- renderPrint({
        str(event_data())
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

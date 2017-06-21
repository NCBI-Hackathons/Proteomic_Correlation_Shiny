# functions to set up the user settings file input and output
# https://shiny.rstudio.com/articles/modules.html
# https://stackoverflow.com/a/40099556/5359531
# https://stackoverflow.com/a/39058108/5359531

# Module UI function
rdsFileInput <- function(id, label = "Choose Settings RDS File (Press 'Load Inputs...' after upload completes") {
    # Create a namespace function using the provided id
    ns <- NS(id)
    
    tagList(
        fileInput(ns("file"), label, accept=c('.RDS'))
    )
}

# Module server function
rdsFile_print <- function(input, output, session) {
    # Print the selected RDS file, if any
    userFile <- reactive({
        # If no file is selected, don't do anything
        validate(need(input$file, message = FALSE))
        input$file
    })
    
    # The user's filepath
    filepath <- reactive({
        cat(sprintf("Magical path to your file: %s", userFile()$datapath))
    })
    
    # We can run observers in here if we want to
    observe({
        msg <- sprintf("File %s was uploaded", userFile()$name)
        cat(msg, "\n")
    })
    
    # Return the reactive that yields the path
    return(filepath)
}

rdsFile_load <- function(input, output, session) {
    # Retrieve saved User UI settings from the RDS file
    userFile <- reactive({
        # If no file is selected, don't do anything
        validate(need(input$file, message = FALSE))
        input$file
    })
    
    savedInputs <- reactive({        
        if(!file.exists(userFile()$datapath)) {return(NULL)}
        readRDS(userFile()$datapath)
    })
    return(savedInputs)
}


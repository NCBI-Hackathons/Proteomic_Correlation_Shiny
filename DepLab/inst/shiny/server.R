
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggplot2)
library(scales)
library(colorspace)
library(scales)
library(shinyBS)
library(shinyFiles)
library(stringr)
library(shinythemes)
library(RSQLite)
library(stringi)
library(NMF)
library(shinyjs)

# load code to be run once (at app launch)


options(stringsAsFactors = FALSE)

benschop_complexes <- DepLab:::read_complexes(filename = system.file("extdata", "defined_complexes", "complexes_benschop.txt", package = "DepLab"),
                                              organism = "yeast")
wodak_complexes <- DepLab:::read_complexes(filename = system.file("extdata", "defined_complexes", "complexes_wodak.txt", package = "DepLab"),
                                           organism = "yeast")
corum_complexes <- DepLab:::read_complexes(filename = system.file("extdata", "defined_complexes", "human_complexes_corum.tab", package = "DepLab"),
                                           organism = "human")


if (file.size(system.file("extdata","path_to_db.txt",package = "DepLab")) == 0){
  if( .Platform$OS.type == "windows"){
    platform_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop/proteomics.db")
  } else {
    platform_path <- "~/Desktop/proteomics.db"
  }
  write.table(platform_path, system.file("extdata","path_to_db.txt",package = "DepLab"), col.names=F, row.names=F) 
}

path_to_db <- read.table(file = system.file("extdata","path_to_db.txt",
                                            package = "DepLab"),
                         sep="\t", stringsAsFactors=FALSE,
                         header=FALSE, quote="\"")

if (file.size(system.file("extdata","path_to_custom_complexes.txt",package = "DepLab")) == 0){
  if( .Platform$OS.type == "windows"){
    platform_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop/complexes_custom.txt")
  } else {
    platform_path <- "~/Desktop/complexes_custom.txt"
  }
  write.table(platform_path, system.file("extdata","path_to_custom_complexes.txt",package = "DepLab"), col.names=F, row.names=F) 
}

path_to_custom_complexes <- read.table(file = system.file("extdata","path_to_custom_complexes.txt",
                                                          package = "DepLab"),
                                       sep="\t", stringsAsFactors=FALSE,
                                       header=FALSE, quote="\"")


custom_complexes_name <- normalizePath(path_to_custom_complexes$V1, winslash = .Platform$file.sep)

database.name <- normalizePath(path_to_db$V1, winslash = .Platform$file.sep)

if(file.info(database.name)$size == 0 || is.na(file.info(database.name)$size)) {
  cur_frac_table_length <- 0
} else {
  cur_frac_table_length <- as.numeric(dbGetQuery(conn = dbConnect(SQLite(), dbname = database.name, cache_size = 5000), "SELECT MAX(_ROWID_) FROM frac_data LIMIT 1;")[1])
}

#database.name <- system.file("extdata", "proteomics.db", package = "DepLab")

radio_list <- list()

if(file.info(custom_complexes_name)$size == 0 || is.na(file.info(custom_complexes_name)$size)) {
  newDF <- data.frame("id"= "id", "gene symbol"="gene_symbol", "complex"="complex")
  write.table(newDF, custom_complexes_name, append = FALSE, row.names = FALSE, col.names = FALSE, na = "NA", quote=F, sep="\t")
} 

if(file.info(database.name)$size == 0 || is.na(file.info(database.name)$size)) {
  initialize.database(database.name, organism = "human", force = FALSE)
  initialize.database(database.name, organism = "yeast", force = FALSE)
} 

shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2) 
  # load code to be run when a page is visited
  volumes <- getVolumes() #c('R Installation'=R.home())
  shinyFileSave(input, 'save', roots=volumes, session=session)
  saveFileName <- renderPrint({parseSavePath(volumes, input$save)})
  observeEvent( input$save, {
    savePath <- strsplit(saveFileName(), " pdf | png ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][2]
    ggsave(savePath)
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("plot", "pdf", sep = ".")
    },
    content = function(file) {
      ggsave(filename = file, device = "pdf", width=10, height=10)
    }
  )
  
  output$downloadPlot2 <- downloadHandler(
    filename = function() {
      paste("plot", "pdf", sep = ".")
    },
    content = function(file) {
      ggsave(filename = file, device = "pdf", width=10, height=10)
    }
  )
  
  output$downloadPlot3 <- downloadHandler(
    filename = function() {
      paste("plot", "pdf", sep = ".")
    },
    content = function(file) {
      ggsave(filename = file, device = "pdf", width=10, height=10)
    }
  )
  
  
  output$downloadPlot_heatmap <- downloadHandler(
    filename = function() {
      paste("plot", "pdf", sep = ".")
    },
    content = function(file) {
      ggsave(filename = file, device = "pdf", width=10, height=10)
    }
  )
  
  
  shinyFileChoose(input,  'db_path' , session=session, roots=volumes, filetypes=c('db' ))
  dbPathName <- renderText({  unlist(parseFilePaths(volumes, input$db_path))   })   #$datapath[1])
  
  database.name.reactive <- reactiveValues(data = database.name)
  
  observeEvent( input$db_path, {
    database.name.reactive$data <- gsub(".*  (.+)","\\1",dbPathName())
    write.table( normalizePath(gsub(".*  (.+)","\\1",dbPathName()), winslash=.Platform$file.sep) , system.file("extdata","path_to_db.txt",package = "DepLab"), row.names=F, col.names=F)
    updateSelectizeInput(session, 'show_expt_id', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
    updateSelectizeInput(session, 'show_expt_id_std', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
    updateSelectizeInput(session, 'show_trypsin_symbol', choices = DepLab:::list.std.gene.symbols(database.name.reactive$data)$id, server = TRUE)
    updateSelectizeInput(session, 'show_expt_id_sum', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
    updateSelectizeInput(session, 'show_expt_id_db_browser', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
    updateSelectizeInput(session, 'show_expt_id_db_viewer', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
    updateSelectizeInput(session, 'show_expt_id_heatmap0', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
    for (i in heatmap.reactive$heatmap_selectize_lists){
      updateSelectizeInput(session, i, choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
    }
    
  })
  
  
  output$db_dir_path <- renderText({ 
    paste("Currently selected database:", database.name.reactive$data, sep = "\n")
  })
  
  output$expt_id_list <- renderText({ 
    expt_id_msg <- paste("The currently selected database contains the following data sets from previous sessions:", paste(as.data.frame(DepLab:::list.expt.ids.w.organism(database.name.reactive$data))$expt_id, collapse="\n"), sep = "\n")
  })
  
  
  shinyFileChoose(input,  'custom_complexes_path' , session=session, roots=volumes, filetypes=c('txt' ))
  customComplexesPathName <- renderText({  unlist(parseFilePaths(volumes, input$custom_complexes_path))   })   #$datapath[1])
  
  custom.complexes.name.reactive <- reactiveValues(data = custom_complexes_name)
  
  observeEvent( input$custom_complexes_path, {
    custom.complexes.name.reactive$data <- gsub(".*  (.+)","\\1",customComplexesPathName())
    write.table( normalizePath(gsub(".*  (.+)","\\1",customComplexesPathName()), winslash = .Platform$file.sep) , system.file("extdata","path_to_custom_complexes.txt",package = "DepLab"), row.names=F, col.names=F)
    custom_complexes <- reactiveValues(data = read.table(file=custom.complexes.name.reactive$data, sep = "\t", stringsAsFactors = FALSE, header = TRUE))
    if (input$complex_source == "custom") {
      updateSelectizeInput(session, 'show_complex', choices = c(unique(custom_complexes$data$complex)))
    } 
    #print(custom_complexes$data)
  })
  
  output$custom_complexes_dir_path <- renderText({ 
    paste("The path to store custom complexes is set to:", custom.complexes.name.reactive$data, sep = "\n")
  })
  
  
  #updateSelectizeInput(session, 'show_complex', choices = c(unique(benschop_complexes$data$complex)))
  
  observeEvent( input$y_axis_choices, {
    if (input$y_axis_choices == "raw.intensity"){
      updateCheckboxGroupInput(session, "y_axis_log", choices = c("log2", "normalize across fractions", "normalize by spike-in","smooth values (SuperSmoother)", "smooth lines (spline)"))
    } else {
      updateCheckboxGroupInput(session, "y_axis_log",  choices = c("log2"))
    }
  })
  
  custom_complexes <- reactiveValues(data = read.table(file=custom_complexes_name, sep = "\t", stringsAsFactors = FALSE, header = TRUE))
  
  fileName<-reactive({
    inFile <- input$chrom_image
    inFile$datapath
  })
  
  
  vars = reactiveValues(protein_counter = 0, pro_info_table = NULL, complex_info_table = NULL, complex_info_table_tmp = NULL)
  
  output$retrieve_protein_button <- renderUI({
    if(!is.null(input$show_expt_id)){
      if (is.null(input$show_gene_symbol) & is.null(input$show_complex)){
        return()
      } else {
        actionLink("retrieve_protein_info", label=retrieve_protein_button_label())
      }
    }
  })
  
  
  rep_list <- list()
  
  ### code to specify reps on individ protein page
  observeEvent(c(input$show_expt_id, input$specify_replicates), {
    if(is.null(input$show_expt_id)){
      rep_list <- NULL
      #   updateSelectizeInput(session, "split_by",  choices = list("None" = "none", "Gene symbol" = "gene_symbol", "ID" = "id", "Expt ID" = "expt_id", "Complex" = "complex"), selected = "expt_id")
      #  updateSelectizeInput(session, "split_by_col",  choices = list("None" = "none", "Gene symbol" = "gene_symbol", "ID" = "id", "Expt ID" = "expt_id", "Complex" = "complex"), selected = "expt_id")
    } else if (input$specify_replicates == FALSE){
      rep_list <- NULL
      #  updateSelectizeInput(session, "split_by",  choices = list("None" = "none", "Gene symbol" = "gene_symbol", "ID" = "id", "Expt ID" = "expt_id", "Complex" = "complex"), selected = "expt_id")
      #  updateSelectizeInput(session, "split_by_col",  choices = list("None" = "none", "Gene symbol" = "gene_symbol", "ID" = "id", "Expt ID" = "expt_id", "Complex" = "complex"), selected = "expt_id")
    } else {
      #  updateSelectizeInput(session, "split_by",  choices = list("None" = "none", "Gene symbol" = "gene_symbol", "ID" = "id", "Expt ID" = "expt_id", "Complex" = "complex", "Replicate" = "replicate"), selected = "expt_id")
      # updateSelectizeInput(session, "split_by_col",  choices = list("None" = "none", "Gene symbol" = "gene_symbol", "ID" = "id", "Expt ID" = "expt_id", "Complex" = "complex",  "Replicate" = "replicate"), selected = "expt_id")
      len_expt <- length(input$show_expt_id)
      for (i in 1:len_expt){
        rep_nm <- paste0("radio_var", i)
        rep_nm_cond <- paste0("rep_nm_cond", i)
        rep_nm_rep <- paste0("rep_nm_rep", i)
        if( i == 1){
          rep_list[[rep_nm]] <-  fluidRow(column(6, offset = 0, style='margin-bottom:0px;', textInput(rep_nm_cond, HTML(paste("Cond.*",paste("<span style='font-weight:normal'>", input$show_expt_id[i], "</span>", sep=""), sep="<br/>")),  value=NULL)),
                                          column(4,offset = 0, style='margin-bottom:0px;margin-top:0px; ', numericInput(rep_nm_rep, label= HTML(paste("Rep.*","<br/>", "<br /> ")), min = 1, value=NULL)))     
        } else {
          rep_list[[rep_nm]] <-  fluidRow(column(6, offset = 0, style='margin-top:0px;margin-bottom:0px;', textInput(rep_nm_cond, label = HTML(paste("<span style='font-weight:normal'>", input$show_expt_id[i], "</span>", sep="")) , value=NULL)),
                                          column(4,offset = 0, style='margin-top:5px;margin-bottom:0px;', numericInput(rep_nm_rep, min = 1, label = '', value=NULL)))     
        }
      }
    }
    output$uiRepList <- renderUI({
      tagList(list(rep_list))
    })
  }, ignoreNULL=FALSE)
  
  ###
  
  
  
  observeEvent( input$show_gene_symbol, {
    if(is.null(input$show_gene_symbol) & is.null(input$show_complex)){
      vars$protein_counter <-  0
    }
  }, ignoreNULL=FALSE)
  
  
  observeEvent(c(input$show_expt_id, input$y_axis_log), {
    if(is.null(input$show_expt_id)){
      radio_list <- NULL
    } else if (!any(grepl("spike", c(input$y_axis_log)))){
      radio_list <- NULL
    } else {
      len_expt <- length(input$show_expt_id)
      for (i in 1:len_expt){
        rad_nm <- paste0("radio_var", i)
        radio_list[[rad_nm]] <- radioButtons(rad_nm, input$show_expt_id[i], choices = (DepLab:::list.selected.std.gene.symbols(database.name.reactive$data, input$show_expt_id[i])$id), inline=TRUE)   
      }
    }
    output$uiRadio <- renderUI({
      tagList(list(radio_list))
    })
  }, ignoreNULL=FALSE)
  
  
  observeEvent(input$retrieve_protein_info, {
    if(!is.null(input$retrieve_protein_info)){
      input$retrieve_protein_info
      isolate({
        vars$protein_counter <-  vars$protein_counter + 1
      })
    }
  })
  
  ## timer
  reac <- reactiveValues(redraw = TRUE, show_gene_symbol = isolate(input$show_gene_symbol), show_expt_id =  isolate(input$show_expt_id), show_complex =  isolate(input$show_complex))
  
  observe({
    input$show_gene_symbol
    input$show_expt_id
    input$show_complex
    reac$redraw <- FALSE
  })
  
  observe({
    invalidateLater(1000, session)
    input$show_gene_symbol
    input$show_expt_id
    input$show_complex
    # input$redraw
    #  isolate(cat(reac$redraw, input$show_gene_symbol, "\n"))
    if (isolate(reac$redraw)) {
      reac$show_gene_symbol <- input$show_gene_symbol
      reac$show_expt_id <- input$show_expt_id
      reac$show_complex <- input$show_complex
    } else {
      isolate(reac$redraw <- TRUE)
    }
  })
  
  ###
  
  ### timer for sum page
  reacSum <- reactiveValues(redrawSum = TRUE, show_expt_id_sum = isolate(input$show_expt_id_sum))
  
  observe({
    input$show_expt_id_sum
    reacSum$redrawSum <- FALSE
  })
  
  observe({
    invalidateLater(1000, session)
    input$show_expt_id_sum
    if (isolate(reacSum$redrawSum)) {
      reacSum$show_expt_id_sum <- input$show_expt_id_sum
    } else {
      isolate(reacSum$redrawSum <- TRUE)
    }
  })
  
  ###
  
  retrieve_protein_button_label <- reactive({
    if(vars$protein_counter >= 1) label <- "Refresh UniProt info..."
    else label <- "Retrieve info from UniProt..."
  })
  
  
  observe({
    if(is.null(fileName())) return(NULL)
    output$preImage <- renderImage({
      filename <- fileName()
      # Return a list containing the filename and alt text
      list(src = filename,
           contentType = 'image/png',
           width = 400,
           height = 300,
           alt = "chromatogram")
      
    }, deleteFile = FALSE)
  })
  
  # ===================
  # Data input tab:
  
  ### dynamic UI for Data input tab
  createAlert(session, "inputalert", "inputAlert", title = "Oops", content = "Experimental ID, Maxquant data, and sample origin meta-data required before saving.", append = FALSE)
  
  mandatory.data.input <- reactiveValues(fieldsMandatory = c("expt.id", "file1", "experimenter.name", "genotype", "cell.type", "buffer.composition", "lysis.method", "digestion.enzyme"))
  observe({
    input$include.prefractionation.metadata 
    input$include.msmethod.metadata
    input$include.data.metadata
    if (input$include.prefractionation.metadata == TRUE){
      shinyjs::show("column.id")
      shinyjs::show("amount.protein.loaded")
      shinyjs::show("sample.vol.loaded")
      shinyjs::show("lc.flow.rate")
      shinyjs::show("lc.fraction.size")
      shinyjs::show("time.per.fraction")
      shinyjs::show("fractions.collected")
      #      shinyjs::show("chrom_image")
    } else {
      shinyjs::hide("column.id")
      shinyjs::hide("amount.protein.loaded")
      shinyjs::hide("sample.vol.loaded")
      shinyjs::hide("lc.flow.rate")
      shinyjs::hide("lc.fraction.size")
      shinyjs::hide("time.per.fraction")
      shinyjs::hide("fractions.collected")
      #      shinyjs::hide("chrom_image")
    }
    
    if (input$include.msmethod.metadata == TRUE){
      shinyjs::show("instrument.id")
      shinyjs::show("run.date")
      shinyjs::show("method.length")
    } else {
      shinyjs::hide("instrument.id")
      shinyjs::hide("run.date")
      shinyjs::hide("method.length")
    }
    
    if (input$include.data.metadata == TRUE){
      shinyjs::show("processing.platform")
      shinyjs::show("search.algorithm")
      shinyjs::show("filtering.algorithm")
      shinyjs::show("filtering.stringency")
    } else {
      shinyjs::hide("processing.platform")
      shinyjs::hide("search.algorithm")
      shinyjs::hide("filtering.algorithm")
      shinyjs::hide("filtering.stringency")
    }
    # check if all mandatory fields have a value
    # ALL TRUE
    if ((input$include.prefractionation.metadata == TRUE) &( input$include.msmethod.metadata == TRUE) & (input$include.data.metadata == TRUE)){
      mandatory.data.input$fieldsMandatory <- c("expt.id", "file1", "experimenter.name", "genotype", "cell.type", "buffer.composition", "lysis.method", "digestion.enzyme",
                                                "column.id","amount.protein.loaded","sample.vol.loaded","lc.flow.rate","lc.fraction.size","time.per.fraction","fractions.collected",
                                                "instrument.id", "method.length", "processing.platform", "search.algorithm", "filtering.algorithm", "filtering.stringency")	
      closeAlert(session, "inputAlert")
      createAlert(session, "inputalert", "inputAlert", title = "Oops", content = "Experimental ID, Maxquant data, sample origin meta-data, prefractionation methods, mass-spec methods, and data processing methods required before saving.", append = FALSE)
    } 
    
    
    if ((input$include.prefractionation.metadata == TRUE) & (input$include.msmethod.metadata == TRUE) & (input$include.data.metadata == FALSE)){
      mandatory.data.input$fieldsMandatory <- c("expt.id", "file1", "experimenter.name", "genotype", "cell.type", "buffer.composition", "lysis.method", "digestion.enzyme",
                                                "column.id","amount.protein.loaded","sample.vol.loaded","lc.flow.rate","lc.fraction.size","time.per.fraction","fractions.collected",
                                                "instrument.id",  "method.length")	
      closeAlert(session, "inputAlert")
      createAlert(session, "inputalert", "inputAlert", title = "Oops", content = "Experimental ID, Maxquant data, sample origin meta-data, prefractionation methods, and mass-spec methods required before saving.", append = FALSE)
    } 
    
    
    
    if ((input$include.prefractionation.metadata == TRUE) & (input$include.msmethod.metadata == FALSE) & (input$include.data.metadata == TRUE)){
      mandatory.data.input$fieldsMandatory <- c("expt.id", "file1", "experimenter.name", "genotype", "cell.type", "buffer.composition", "lysis.method", "digestion.enzyme",
                                                "column.id","amount.protein.loaded","sample.vol.loaded","lc.flow.rate","lc.fraction.size","time.per.fraction","fractions.collected",
                                                "processing.platform", "search.algorithm", "filtering.algorithm", "filtering.stringency")	
      closeAlert(session, "inputAlert")
      createAlert(session, "inputalert", "inputAlert", title = "Oops", content = "Experimental ID, Maxquant data, sample origin meta-data, prefractionation methods, and data processing methods required before saving.", append = FALSE)
      
    } 
    
    
    if ((input$include.prefractionation.metadata == FALSE) & (input$include.msmethod.metadata == TRUE) & (input$include.data.metadata == TRUE)){
      mandatory.data.input$fieldsMandatory <- c("expt.id", "file1", "experimenter.name", "genotype", "cell.type", "buffer.composition", "lysis.method", "digestion.enzyme",
                                                "instrument.id",  "method.length",
                                                "processing.platform", "search.algorithm", "filtering.algorithm", "filtering.stringency")	
      closeAlert(session, "inputAlert")
      createAlert(session, "inputalert", "inputAlert", title = "Oops", content = "Experimental ID, Maxquant data, sample origin meta-data, mass-spec methods, and data processing methods required before saving.", append = FALSE)
      
    } 
    
    if ((input$include.prefractionation.metadata == TRUE) & (input$include.msmethod.metadata == FALSE) & (input$include.data.metadata == FALSE)){
      mandatory.data.input$fieldsMandatory <- c("expt.id", "file1", "experimenter.name", "genotype", "cell.type", "buffer.composition", "lysis.method", "digestion.enzyme",
                                                "column.id","amount.protein.loaded","sample.vol.loaded","lc.flow.rate","lc.fraction.size","time.per.fraction","fractions.collected")	
      closeAlert(session, "inputAlert")
      createAlert(session, "inputalert", "inputAlert", title = "Oops", content = "Experimental ID, Maxquant data, sample origin meta-data, and prefractionation methods required before saving.", append = FALSE)
    } 
    
    
    
    if ((input$include.prefractionation.metadata == FALSE )& (input$include.msmethod.metadata == TRUE) & (input$include.data.metadata == FALSE)){
      mandatory.data.input$fieldsMandatory <- c("expt.id", "file1", "experimenter.name", "genotype", "cell.type", "buffer.composition", "lysis.method", "digestion.enzyme",
                                                "instrument.id",  "method.length")
      closeAlert(session, "inputAlert")
      createAlert(session, "inputalert", "inputAlert", title = "Oops", content = "Experimental ID, Maxquant data, sample origin meta-data, and mass-spec methods required before saving.", append = FALSE)
      
    } 
    
    
    if (input$include.prefractionation.metadata == FALSE & input$include.msmethod.metadata == FALSE & input$include.data.metadata == TRUE){
      mandatory.data.input$fieldsMandatory <- c("expt.id", "file1", "experimenter.name", "genotype", "cell.type", "buffer.composition", "lysis.method", "digestion.enzyme",
                                                "processing.platform", "search.algorithm", "filtering.algorithm", "filtering.stringency")	
      closeAlert(session, "inputAlert")
      createAlert(session, "inputalert", "inputAlert", title = "Oops", content = "Experimental ID, Maxquant data, sample origin meta-data, and data processing methods required before saving.", append = FALSE)
      
    }
    
    #ALL FALSE
    if ((input$include.prefractionation.metadata == FALSE) & (input$include.msmethod.metadata == FALSE) & (input$include.data.metadata == FALSE)){
      mandatory.data.input$fieldsMandatory <- c("expt.id", "file1", "experimenter.name", "genotype", "cell.type", "buffer.composition", "lysis.method", "digestion.enzyme")	
      closeAlert(session, "inputAlert")
      createAlert(session, "inputalert", "inputAlert", title = "Oops", content = "Experimental ID, Maxquant data, and sample origin meta-data required before saving.", append = FALSE)
      
    } 
    
    mandatoryFilledDataInput <-
      vapply(mandatory.data.input$fieldsMandatory,
             function(x) {
               !is.null(input[[x]]) && input[[x]] != ""
             },
             logical(1))
    mandatoryFilledDataInput <- all(mandatoryFilledDataInput)
    # enable/disable the submit button
    shinyjs::toggleState(id = "saveButton", condition = mandatoryFilledDataInput)
    if (mandatoryFilledDataInput){ 
      closeAlert(session, "inputAlert")
      return()
    }
  })  
  
  
  
  # makes sure the currently selected variables are used to populate the vbnc data frame
  datasetInput <- reactive({
    num.replicates <- input$replicate.number
    if (input$stress.deviation == "None") {
      stress.deviation <- list(stress.deviation = "None",
                               stress.deviation.value = NA,
                               stress.deviation.unit = NA)
    } else {
      stress.deviation <- list(stress.deviation = input$stress.deviation,
                               stress.deviation.value = input$stress.deviation.value,
                               stress.deviation.unit = input$stress.deviation.unit)
    }
    
    if (input$resuscitation.media == "None") {
      resuscitation <- list(resuscitation.media = NA,
                            resuscitation.media.supplement = NA,
                            resuscitation.media.treatment = NA,
                            resuscitation.media.treatment.value = NA,
                            resuscitation.media.treatment.unit = NA)
    } else {
      if (input$resuscitation.media.supplement == "None") {
        resuscitation <- list(resuscitation.media = input$resuscitation.media,
                              resuscitation.media.supplement = NA,
                              resuscitation.media.treatment = NA,
                              resuscitation.media.treatment.value = NA,
                              resuscitation.media.treatment.unit = NA)
      } else {
        if (input$resuscitation.media.treatment == "None") {
          resuscitation <- list(resuscitation.media = input$resuscitation.media,
                                resuscitation.media.supplement = input$resuscitation.media.supplement,
                                resuscitation.media.treatment = NA,
                                resuscitation.media.treatment.value = NA,
                                resuscitation.media.treatment.unit = NA)
        } else {
          resuscitation <- list(resuscitation.media = input$resuscitation.media,
                                resuscitation.media.supplement = input$resuscitation.media.supplement,
                                resuscitation.media.treatment = input$resuscitation.media.treatment,
                                resuscitation.media.treatment.value = input$resuscitation.media.treatment.value,
                                resuscitation.media.treatment.unit = input$resuscitation.media.treatment.unit)
        }
      }
    }
    
    new.row <- c(read.time = input$read.time,
                 culture.strain = input$culture.strain,
                 stress.media = input$stress.media,
                 stress.duration = input$stress.duration,
                 stress.deviation,
                 specimen.vessel = input$specimen.vessel,
                 treatment = input$treatment,
                 treatment.value = input$treatment.value,
                 treatment.unit = input$treatment.unit,
                 treatment.exposure.time = input$treatment.exposure.time,
                 times.washed = input$times.washed,
                 vol.washing.buffer = input$vol.washing.buffer,
                 resuscitation)
    
    for (i in 1:num.replicates) {
      expt.id <- paste(input$stress.media,
                       input$treatment,
                       paste(input$stress.duration, "d", sep=""),
                       input$resuscitation.media,
                       LETTERS[i], sep="_")
      row <- c(expt.id = expt.id, replicate.number = LETTERS[i], new.row)
      vbnc <- rbind(vbnc, row)
    }
    return(vbnc)
  })
  
  output$newdata <- renderTable({
    dataset <- datasetInput()
    dataset
  })
  
  
  output$helptextsave <- renderUI({
    if(is.null(input$expt.id)){
      helpText('Select expt ID')
    }
  })
  
  ### dynamic UI for heatmap tab
  
  heatmap.reactive <- reactiveValues(data = NULL, fieldsMandatory = c("heatmap.cond0", "heatmap.rep0", "show_expt_id_heatmap0"), heatmap_selectize_lists = "show_expt_id_heatmap0", pairwise_index = NULL, conditionwise_index = NULL)
  heatmap_inserted <- character(0)
  
  observeEvent(input$heatmap.add, {
    btn <- input$heatmap.add
    id <- paste0('txt', btn)
    insertUI(
      selector = '#placeholderHeatmapButtons',
      ui = tags$div(
        
        fluidRow(column(2, offset = 0, style='padding:0px; margin-top:-20px;', textInput(paste0("heatmap.cond", input$heatmap.add), label = "",  value=NULL ))   , #  value = paste0("cond", input$heatmap.add) 
                 column(2, offset = 0, style='padding:0px; margin-top:-20px;',  numericInput(paste0("heatmap.rep", input$heatmap.add), min=1, label = "", value = NULL)),
                 column(6, offset = 0, style='padding:0px; margin-top:-20px;',  selectizeInput(paste0("show_expt_id_heatmap", input$heatmap.add), "", choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, selected=NULL, multiple=TRUE,      options = list(
                   maxItems = 1)))),
        id = id
      )
    )
    heatmap_inserted <<- c(heatmap_inserted, id)
    
    len <- c(gsub("txt", "", heatmap_inserted), 0)
    exps <- paste0("show_expt_id_heatmap", len)
    replicate <- paste0("heatmap.rep",len)
    condition <- paste0("heatmap.cond", len)
    heatmap.reactive$fieldsMandatory <- c(exps, replicate, condition)
    heatmap.reactive$heatmap_selectize_lists <- c(exps)
  })
  
  
  observeEvent(input$heatmap.delete, {
    removeUI(
      ## pass in appropriate div id
      selector = paste0('#', heatmap_inserted[length(heatmap_inserted)])
    )
    heatmap_inserted <<- heatmap_inserted[1:length(heatmap_inserted)-1]
    len <- c(gsub("txt", "", heatmap_inserted), 0)
    exps <- paste0("show_expt_id_heatmap", len)
    replicate <- paste0("heatmap.rep",len)
    condition <- paste0("heatmap.cond", len)
    heatmap.reactive$fieldsMandatory <- c(exps, replicate, condition)
    heatmap.reactive$heatmap_selectize_lists <- c(exps)
  })
  
  hmres <- list()
  hmtitle.pw <- "Correlations between individual samples"
  hmtitle.cond <-  "Correlations between the median intensities of the conditions"
  
  ## Pairwise heatmaps
  output$heatmap_pairwise_all <- renderPlot({
    heatmap.reactive$data
    input$heatmap.pairwise.corr.range
    if(!is.null(heatmap.reactive$data$pairwise)){
      plotdat <- heatmap.reactive$data$pairwise
      plotdat_col_names <- colnames(plotdat)
      min_row <- apply(plotdat,1,min)
      max_row <- apply(plotdat,1,max)
      plotdat <-  as.data.frame(plotdat[(min_row >= input$heatmap.pairwise.corr.range[1]) & (max_row <= input$heatmap.pairwise.corr.range[2]),])
      colnames(plotdat) <- plotdat_col_names
      # pairwise heatmap
      if(dim(plotdat)[2] == 1){
        hmres$pairwise <- aheatmap(plotdat, Colv = NA, Rowv = NA, breaks = 0, main = hmtitle.pw)
        heatmap.reactive$pairwise_index <- hmres$pairwise$rowInd
      }else{
        hmres$pairwise <- aheatmap(plotdat, Colv = NA, Rowv = FALSE, breaks = 0, info = TRUE,main = hmtitle.pw) }
      heatmap.reactive$pairwise_index <- hmres$pairwise$rowInd
    }
  })
  
  ### PW subset
  output$heatmap_pairwise_subset <- renderPlot({
    heatmap.reactive$data
    input$heatmap.pairwise.corr.range
    if(!is.null(heatmap.reactive$data$pairwise)){
      plotdat <- heatmap.reactive$data$pairwise
      plotdat_col_names <- colnames(plotdat)
      min_row <- apply(plotdat,1,min)
      max_row <- apply(plotdat,1,max)
      plotdat <-  as.data.frame(plotdat[(min_row >= input$heatmap.pairwise.corr.range[1]) & (max_row <= input$heatmap.pairwise.corr.range[2]),])
      colnames(plotdat) <- plotdat_col_names
      plotdat <- as.matrix(plotdat)
      # get data for pairwise heatmap top/bottom
      if(input$heat.n.choice == "top"){
        hmdat.ss.pw <- tail( plotdat[heatmap.reactive$pairwise_index,], input$heatmap.n.genes)
      }else{
        hmdat.ss.pw <- head( plotdat[heatmap.reactive$pairwise_index,], input$heatmap.n.genes)
      }
      
      # need to manually fix the format if there's just one column present
      if(dim(plotdat)[2] == 1){
        hmdat.ss.pw <- as.matrix(hmdat.ss.pw, ncol = 1)
        colnames(hmdat.ss.pw) <- colnames(plotdat)
      }
      
      # plot heatmap
      hmres$pairwise <- aheatmap(hmdat.ss.pw, Colv = NA, Rowv = NA, breaks = 0, main = hmtitle.pw, info = TRUE)
      
      # define table output
      tbl <- as.data.frame(row.names(hmdat.ss.pw))
      colnames(tbl) <- "Protein ID"
      output$top_n_pair <- renderTable({tbl})
    }
  })
  
  ## Condition-wise heatmaps
  output$heatmap_conditionwise_all <- renderPlot({
    heatmap.reactive$data
    input$heatmap.conditionwise.corr.range
    if(!is.null(heatmap.reactive$data$condition_comp)){
      plotdat <- heatmap.reactive$data$condition_comp
      plotdat_col_names <- colnames(plotdat)
      min_row <- apply(plotdat,1,min)
      max_row <- apply(plotdat,1,max)
      plotdat <-  as.data.frame(plotdat[(max_row >= input$heatmap.conditionwise.corr.range[1]) & (min_row <= input$heatmap.conditionwise.corr.range[2]),])
      colnames(plotdat) <- plotdat_col_names
      plotdat <- as.matrix(plotdat)
      #condition-wise heatmap
      if(dim(plotdat)[2] == 1){
        hmres$condition_comp <- aheatmap(plotdat, Colv = NA, Rowv = NA, breaks = 0, main = hmtitle.cond, info = TRUE)
        heatmap.reactive$conditionwise_index<- hmres$condition_comp$rowInd
      }else{
        hmres$condition_comp <-  aheatmap(plotdat, Colv = NA, Rowv = FALSE, breaks = 0, info = TRUE, main = hmtitle.cond) 
        heatmap.reactive$conditionwise_index<- hmres$condition_comp$rowInd
      }
    }
  })
  
  ### Condition-wise subset heatmaps
  output$heatmap_conditionwise_subset <- renderPlot({
    heatmap.reactive$data
    input$heatmap.conditionwise.corr.range
    if(!is.null(heatmap.reactive$data$condition_comp)){
      plotdat <- heatmap.reactive$data$condition_comp
      plotdat_col_names <- colnames(plotdat)
      min_row <- apply(plotdat,1,min)
      max_row <- apply(plotdat,1,max)
      plotdat <-  as.data.frame(plotdat[(max_row >= input$heatmap.conditionwise.corr.range[1]) & (min_row <= input$heatmap.conditionwise.corr.range[2]),])
      colnames(plotdat) <- plotdat_col_names
      plotdat <- as.matrix(plotdat)   
      # get data for condition-wse heatmap top/bottom
      if(input$heat.n.choice == "top"){
        hmdat.ss.cond <- tail( plotdat[heatmap.reactive$conditionwise_index,], input$heatmap.n.genes)
      }else{
        hmdat.ss.cond <- head( plotdat[heatmap.reactive$conditionwise_index,], input$heatmap.n.genes)
      }
      
      # need to manually fix the format if there's just one column present
      if(dim(plotdat)[2] == 1){
        hmdat.ss.cond <- as.matrix(hmdat.ss.cond, ncol = 1)
        colnames(hmdat.ss.cond) <- colnames(plotdat)
      }
      
      # plot heatmap
      hmres$condition_comp <- aheatmap(hmdat.ss.cond, Colv = NA, Rowv = NA, breaks = 0, main = hmtitle.pw, info = TRUE)
      
      # define table output
      tbl <- as.data.frame(row.names(hmdat.ss.cond))
      colnames(tbl) <- "Protein ID"
      output$top_n_conditionwise <- renderTable({tbl})
    }
  })
  
  
  
  observe({
    # check if all mandatory fields have a value
    mandatoryFilled <-
      vapply(heatmap.reactive$fieldsMandatory,
             function(x) {
               !is.null(input[[x]]) && input[[x]] != ""
             },
             logical(1))
    mandatoryFilled <- all(mandatoryFilled)
    # enable/disable the submit button
    shinyjs::toggleState(id = "plot_heatmap", condition = mandatoryFilled)
  })  
  
  shinyjs::hide("heatmap_fluid")
  shinyjs::hide("heatmap_sub_fluid")
  
  observeEvent(input$plot_heatmap, {
    shinyjs::hide("heatmap_fluid")
    shinyjs::hide("heatmap_sub_fluid")
    closeAlert(session, "heatAlert")
    
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating.... This may take a while.", value = 0)
    
    heatmap.reactive$data <- NULL
    hmres <- NULL
    
    #define experiments to be retrieved from the DB
    # len <- length(heatmap_inserted)+1
    len <- c(gsub("txt", "", heatmap_inserted), 0)
    exps <- paste0("input$", paste0("show_expt_id_heatmap", len ))
    exps <- sapply(exps, function(x) eval(parse(text=x)))
    names(exps) <- NULL      
    
    uniquifiers <- c("gene_symbol","condition","replicate")
    
    # add information --> ideally, these should be supplied by the user during the upload/data selection!
    replicate <- paste0("input$", paste0("heatmap.rep", len))
    replicate <- sapply(replicate, function(x) eval(parse(text=x)))
    names(replicate) <- NULL      
    
    condition <- paste0("input$", paste0("heatmap.cond", len))
    condition <- sapply(condition, function(x) eval(parse(text=x)))
    names(condition) <- NULL  
    
    n.con <- length(unique(condition))
    n.rep <- length(unique(replicate))
    
    if(n.con == 1 & n.rep == 1){
      createAlert(session, "heatmapAlert", "heatAlert", title = "Error",
                  content = "You need to specify at least either two different replicates or two different conditions, otherwise there is no way to calculate and display correlations.", append = FALSE)  
      shinyjs::hide("heatmap_fluid")
      shinyjs::hide("heatmap_sub_fluid")
    }else{
      closeAlert(session, "heatAlert")
      
      # retrieve data from data.base
      hmdat <- as.data.table( DepLab:::query.measurements.by.expt.with.gene.symbol.v2(database.name.reactive$data, exps, "raw.intensity"))
      setkey(hmdat, expt_id)
      old.dim <- dim(hmdat)[1]
      
      meta.dt <- data.table(expt_id = exps, replicate = replicate, condition = condition)
      setkey(meta.dt, expt_id)
      
      hmdat <- hmdat[meta.dt]
      if(old.dim != dim(hmdat)[1]){print("Oh, oh! Something was changed during the data.table merging.")}
      
      setorder(hmdat, gene_symbol, expt_id,fraction)
      setnames(hmdat, "value","raw")
      
      # add normalized values
      hmdat[, superSmu := supsmu(x = seq(1, length(raw)), y = raw )$y , by = uniquifiers]
      hmdat[, norm.value:= superSmu * length(unique(fraction)) / sum(superSmu), by = uniquifiers]
      
      calc.mode <- ifelse( n.con == 1 | length(interaction(condition, replicate)) == n.con , "pairwise", "both" )
      
      # generate data for heatmaps
      heatmap.reactive$data  <- DepLab:::corrHM_wrap(hmdat, measurement = input$heatmap_measurement ,
                                                     corr.method = input$heatmap_cor_method,
                                                     calc.mode = calc.mode,
                                                     uniquifiers = uniquifiers)
      shinyjs::show("heatmap_fluid")
      shinyjs::show("heatmap_sub_fluid")
    }
  })
  
  
  ## UI for DB viewer
  #PAUL
  observeEvent(input$show_expt_id_db_viewer, {
    if(input$show_expt_id_db_viewer != ""){
      output$dv_viewer_summary <- renderUI({
        #    HTML(
        origin_header <- paste0("<b><font size='4'>", "Sample origin and preparation", "</font></b>")
        expt_id <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Experimental ID: ", "</b>"), input$show_expt_id_db_viewer)
        organism <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Organism: ", "</b>"),  DepLab:::list.organism.by.expt.v2(database.name.reactive$data, input$show_expt_id_db_viewer))
        origin_df <-  DepLab:::list.origin.data.by.expt.v2(database.name.reactive$data, input$show_expt_id_db_viewer)
        expterimenter <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Name of experimenter: ", "</b>"), origin_df$experimenter)
        genotype <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Genotype: ", "</b>"), origin_df$genotype)
        celltype <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Cell type: ", "</b>"), origin_df$cell_type)
        harvest_date <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Harvest date: ", "</b>"), origin_df$harvest_date)
        buffer_composition <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Buffer composition: ", "</b>"), origin_df$buffer_composition)
        lysis_method <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Lysis method: ", "</b>"), origin_df$lysis_method)
        digestion_enzyme <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Digestion enzyme: ", "</b>"), origin_df$digestion_enzyme)
        notes <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Notes: ", "</b>"), origin_df$notes)
        
        
        prefractionation_df <-  DepLab:::list.prefractionation.data.by.expt.v2(database.name.reactive$data, input$show_expt_id_db_viewer)
        prefractionation_header <- paste0("<b><font size='4'>", "Prefractionation method", "</font></b>")
        column_id <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Column ID: ", "</b>"), prefractionation_df$column_id)
        amount_protein_loaded <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Amount of protein loaded (ug): ", "</b>"), prefractionation_df$amount_protein_loaded)
        sample_vol_loaded <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Sample volume loaded (ul): ", "</b>"), prefractionation_df$sample_vol_loaded)
        lc_flow_rate <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "LC flow rate: ", "</b>"), prefractionation_df$lc_flow_rate)
        lc_fraction_size <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "LC fraction size (ml): ", "</b>"), prefractionation_df$lc_fraction_size)
        time_per_fraction <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Time per fraction (s): ", "</b>"), prefractionation_df$time_per_fraction)
        fractions_collected <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Number of fractions collected: ", "</b>"), prefractionation_df$fractions_collected)
        
        
        msmethod_df <-  DepLab:::list.msmethod.data.by.expt.v2(database.name.reactive$data, input$show_expt_id_db_viewer)
        msmethod_header <- paste0("<b><font size='4'>", "Mass spectrometry method", "</font></b>")
        instrument_id <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Instrument ID: ", "</b>"), msmethod_df$instrument_id)
        run_date <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Run date: ", "</b>"), msmethod_df$run_date)
        method_length <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Length of method (m): ", "</b>"), msmethod_df$method_length)
        
        
        dataproc_df <-  DepLab:::list.data.processing.data.by.expt.v2(database.name.reactive$data, input$show_expt_id_db_viewer)
        dataproc_header <- paste0("<b><font size='4'>", "Data processing", "</font></b>")
        processing_platform <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Processing platform: ", "</b>"), dataproc_df$processing_platform)
        search_algorithm <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Search algorithm: ", "</b>"), dataproc_df$search_algorithm)
        filtering_algorithm <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Filtering algorithm: ", "</b>"), dataproc_df$filtering_algorithm)
        filtering_stringency <- paste(paste0("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>", "Filtering stringency: ", "</b>"), dataproc_df$filtering_stringency)
        
        
        HTML(paste(origin_header, expt_id, organism, expterimenter, genotype, celltype, harvest_date, buffer_composition, lysis_method, digestion_enzyme, notes, 
                   prefractionation_header, column_id, amount_protein_loaded, amount_protein_loaded, sample_vol_loaded, lc_flow_rate, lc_fraction_size, time_per_fraction, fractions_collected,
                   msmethod_header, instrument_id, run_date, method_length,
                   dataproc_header, processing_platform, search_algorithm, filtering_algorithm, filtering_stringency, sep = '<br/>'))
        
      })
    }
  })
  
  
  ## dynamic UI for DB editor
  inserted <- character(0)
  
  observeEvent(input$EditButton, {
    if(input$show_expt_id_db_browser != ""){
      btn <- input$EditButton
      id <- paste0('txt', btn)
      updateButton(session, "EditButton", disabled = TRUE)
      updateButton(session, "deleteButton", disabled = TRUE)
      
      insertUI(
        selector = '#placeholderEditButtons',
        ui = tags$div(
          textInput("edit.expt.id", label = "Expt ID", value = input$show_expt_id_db_browser) ,
          selectInput("edit.organism", label = "Organism", choices = list("human", "yeast") ,selected = DepLab:::list.organism.by.expt.v2(database.name.reactive$data, input$show_expt_id_db_browser)  ),
          actionButton("saveEditsButton", "Save changes"),  actionButton("cancelEditsEbutton", "Cancel"),
          bsAlert("saveEditAlert"),
          id = id
        )
      )
      inserted <<- c(id, inserted)
    }
  })
  
  observeEvent(input$show_expt_id_db_browser, {
    updateTextInput(session, 'edit.expt.id', value = input$show_expt_id_db_browser)
    updateSelectInput(session, 'edit.organism', selected = DepLab:::list.organism.by.expt.v2(database.name.reactive$data, input$show_expt_id_db_browser))
  })
  
  observeEvent(input$saveEditsButton, {
    createAlert(session, "saveEditAlert", "savEdAlert", title = "Saving",
                content = "Saving changes to database.. This might take a while....", append = FALSE)
    DepLab:::update.expt(database.name.reactive$data, input$show_expt_id_db_browser, input$edit.expt.id, input$edit.organism )
    updateSelectizeInput(session, 'show_expt_id', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
    updateSelectizeInput(session, 'show_expt_id_std', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
    updateSelectizeInput(session, 'show_trypsin_symbol', choices = DepLab:::list.std.gene.symbols(database.name.reactive$data)$id, server = TRUE)
    updateSelectizeInput(session, 'show_expt_id_sum', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
    updateSelectizeInput(session, 'show_expt_id_db_browser', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
    updateSelectizeInput(session, 'show_expt_id_db_viewer', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
    updateSelectizeInput(session, 'show_expt_id_heatmap0', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
    closeAlert(session, "savEdAlert")
    createAlert(session, "saveEditAlert", "savEdAlert", title = "Saving",
                content = "Saving...Success!", append = FALSE)  
    output$expt_id_list <- renderText({ 
      expt_id_msg <- paste("The currently selected database contains the following data sets from previous sessions:", paste(as.data.frame(DepLab:::list.expt.ids.w.organism(database.name.reactive$data))$expt_id, collapse="\n"), sep = "\n")
    })
    updateButton(session, "deleteButton", disabled = FALSE)
    updateButton(session, "EditButton", disabled = FALSE)
  })
  
  
  observeEvent(input$cancelEditsEbutton, {
    removeUI(
      selector = paste0('#', inserted[length(inserted)])
      # selector = '#placeholderEditButtons'
    )
    inserted <<- inserted[-length(inserted)]
    updateSelectizeInput(session, 'show_expt_id_db_browser', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE, selected = NULL)
    updateButton(session, "deleteButton", disabled = FALSE)
    updateButton(session, "EditButton", disabled = FALSE)
  })
  
  ## dynamic UI for DB editor
  
  
  observeEvent(input$file1$name, {
    closeAlert(session, "exampleAlert")
  })
  
  observeEvent(input$expt.id, {
    if(input$expt.id != ""){
      closeAlert(session, "exampleAlert")
    }
  })
  
  observeEvent(input$deleteButton, {
    if(input$show_expt_id_db_browser != ""){
      createAlert(session, "deleteAlert", "delAlert", title = "Deleting",
                  content = "Deleting experiment from database... This might take a while....", append = FALSE)
      DepLab:::delete.expt(database.name.reactive$data, input$show_expt_id_db_browser)
      
      updateSelectizeInput(session, 'show_expt_id', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
      updateSelectizeInput(session, 'show_expt_id_std', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
      updateSelectizeInput(session, 'show_trypsin_symbol', choices = DepLab:::list.std.gene.symbols(database.name.reactive$data)$id, server = TRUE)
      updateSelectizeInput(session, 'show_expt_id_sum', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
      updateSelectizeInput(session, 'show_expt_id_db_browser', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
      updateSelectizeInput(session, 'show_expt_id_db_viewer', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
      updateSelectizeInput(session, 'show_expt_id_heatmap0', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
      
      closeAlert(session, "delAlert")
      createAlert(session, "deleteAlert", "delAlert", title = "Deleting",
                  content = "Deleting...Success!", append = FALSE)  
      
      output$expt_id_list <- renderText({ 
        expt_id_msg <- paste("The currently selected database contains the following data sets from previous sessions:", paste(as.data.frame(DepLab:::list.expt.ids.w.organism(database.name.reactive$data))$expt_id, collapse="\n"), sep = "\n")
      })
    }
  })
  
  
  observeEvent(input$saveButton, {
    closeAlert(session, "exampleAlert")
    
    if(is.null(input$file1$name)){
      createAlert(session, "alert", "exampleAlert", title = "Oops",
                  content = "Select Maxquant data file before saving.", append = FALSE)
      return()
    } else if(input$expt.id == ""){
      createAlert(session, "alert", "exampleAlert", title = "Oops",
                  content = "Experimental ID required in order to save.", append = FALSE)
      return()
    }
    
    createAlert(session, "alert", "exampleAlert", title = "Saving",
                content = "Saving... This might take a while....", append = FALSE)
    if(file.info(database.name.reactive$data)$size == 0 || is.na(file.info(database.name.reactive$data)$size))
    {
      initialize.database(database.name.reactive$data, organism = "human", force = FALSE)
      initialize.database(database.name.reactive$data, organism = "yeast", force = FALSE)
    } 
    
    # read in proteins of interest for the corresponding organism
    x <- read.MQ.data(normalizePath(paste(input$file1$datapath),winslash=.Platform$file.sep), input$expt.id, data.subset = "poi", input$organism)
    
    # read in the proteins used for standardization
    if(input$id_for_standard == "") {
      # if nothing else is given by the user, extract values for trypsin from pig,
      # which is generally used to digest the proteins before MS
      x.std <- read.MQ.data(normalizePath(paste(input$file1$datapath),winslash=.Platform$file.sep), input$expt.id, data.subset = "trypsin", organism = NULL)
    } else {
      std_list <- input$id_for_standard
      std_list <- unlist(strsplit(std_list,","))
      std_list <- gsub(" ", "", std_list, fixed = TRUE)
      x.std <- read.MQ.data(normalizePath(paste(input$file1$datapath),winslash=.Platform$file.sep), input$expt.id, data.subset = c(std_list), organism = NULL) 
    }
    
    
    origin_df <- data.frame(expt_id = input$expt.id, experimenter = input$experimenter.name, genotype = input$genotype, cell_type = input$cell.type, harvest_date = as.character(input$harvest.date), buffer_composition = input$buffer.composition, lysis_method = input$lysis.method, digestion_enzyme = input$digestion.enzyme, notes = input$notes)
    
    if (input$include.prefractionation.metadata == TRUE){
      prefractionation_df <- data.frame(expt_id = input$expt.id, column_id = input$column.id, amount_protein_loaded = input$amount.protein.loaded, sample_vol_loaded = input$sample.vol.loaded, lc_flow_rate = input$lc.flow.rate, lc_fraction_size = input$lc.fraction.size, time_per_fraction = input$time.per.fraction, fractions_collected = input$fractions.collected)
    } else {
      prefractionation_df = NULL
    }
    
    if (input$include.msmethod.metadata == TRUE){
      msmethods_df <- data.frame(expt_id = input$expt.id, instrument_id = input$instrument.id, run_date = as.character(input$run.date), method_length = input$method.length)
    } else {
      msmethods_df = NULL
    }
    
    if (input$include.data.metadata == TRUE){
      dataproc_df <- data.frame(expt_id = input$expt.id, processing_platform = input$processing.platform, search_algorithm = input$search.algorithm, filtering_algorithm = input$filtering.algorithm, filtering_stringency = input$filtering.stringency)
    } else {
      dataproc_df = NULL
    }
    
    df_readIn = try(add.expt.to.database(database.name.reactive$data, data.frame(expt_id = input$expt.id, organism = input$organism), NULL, x, x.std, origin_df, prefractionation_df, msmethods_df, dataproc_df))
    
    
    if(class(df_readIn) == "try-error"){
      print(database.name.reactive$data)
      closeAlert(session, "exampleAlert")
      createAlert(session, "alert", "exampleAlert", title = "Oops",
                  content = "Experimental ID not unique!  Saving failed.", append = FALSE)
    } else {
      updateSelectizeInput(session, 'show_expt_id', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
      updateSelectizeInput(session, 'show_expt_id_std', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
      updateSelectizeInput(session, 'show_trypsin_symbol', choices = DepLab:::list.std.gene.symbols(database.name.reactive$data)$id, server = TRUE)
      updateSelectizeInput(session, 'show_expt_id_sum', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
      updateSelectizeInput(session, 'show_expt_id_db_browser', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
      updateSelectizeInput(session, 'show_expt_id_db_viewer', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
      updateSelectizeInput(session, 'show_expt_id_heatmap0', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, server = TRUE)
      
      closeAlert(session, "exampleAlert")
      createAlert(session, "alert", "exampleAlert", title = "Saving",
                  content = "Saving...Success!", append = FALSE)  
      database.name.reactive <- reactiveValues(data = database.name)
      
      output$expt_id_list <- renderText({ 
        expt_id_msg <- paste("The currently selected database contains the following data sets from previous sessions:", paste(as.data.frame(DepLab:::list.expt.ids.w.organism(database.name.reactive$data))$expt_id, collapse="\n"), sep = "\n")
      })
      
    }
  })
  
  # ===================
  # Plot tab:
  
  
  observeEvent(input$new.complex.name, ({
    txt <- input$new.complex.name
    disabled = NULL
    style = "default"
    icon = ""
    if(txt == "") {
      disabled = TRUE
      icon <- icon("ban")
    } else {
      disabled = FALSE
      icon = icon("check")
      
    }
    updateButton(session, "saveModalButton", disabled = disabled, style = style, icon = icon)
  }))  
  
  observeEvent(input$saveModalButton, {
    newDF <- data.frame("id"= input$new.complex.name, "gene symbol"=input$show_gene_symbol, "complex"=input$new.complex.name)
    write.table(newDF, custom.complexes.name.reactive$data, append = TRUE, row.names = FALSE, col.names = FALSE, na = "NA", quote=F, sep="\t")
    toggleModal(session, "save_as_complex_modal",close)
    updateTextInput(session, "new.complex.name", value = "")     
    custom_complexes <- reactiveValues(data = read.table(file=custom.complexes.name.reactive$data,sep="\t", stringsAsFactors=FALSE, header=TRUE))
    
    if(input$complex_source == "benschop"){
      cur_genes <- DepLab:::list.all.gene.symbols.v2(database.name.reactive$data, "yeast", cur_len = cur_frac_table_length)$gene
      cur_inter <- intersect(cur_genes, benschop_complexes$data$gene_symbol)
      updateSelectizeInput(session, 'show_complex', choices = c(unique(benschop_complexes$data[match(cur_inter, benschop_complexes$data$gene_symbol),]$complex)))
    } else if (input$complex_source == "wodak"){
      cur_genes <- DepLab:::list.all.gene.symbols.v2(database.name.reactive$data, "yeast", cur_len = cur_frac_table_length)$gene
      cur_inter <- intersect(cur_genes, wodak_complexes$data$gene_symbol)    
      updateSelectizeInput(session, 'show_complex', choices = c(unique(wodak_complexes$data[match(cur_inter, wodak_complexes$data$gene_symbol),]$complex)))
    } else if (input$complex_source == "custom") {
      updateSelectizeInput(session, 'show_complex', choices = c(unique(custom_complexes$data$complex)))
    } else if (input$complex_source == "corum") {
      cur_genes <- DepLab:::list.all.gene.symbols.v2(database.name.reactive$data, "human", cur_len = cur_frac_table_length)$gene
      cur_inter <- intersect(cur_genes, corum_complexes$data$gene_symbol)
      updateSelectizeInput(session, 'show_complex', choices = c(unique(corum_complexes$data[match(cur_inter, corum_complexes$data$gene_symbol),]$complex)))
    }
    
  })
  
  
  ##
  values <- reactiveValues()
  proteasome_data <- reactiveValues(data=NULL)
  summary_data <- reactiveValues(data=NULL)
  standards_data <- reactiveValues(data=NULL)
  
  ### open palette selected upon button click
  observe({
    if (is.null(input$colorButton) || input$colorButton == 0){return()}
    values$pal <- choose_palette()
  })
  
  
  observe({
    if (is.null(input$colorButton2) || input$colorButton2 == 0){return()}
    values$pal2 <- choose_palette(n=50)
  })
  
  
  ### selectize complexes
  output$complexChoices <- renderUI({
    custom_complexes <- reactiveValues(data = read.table(file=custom.complexes.name.reactive$data,sep="\t", stringsAsFactors=FALSE, header=TRUE))
    if(input$complex_source == "benschop"){
      cur_genes <- DepLab:::list.all.gene.symbols.v2(database.name.reactive$data, "yeast", cur_len = cur_frac_table_length)$gene
      cur_inter <- intersect(cur_genes, benschop_complexes$gene_symbol)
      selectizeInput('show_complex', '', choices = c(unique(benschop_complexes[match(cur_inter, benschop_complexes$gene_symbol),]$complex)) , multiple=TRUE)
    } else if (input$complex_source == "wodak"){
      cur_genes <- DepLab:::list.all.gene.symbols.v2(database.name.reactive$data, "yeast", cur_len = cur_frac_table_length)$gene
      cur_inter <- intersect(cur_genes, wodak_complexes$gene_symbol)
      selectizeInput('show_complex', '', choices = c(unique(wodak_complexes[match(cur_inter, wodak_complexes$gene_symbol),]$complex)) , multiple=TRUE) 
    } else  if (input$complex_source == "custom"){
      selectizeInput('show_complex', '', choices = c(unique(custom_complexes$data$complex)), multiple=TRUE) 
    } else if (input$complex_source == "corum") {
      cur_genes <- DepLab:::list.all.gene.symbols.v2(database.name.reactive$data, "human", cur_len = cur_frac_table_length)$gene
      cur_inter <- intersect(cur_genes, corum_complexes$gene_symbol)
      selectizeInput('show_complex', '', choices = c(unique(corum_complexes[match(cur_inter, corum_complexes$gene_symbol),]$complex)) , multiple=TRUE) 
    }
  })
  
  ### zooming function
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  observeEvent(input$proteasomePlot_dblclick, {
    brush <- input$proteasomePlot_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  ranges_2 <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$proteasomeSummaryPlot_dblclick, {
    brush <- input$proteasomeSummaryPlot_brush
    if (!is.null(brush)) {
      ranges_2$x <- c(brush$xmin, brush$xmax)
      ranges_2$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges_2$x <- NULL
      ranges_2$y <- NULL
    }
  })
  
  ###################################################### START TRYPSIN STUFF
  vars = reactiveValues(std_protein_counter = 0, std_pro_info_table = NULL)
  
  output$retrieve_protein_std_button <- renderUI({
    if(!is.null(input$show_expt_id_std)){
      if (is.null(input$show_trypsin_symbol)){
        return()
      } else {
        actionLink("retrieve_protein_std_info", label=retrieve_protein_std_button_label())
      }
    }
  })
  
  
  
  observeEvent( input$show_trypsin_symbol, {
    if(is.null(input$show_trypsin_symbol)){
      vars$std_protein_counter <-  0
    }
  }, ignoreNULL=FALSE)
  
  
  observeEvent(input$retrieve_protein_std_info, {
    if(!is.null(input$retrieve_protein_std_info)){
      input$retrieve_protein_std_info
      isolate({
        vars$std_protein_counter <-  vars$std_protein_counter + 1
      })
    }
  })
  
  retrieve_protein_std_button_label <- reactive({
    if(vars$std_protein_counter >= 1) label <- "Refresh UniProt info..."
    else label <- "Retrieve info from UniProt..."
  })
  
  
  observeEvent( input$show_trypsin_symbol, {
    if(is.null(input$show_trypsin_symbol)){
      df <- data.frame(NULL)
      vars$std_pro_info_table <- NULL
    }
  }, ignoreNULL=FALSE)
  
  output$std_prot_table = renderDataTable({
    vars$std_pro_info_table
  },  escape = FALSE, options = list(searching = FALSE, paging = FALSE, ordering = FALSE, info = FALSE, language = list(emptyTable = " ")  )  )
  
  
  
  ranges_tryp <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$proteasomeTrypsinPlot_dblclick, {
    brush <- input$proteasomeTrypsinPlot_brush
    if (!is.null(brush)) {
      ranges_tryp$x <- c(brush$xmin, brush$xmax)
      ranges_tryp$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges_tryp$x <- NULL
      ranges_tryp$y <- NULL
    }
  })
  
  observeEvent(input$show_expt_id_std, {
    updateSelectizeInput(session, 'show_trypsin_symbol', selected = c(DepLab:::list.selected.std.gene.symbols(database.name.reactive$data, input$show_expt_id_std)$id, input$show_trypsin_symbol), server = FALSE)
  })
  
  observeEvent( input$y_axis_choices_tryp, {
    if (input$y_axis_choices_tryp == "raw.intensity"){
      updateCheckboxGroupInput(session, "y_axis_log_tryp",  choices = c("log2", "normalize across fractions"))
    } else {
      updateCheckboxGroupInput(session, "y_axis_log_tryp",  choices = c("log2"))
    }
  })
  
  
  #    output$plot_trypsin <- renderPlotly({
  output$plot_trypsin <- renderPlot({
    # only render plot if gene symbol selected
    if(!is.null(input$show_expt_id_std)){
      if (is.null(input$show_trypsin_symbol)){
        return()
      } else {
        trypsin.dat <- as.data.frame(DepLab:::query.std.measurements.v2(database.name.reactive$data, c(input$show_expt_id_std), c(input$show_trypsin_symbol), c(input$y_axis_choices_tryp)))
        validate(
          need(nrow(trypsin.dat) > 0 , "UniProt ID(s) not found in any experiment.")
        )
        #  standards_data$data <- trypsin.dat
        if (is.null(values$pal_tryp)){
          color.palette=NULL
        } else {
          color.palette = values$pal_tryp
        }
        if (input$split_by_std == "none" ){
          split.by = NULL
        } else {
          split.by = input$split_by_std
        }
        if (input$split_by_std_col == "none" ){
          split.by.col = "."
        } else {
          split.by.col = input$split_by_std_col
        }
        if(any(grepl("fractions", c(input$y_axis_log_tryp)))) {
          trypsin.dat <- normalize_values(long.df = trypsin.dat,norm.type = "fraction", prot.identifier = "id")
        }
        if(any(grepl("log2", c(input$y_axis_log_tryp)))) {
          trypsin.dat$value <- log2( trypsin.dat$value +1)
          trypsin.dat$measurement <-  paste(trypsin.dat$measurement , "log2", sep= "_")
        }
        standards_data$data <- trypsin.dat
        plot_profile(trypsin.dat,  what = c("id","expt_id"),
                     color.palette=color.palette,split.by.col=split.by.col,
                     y.lab=head(trypsin.dat$measurement, 1),
                     split.by=split.by, color.by=input$color_by_choices_std, title=input$plot_tryp_title)  + coord_cartesian(xlim = ranges_tryp$x, ylim = ranges_tryp$y)
        #      #   p %>% ggplotly()
        #         p <- ggplotly(p)
        #         p %>% 
        #           layout(showlegend = T)
      } 
    }
  }, height=function() { session$clientData$output_plot_trypsin_width * 0.8 })
  
  
  observeEvent( input$retrieve_protein_std_info, {
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Retrieving info from UniProt....", value = 0)
    protDf <- read.table(text = unlist(lapply(as.list(isolate(unique(input$show_trypsin_symbol))), function(x) get_UniProt_info(x))), sep="")
    urls <- paste("http://www.uniprot.org/uniprot/", protDf$V1, sep="")
    refs <- paste0(paste0(paste0("<a href='",  urls, "' target='_blank'>"), "", protDf$V1), "", "</a>")
    protDf$V1 <- refs
    protDf$V5 <- paste(protDf$V5, "AA")
    protDf <- protDf[,c(-2, -6)]
    colnames(protDf) <- c("protein", "organism", "status", "size")
    vars$std_pro_info_table <- protDf[order(protDf$protein),]
  })
  
  
  observe({
    if (is.null(input$colorButton_tryp) || input$colorButton_tryp == 0){return()}
    values$pal_tryp <- choose_palette(n=3)
  })
  
  ###################################################### END TRYPSIN STUFF
  
  updateSelectizeInput(session, 'show_expt_id', choices = DepLab:::list.expt.ids.v2(database.name)$expt_id, server = TRUE)
  updateSelectizeInput(session, 'show_expt_id_sum', choices = DepLab:::list.expt.ids.v2(database.name)$expt_id, server = TRUE)
  updateSelectizeInput(session, 'show_expt_id_std', choices = DepLab:::list.expt.ids.v2(database.name)$expt_id, server = TRUE)
  updateSelectizeInput(session, 'show_trypsin_symbol', choices = DepLab:::list.std.gene.symbols(database.name)$id, server = TRUE)
  updateSelectizeInput(session, 'show_expt_id_db_browser', choices = DepLab:::list.expt.ids.v2(database.name)$expt_id, server = TRUE)
  updateSelectizeInput(session, 'show_expt_id_db_viewer', choices = DepLab:::list.expt.ids.v2(database.name)$expt_id, server = TRUE)
  updateSelectizeInput(session, 'show_expt_id_heatmap0', choices = DepLab:::list.expt.ids.v2(database.name)$expt_id, server = TRUE)
  
  
  current_org <- reactiveValues(counter=0 )
  
  observeEvent(input$show_expt_id, {
    if(is.null(input$show_expt_id)){
      output$current_organism <- renderText({paste("")})
      isolate(updateSelectizeInput(session, 'show_expt_id', choices = DepLab:::list.expt.ids.v2(database.name.reactive$data)$expt_id, selected=NULL, server = TRUE))
      isolate(updateSelectizeInput(session, 'show_gene_symbol', choices = NULL, selected=NULL, server = TRUE))
      current_org$counter = 0
    } 
  },ignoreNULL=FALSE)  
  
  observeEvent(input$show_expt_id, {
    if (current_org$counter == 0){
      current_org$counter =  1
      cur_sel = input$show_expt_id
      cur_sel_genes =  input$show_gene_symbol
      cur_org <- head(DepLab:::list.organism.by.expt.v2(database.name.reactive$data, input$show_expt_id)$organism, 1)
      output$current_organism <- renderText({paste("Data from:", cur_org)})
      
      isolate(updateSelectizeInput(session, 'show_expt_id', choices = DepLab:::list.expt.by.organism.v2(database.name.reactive$data, cur_org)$expt_id, selected=cur_sel, server = FALSE))
      
      isolate(updateSelectizeInput(session, 'show_gene_symbol', choices = DepLab:::list.all.gene.symbols.v2(database.name.reactive$data, cur_org, cur_len = cur_frac_table_length)$gene, selected=cur_sel_genes, server = TRUE))
      # isolate(updateSelectizeInput(session, 'show_gene_symbol', choices =setNames( DepLab:::list.all.gene.symbols.v2(database.name.reactive$data, cur_org, cur_len = cur_frac_table_length)$gene, paste("hey",  DepLab:::list.all.gene.symbols.v2(database.name.reactive$data, cur_org, cur_len = cur_frac_table_length)$gene, sep="-")), selected=cur_sel_genes, server = TRUE))
      
    }
  },ignoreNULL=TRUE)  
  
  
  
  
  
  output$plot_proteasome <- renderPlot({
    # only render plot if gene symbol selected
    if(!is.null(reac$show_expt_id)){
      if (is.null(reac$show_gene_symbol) & is.null(reac$show_complex)){
        return()
      } else {
        custom_complexes <- reactiveValues(data = read.table(file=custom.complexes.name.reactive$data,sep="\t", stringsAsFactors=FALSE, header=TRUE))
        if(input$complex_source == "benschop"){
          complex_select <- grep(  paste(paste(paste("^", make.names(reac$show_complex), sep=""), "$", sep=""), collapse="|") , make.names(benschop_complexes$complex))
          proteasome.dat <- as.data.frame(DepLab:::query.measurements.v2(database.name.reactive$data, reac$show_expt_id, c(reac$show_gene_symbol, benschop_complexes[complex_select,]$gene_symbol), input$y_axis_choices))
          complex_list <- benschop_complexes[complex_select,]
        } else if (input$complex_source == "wodak"){
          complex_select <- grep(  paste(paste(paste("^", make.names(reac$show_complex), sep=""), "$", sep=""), collapse="|") , make.names(wodak_complexes$complex))
          proteasome.dat <- as.data.frame(DepLab:::query.measurements.v2(database.name.reactive$data, reac$show_expt_id, c(reac$show_gene_symbol, wodak_complexes[complex_select,]$gene_symbol), input$y_axis_choices))
          complex_list <- wodak_complexes[complex_select,]
        } else if (input$complex_source == "corum"){
          complex_select <- grep(  paste(paste(paste("^", make.names(reac$show_complex), sep=""), "$", sep=""), collapse="|") , make.names(corum_complexes$complex))
          proteasome.dat <- as.data.frame(DepLab:::query.measurements.v2(database.name.reactive$data, reac$show_expt_id, c(reac$show_gene_symbol, corum_complexes[complex_select,]$gene_symbol), input$y_axis_choices))
          complex_list <- corum_complexes[complex_select,]
        } else {
          complex_select <- grep(  paste(paste(paste("^", make.names(reac$show_complex), sep=""), "$", sep=""), collapse="|") , make.names(custom_complexes$data$complex))
          proteasome.dat <- as.data.frame(DepLab:::query.measurements.v2(database.name.reactive$data, reac$show_expt_id, c(reac$show_gene_symbol, custom_complexes$data[complex_select,]$gene_symbol), input$y_axis_choices))
          complex_list <- custom_complexes$data[complex_select,]
        }
        validate(
          need(nrow(proteasome.dat) > 0 , "Gene(s) or complex(es) not found in any experiment.")
        )
        if (is.null(values$pal)){
          color.palette=NULL
        } else {
          color.palette = values$pal
        }
        if (input$split_by == "none" ){
          split.by = NULL
        } else {
          split.by = input$split_by
        }
        if (input$split_by_col == "none" ){
          split.by.col = "."
        } else {
          split.by.col = input$split_by_col
        }
        if (input$split_by == "complex" || input$split_by_col == "complex" || input$color_by_choices == "complex" || nrow(complex_list) > 0 ){
          validate(
            need(nrow(complex_list) > 0 , "Can't split by complex unless one or more complexes is selected. ")
          )
          proteasome.dat <- merge(proteasome.dat, complex_list, by=c("gene_symbol"), all=TRUE)
          
          #expt
          not_pres <- proteasome.dat[is.na(proteasome.dat$expt_id),] 
          not_pres <- rbind(not_pres, data.frame(gene_symbol="filler", fraction="", value="", measurement="", expt_id="", organism = "", id.x="", id.y="", complex=unique(complex_list$complex)))
          not_pres <- not_pres[order(not_pres$complex),]
          missing <- unlist(lapply(unique(not_pres$complex), function(x)  paste(not_pres[not_pres$complex == x, ]$gene_symbol, collapse=", ")))
          complex_list <- complex_list[order(complex_list$complex),]
          num_total <- dcast(complex_list, complex ~ complex, value.var="gene_symbol",  fun.aggregate=length, fill = 0)
          rownames(num_total) <- colnames(num_total[-1])
          num_total <- num_total[-1]
          num_total <- as.matrix(num_total)
          complex_stats <- data.frame(complex=unique(not_pres$complex), total = diag(num_total), missing = missing)
          complex_stats$present <- complex_stats$total - str_count(complex_stats$missing, "\\)") 
          complex_stats$missing <- gsub(", filler", "", complex_stats$missing)
          complex_stats$missing <- gsub("filler", "NA", complex_stats$missing)
          vars$complex_info_table <- complex_stats
          
          
          proteasome.dat <- proteasome.dat[!is.na(proteasome.dat$expt_id),] 
          proteasome.dat$id.y <- NULL
          colnames(proteasome.dat) <- c("gene_symbol", "fraction", "value", "measurement", "expt_id", "organism","id", "complex")
          
        }
        if(input$pointShape == "complex"){
          validate(
            need(nrow(complex_list) > 0 , "Can't assign point shape to complex unless one or more complexes is selected. ")
          )
        }
        
        if(input$lineType == "complex"){
          validate(
            need(nrow(complex_list) > 0 , "Can't assign line type to complex unless one or more complexes is selected. ")
          )
        }
        if(any(grepl("SuperSmoother", c(input$y_axis_log)))){
          # then get superSmu value
          #uniquifiers <- c("gene_symbol","condition","replicate")
          #
          proteasome.dat <- superSmooth_values(long.df = proteasome.dat, prot.identifier = "gene_symbol")
        }
        if(any(grepl("spike", c(input$y_axis_log)))) {
          selected_ctrl_names <- paste("input$radio_var", seq(1,length(reac$show_expt_id)), sep="")
          selected_ctrl_list <- NULL
          for (i in 1:length(selected_ctrl_names) ) {
            selected_ctrl_list <- c(selected_ctrl_list, eval(parse(text=selected_ctrl_names[i])))
          }
          prot.stds <- as.data.frame(DepLab:::query.std.measurements.v2(database.name.reactive$data, reac$show_expt_id, selected_ctrl_list, input$y_axis_choices) )
          proteasome.dat <- normalize_values(long.df = proteasome.dat, norm.type = "spike-in", std.df = prot.stds)
        }
        if(any(grepl("fractions", c(input$y_axis_log)))) {
          proteasome.dat <- normalize_values(long.df = proteasome.dat,norm.type = "fraction", prot.identifier = "gene_symbol")
        }
        if(any(grepl("log2", c(input$y_axis_log)))) {
          proteasome.dat$value <- log2( proteasome.dat$value +1)
          proteasome.dat$measurement <-  paste(proteasome.dat$measurement, "log2", sep= "_")
        }
        if(input$specify_replicates == FALSE & (input$split_by == "replicate" | input$split_by_col == "replicate" | input$color_by_choices == "replicate" | input$pointShape == "replicate"  | input$lineType == "replicate" | input$split_by == "condition" | input$split_by_col == "condition" | input$color_by_choices == "condition" | input$pointShape == "condition"  | input$lineType == "condition") ){
          validate(
            need(input$specify_replicates == TRUE , "Must specify conditions and replicates before faceting by them.  Check 'Specify replicates' and fill in the boxes before proceeding. ")
          )
        }
        if(input$specify_replicates == TRUE & (input$split_by == "replicate" | input$split_by_col == "replicate" | input$color_by_choices == "replicate" | input$pointShape == "replicate"  | input$lineType == "replicate" | input$split_by == "condition" | input$split_by_col == "condition" | input$color_by_choices == "condition" | input$pointShape == "condition"  | input$lineType == "condition") ){
          selected_cond_idx <- paste("input$rep_nm_cond", seq(1,length(reac$show_expt_id)), sep="")
          selected_rep_idx <- paste("input$rep_nm_rep", seq(1,length(reac$show_expt_id)), sep="")
          selected_cond_names <- NULL
          selected_rep_names <- NULL
          
          for (i in 1:length(selected_cond_idx) ) {
            selected_cond_names <- c(selected_cond_names, eval(parse(text=selected_cond_idx[i])))
          }
          
          validate(
            need(any(selected_cond_names == "") == FALSE , "Missing condition names.  Please specify a condition for each experimental ID.")
          )
          
          for (i in 1:length(selected_rep_idx) ) {
            selected_rep_names <- c(selected_rep_names, eval(parse(text=selected_rep_idx[i])))
          }
          validate(
            need(any(is.na(selected_rep_names)) == FALSE, "Missing replicate number.  Please specify a replicate numberfor each condition.")
          )
          rep_cond_df <- data.frame(expt_id=reac$show_expt_id, condition=selected_cond_names, replicate=as.character(selected_rep_names))
          proteasome.dat <- merge(proteasome.dat,rep_cond_df, by.x="expt_id", by.y="expt_id")
        }
        
        proteasome_data$data <- proteasome.dat
        plot_profile(proteasome.dat, line.smooth=any(grepl("spline", c(input$y_axis_log))),
                     color.palette=color.palette, y.lab = head(proteasome.dat$measurement, 1), 
                     split.by=split.by , split.by.col=split.by.col, 
                     color.by=input$color_by_choices, line.type=input$lineType, 
                     line.size=as.numeric(input$lineSize), point.shape=input$pointShape, 
                     point.size=as.numeric(input$pointSize), title=input$plot_title)  + coord_cartesian(xlim = ranges$x, ylim = ranges$y)
      }
    }
  }, height=function() { session$clientData$output_plot_proteasome_width * 0.8 })
  
  
  observeEvent( input$show_gene_symbol, {
    if(is.null(input$show_gene_symbol) & is.null(input$show_complex)){
      #   df <- data.frame(NULL)
      vars$pro_info_table <- NULL
      vars$complex_info_table <- NULL
    }
  }, ignoreNULL=FALSE)
  
  
  observeEvent( input$show_complex, {
    if(is.null(input$show_gene_symbol) & is.null(input$show_complex)){
      vars$pro_info_table <- NULL
      vars$complex_info_table <- NULL
    }
  }, ignoreNULL=FALSE)
  
  output$prot_table = renderDataTable({
    vars$pro_info_table
  },  escape = FALSE, options = list(searching = FALSE, paging = FALSE, ordering = FALSE, info = FALSE, language = list(emptyTable = " ")  )  )
  
  output$complex_table = renderDataTable({
    if(!is.null(input$show_complex)){
      vars$complex_info_table
    }
  },  escape = FALSE, options = list(searching = FALSE, paging = FALSE, ordering = FALSE, info = FALSE, language = list(emptyTable = " ")  )  )
  
  observeEvent( input$retrieve_protein_info, {
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Retrieving info from UniProt....", value = 0)
    dfProt <- isolate(proteasome_data$data[,c(1,6,7)])
    dfProt <- unique(dfProt)
    dfProt$gene_symbol <- stri_extract_last_words(dfProt$gene_symbol)
    if(dfProt$organism == "yeast") {
      unprot_df <- readLines("./yeast.txt")
      unprot_df <- gsub("; ", "", unprot_df) 
      unprot_df <- gsub("\\s+", " ", str_trim(unprot_df))
      unprot_df <- strsplit(unprot_df, " ")
      indx <- sapply(unprot_df, length)
      res <- as.data.frame(do.call(rbind,lapply(unprot_df, `length<-`,max(indx))))
      findUniProt <- subset(res, V2 %in% c(dfProt$id))[,c(3,2)]
      dfProt <- merge(dfProt,findUniProt, by.x="id", by.y="V2")
      colnames(dfProt) <- c("id", "gene_symbol","organism", "protein")
      protDf <- read.table(text = unlist(lapply(as.list(unique(findUniProt$V3)), function(x) get_UniProt_info(x))), sep="")
      protDf <- merge(dfProt, protDf, by.x="protein", by.y="V1")
      urls <- paste("http://www.uniprot.org/uniprot/", protDf$protein, sep="")
      refs <- paste0(paste0(paste0("<a href='",  urls, "' target='_blank'>"), "", protDf$protein), "", "</a>")
      protDf$protein <- refs
      protDf$V5 <- paste(protDf$V5, "AA")
      protDf <- protDf[,c(-4,-5,-9)]
      protDf <- protDf[,c(3,2,1,4,5,6)]
      colnames(protDf) <- c("gene", "id", "protein", "organism", "status", "size")
      vars$pro_info_table <- protDf[order(protDf$gene),]
    } else if(dfProt$organism == "human"){
      protDf <- read.table(text = unlist(lapply(as.list(isolate(unique(dfProt$id))), function(x) get_UniProt_info(x))), sep="")
      dfProt <- merge(dfProt,protDf, by.x="id", by.y="V1")
      urls <- paste("http://www.uniprot.org/uniprot/", protDf$V1, sep="")
      refs <- paste0(paste0(paste0("<a href='",  urls, "' target='_blank'>"), "", protDf$V1), "", "</a>")
      dfProt$id <- refs
      dfProt$V5 <- paste(dfProt$V5, "AA")
      dfProt <- dfProt[,c(-3, -4, -8)]
      colnames(dfProt) <- c("protein", "gene", "organism", "status", "size")
      dfProt$organism <- sub(".*_", "", dfProt$organism)
      dfProt <- dfProt[,c(2,1,3,4,5)]
      vars$pro_info_table <- dfProt[order(dfProt$gene),]
    }
    #  vars$complex_info_table <- vars$complex_info_table_tmp
    
  })
  
  
  output$plot_summary <- renderPlot({
    if(!is.null(reacSum$show_expt_id_sum)){
      #   if(isolate(reacSum$redrawSum)){
      if (is.null(values$pal2)){
        color.palette=NULL
      } else {
        color.palette = values$pal2
      }
      if (input$split_by_summary == "none" ){
        split.by = NULL
      } else {
        split.by = input$split_by_summary
      } 
      long.df <- as.data.frame(DepLab:::query.measurements.by.expt.v2(database.name.reactive$data, reacSum$show_expt_id_sum, input$y_axis_choices_sum))#[,c(2,3,5)]
      
      if(any(grepl("log2", c(input$y_axis_log_2)))) {
        long.df$value <- log2( long.df$value +1)
        #   long.df$measurement <-  paste(long.df$measurement , "log2", sep= "_")
      }
      
      summary_data$data <- long.df
      plot_cumulativeValues(long.df = long.df, split.by = split.by, color.palette=color.palette, y.lab = input$y_axis_choices_sum ,title=input$plot_summary_title, x.interval=input$summary_x_interval) + coord_cartesian(xlim = ranges_2$x, ylim = ranges_2$y)
    }
    # }
  })
  
  
  # ===================
  # Table tab: Render searchable/subsettable tables
  output$proteasomedata <- renderDataTable({
    validate(
      need(nrow(proteasome_data$data) > 0 , "Select some data on corresponding plot page.")
    )
    proteasome_data$data 
  })
  
  output$summarydata <- renderDataTable({
    validate(
      need(nrow(summary_data$data ) > 0 , "Select some data on corresponding plot page.")
    )
    summary_data$data 
  })
  
  #   
  output$save_sum_table <- downloadHandler(
    filename = function() {
      paste("summary", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(summary_data$data , file, col.names=TRUE, row.names=FALSE, quote=FALSE)
    }
  )
  
  #  shinyFileSave(input, 'save_sum_table',  roots=volumes, session=session)
  #  saveFileName <- renderPrint({parseSavePath(volumes, input$save_sum_table)})
  #  observeEvent( input$save_sum_table, {
  #    savePath <- strsplit(saveFileName(), " csv ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][2]
  #    write.csv(summary_data$data , savePath, col.names=TRUE, row.names=FALSE, quote=FALSE)
  #  })
  
  output$save_indiv_table <- downloadHandler(
    filename = function() {
      paste("individual_profile", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(proteasome_data$data , file, col.names=TRUE, row.names=FALSE, quote=FALSE)
    }
  )
  
  
  #  shinyFileSave(input, 'save_indiv_table', roots=volumes, session=session)
  #  saveFileName <- renderPrint({parseSavePath(volumes, input$save_indiv_table)})
  #  observeEvent( input$save_indiv_table, {
  #    savePath <- strsplit(saveFileName(), " csv ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][2]
  #    write.csv(proteasome_data$data, savePath, col.names=TRUE, row.names=FALSE, quote=FALSE)
  #  })
  
  output$save_std_table <- downloadHandler(
    filename = function() {
      paste("standards", "csv", sep = ".")
    },
    content = function(file) {
      write.csv(standards_data$data , file, col.names=TRUE, row.names=FALSE, quote=FALSE)
    }
  )
  #  shinyFileSave(input, 'save_std_table', roots=volumes, session=session)
  #  saveFileName <- renderPrint({parseSavePath(volumes, input$save_std_table)})
  #  observeEvent( input$save_std_table, {
  #    savePath <- strsplit(saveFileName(), " csv ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][2]
  #    write.csv(standards_data$data, savePath, col.names=TRUE, row.names=FALSE, quote=FALSE)
  #  })
  
  output$standardsedata <- renderDataTable({
    validate(need(
      nrow(standards_data$data) > 0 , "Select some data on corresponding plot page."
    ))
    standards_data$data
  })
  
  
})

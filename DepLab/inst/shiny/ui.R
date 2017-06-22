# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggplot2)
library(scales)
library(colorspace)
#library(plotly)
library(shinyBS)
library(shinyFiles)
library(shinythemes)
library(plyr)
library(shinyjs)


shinyUI(  #            navbarPage(theme = shinytheme("Flatly"), "PCP",
  navbarPage( "PCP",
              
              tabPanel("Data Import",
                       h1("SELECT EXISTING DATABASE"),
                       verbatimTextOutput('db_dir_path'),
                       verbatimTextOutput('expt_id_list'),
                       shinyjs::useShinyjs(),
                       shinyFilesButton( 'db_path' ,  'Change database' ,  'Please select database file' , FALSE),
                       br(), br(), 
                       #textOutput('expt_id_list'),
                       verbatimTextOutput('custom_complexes_dir_path'),
                       shinyFilesButton( 'custom_complexes_path' ,  'Change file for custom complexes' ,  'Please select file containing custom complexes' , FALSE),
                       br(), br(), 
                       helpText("If you're fine with these selections and you do not wish to upload new data, simply continue to the visualization tab."),
                       br(),
                       h1("UPLOAD NEW EXPERIMENT"),
                       h2("Upload Maxquant data"),
                       br(),
                       fluidRow(column(4, textInput("expt.id", label = "Expt ID", value = "")),
                                column(4,selectInput("organism", label = "Organism", choices = list("human", "yeast"),selected = "human"))
                       ),
                       fileInput('file1', '',
                                 accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                       helpText("The expected input file is the proteinGroups.txt output from MaxQuant."),
                       br(),
                       fluidRow(
                         column(12,
                                textInput("id_for_standard", 
                                          label = "Enter comma-separated list of UniProt ID(s) to use as standard (optional).", value = ""),
                                helpText("E.g., P00761 for trypsin (pig), P02769 for BSA."),
                                tags$head(tags$style(HTML('#saveButton{background-color: #E77471}')) ),
                                actionButton("saveButton", "Save") ),
                         column(5, bsAlert("inputalert") , bsAlert("alert") )),
                       
                       #br(),
                       #br(),
                       #br(),
                       h2("Sample origin and preparation"),
                       selectInput("experimenter.name",
                                   choices  = list(Paola="Paola Cavaliere", Noah="Noah Dephoure", Vijay="Vijay Raja", Nadia="Nadia Iqbal"),
                                   label = "Name of experimenter",
                                   selected = "Paola Cavaliere"),
                       fluidRow(
                         column(4,selectInput("genotype",
                                              label = "Genotype",
                                              choices = list(EV="Empty vector", WT="Wild type", D109N="D109N"),
                                              selected = "Wild type")),
                         column(4,selectInput("cell.type",
                                              label = "Cell type",
                                              choices = list(MCF10A="MCF10A", MDA231="MDA231"),
                                              selected = "MCF10A"))
                       ),
                       fluidRow(
                         column(4,dateInput("harvest.date",
                                            label="Harvest date")),
                         column(4,selectInput("buffer.composition",
                                              label = "Buffer composition",
                                              choices = list(Tris_HCl="TrisHCl 50mM pH 7.5 KCl 150mM",
                                                             Hepes = "Hepes 150mM pH 7.9 MgCl2 1.5mM Ammonium acetate 150mM",
                                                             None = "NA"),
                                              selected = "TrisHCl 50mM pH 7.5 KCl 150mM"))),
                       
                       fluidRow(
                         column(4,selectInput("lysis.method",
                                              label = "Lysis method",
                                              choices = list("Dounce","NP40","Sonication", "None"),
                                              selected = "Dounce")),
                         # column(4,selectInput("protein.assay.method",
                         #                        label = "Protein assay method",
                         #                        choices = list("none"),
                         #                        selected = "None")),
                         # column(4,numericInput("protein.concentration",
                         #                       label = "Concentration (uM)",
                         #                       value = 1)),
                         column(4,selectInput("digestion.enzyme",
                                              label = "Digestion enzyme",
                                              choices = list("Trypsin", "Lysine C","Trypsin_LysC"),
                                              selected = "Trypsin"))
                       ),
                       
                       
                       fluidRow(
                         column(4,textInput("notes",
                                            label = "Notes",
                                            value = ""))
                       ),
                       
                       h2("Prefractionation method"), checkboxInput("include.prefractionation.metadata", "Include?", value = FALSE),
                       textInput("column.id",label = "Column ID",value = ""),
                       fluidRow(
                         column(4,numericInput("amount.protein.loaded",
                                               label = "Amount of protein loaded (ug)",
                                               value = 1)),
                         column(4,numericInput("sample.vol.loaded",
                                               label = "Sample volume loaded (ul)",
                                               value = 1)),
                         column(4,numericInput("lc.flow.rate",
                                               label = "LC flow rate",
                                               value = 1)),
                         column(4,numericInput("lc.fraction.size",
                                               label = "LC fraction size (ml)",
                                               value = 1)),
                         column(4,numericInput("time.per.fraction",
                                               label = "Time per fraction (s)",
                                               value = 1)),
                         column(4,numericInput("fractions.collected",
                                               label = "Number of fractions collected",
                                               value = 1))
                       ),
                       #          fileInput('chrom_image', 'Choose Chromatogram Image', accept=c('image/png')),
                       #          imageOutput("preImage", inline = TRUE),
                       
                       h2("Mass spectrometry method"), checkboxInput("include.msmethod.metadata", "Include?", value = FALSE),
                       textInput("instrument.id",label = "Instrument ID",value = "Orbitrap fusion"),
                       fluidRow(
                         column(4,dateInput("run.date",
                                            label = "Run date")),
                         # column(4,textInput("data.file.name",
                         #                    label = "Data file",
                         #                    value = "")),
                         # column(4,textInput("method.file.name",
                         #                    label = "Method file",
                         #                    value = "")),
                         column(4,numericInput("method.length",
                                               label = "Length of method (m)",
                                               value = 65))#,
                         # column(4,selectInput("quantitative.method",
                         #                      label = "Quantitative method",
                         #                      choices = list("LFQ"),
                         #                      selected = "LFQ"))
                       ),
                       
                       h2("Data processing"), checkboxInput("include.data.metadata", "Include?", value = FALSE),
                       fluidRow(
                         column(4,selectInput("processing.platform",
                                              label = "Processing platform",
                                              choices = list("core", "maxquant", "OPEN-MS"),
                                              selected = "maxquant")),
                         column(4,textInput("search.algorithm",
                                            label = "Search algorithm",
                                            value = "")),
                         column(4,textInput("filtering.algorithm",
                                            label = "Filtering algorithm",
                                            value = "")),
                         column(4,textInput("filtering.stringency",
                                            label = "Filtering stringency",
                                            value = ""))
                       )
                       
                       # contains the reactively rendered output, dependent on selected widget elements
                       #     h2("Summary"),
                       #    br()
                       
              ),
              tabPanel("Visualization",
                       tabsetPanel("Panel 1.x",
                                   
                                   tabPanel("Summary plots",
                                            sidebarLayout(
                                              
                                              # Sidebar with a slider input
                                              wellPanel(
                                                tags$style(type="text/css", '#leftPanel { width:250px; float:left;}'),
                                                id = "leftPanel",
                                                helpText("Select the relevant variable conditions"),
                                                selectizeInput('show_expt_id_sum', 'Select expt ID', choices = NULL,  selected =  NULL, multiple=TRUE),  
                                                
                                                selectInput("y_axis_choices_sum", label = ("y variable"), 
                                                            choices = list("Raw intensity" = "raw.intensity",
                                                                           #"LFQ intensity" = "LFQ.intensity", 
                                                                           "MS count" = "MS.MS.count", 
                                                                           "Peptides count"= "peptides.count", 
                                                                           "Unique peptides only" ="unique.peptides.only",
                                                                           "Razor and unique peptides" = "razor.and.unique.peptides" #, "Sequence coverage" ="sequence.coverage"
                                                            ), 
                                                            selected = "raw.intensity") ,
                                                checkboxGroupInput("y_axis_log_2",
                                                                   label = "y transformation",
                                                                   choices = "log2"),
                                                selectInput("split_by_summary", label = ("Split into rows by"), 
                                                            choices = list("None" = "none", "Expt ID" = "expt_id"), 
                                                            selected = "expt_id"),
                                                sliderInput("summary_x_interval", "x interval", 1, 50, 5),
                                                textInput("plot_summary_title", label = "Plot title", value = ""),
                                                actionButton("colorButton2", "Change colors"),
                                                downloadButton('downloadPlot', "Save plot")
                                                # shinySaveButton('save', 'Save plot', 'Save plot as...', filetype=list(pdf=c('pdf'), png=c('png')))
                                                
                                              ),
                                              # Show a plot of the generated distribution
                                              mainPanel(
                                                plotOutput("plot_summary",
                                                           #  height ="auto",
                                                           height = 600,
                                                           width="100%",
                                                           dblclick = "proteasomeSummaryPlot_dblclick",
                                                           brush = brushOpts(
                                                             id = "proteasomeSummaryPlot_brush",
                                                             resetOnNew = TRUE))
                                              )
                                            )),
                                   
                                   
                                   
                                   
                                   
                                   tabPanel("Individual proteins",            
                                            sidebarLayout(
                                              # side panel contains the variable selection widgets
                                              wellPanel(
                                                tags$style(type="text/css", '#leftPanel { width:250px; float:left;}'),
                                                id = "leftPanel",
                                                helpText("Select the relevant variable conditions"),
                                                selectizeInput('show_expt_id', 'Select expt ID', choices = NULL,  selected =  NULL, multiple=TRUE),
                                                textOutput("current_organism"),
                                                checkboxInput("specify_replicates", "Specify replicates", value = FALSE),
                                                uiOutput("uiRepList"),
                                                # gene symbol to view
                                                selectizeInput('show_gene_symbol', 'Select gene symbol', choices = NULL,  selected =  NULL, multiple=TRUE),  
                                                actionLink("save_as_complex", "save as complex..."),
                                                br(),
                                                actionLink("paste", "                paste from clipboard..."),
                                                
                                                div(style="margin-top:10px; width:100%;", 
                                                    radioButtons(inputId="complex_source", label="Select complex from", choices=list("Benschop (Sc)" = "benschop","Wodak (Sc)" = "wodak", "CORUM core (Hs)" = "corum", "Custom" = "custom"), inline=FALSE)
                                                ),
                                                ### bs modal for save group of genes as complex
                                                bsModal("save_as_complex_modal", "Save as...", "save_as_complex", size="small", 
                                                        textInput("new.complex.name",
                                                                  label = "Save complex as ",
                                                                  value = ""),
                                                        bsButton("saveModalButton", label = "Save",  icon = icon("ban"))
                                                ),
                                                htmlOutput("complexChoices"),
                                                selectInput("y_axis_choices", label = ("y variable"), 
                                                            choices = list("Raw intensity" = "raw.intensity",
                                                                           #"LFQ intensity" = "LFQ.intensity", 
                                                                           "MS count" = "MS.MS.count", 
                                                                           "Peptides count"= "peptides.count", 
                                                                           "Unique peptides only" ="unique.peptides.only",
                                                                           "Razor and unique peptides" = "razor.and.unique.peptides"#,
                                                                           #"Sequence coverage" ="sequence.coverage"
                                                            ), 
                                                            selected = "raw.intensity"),
                                                checkboxGroupInput("y_axis_log",
                                                                   label = "y transformation",
                                                                   choices = c("log2", "normalize across fractions", "normalize by spike-in", "smooth values (SuperSmoother)", "smooth lines (spline)")),
                                                uiOutput("uiRadio"),
                                                selectInput("split_by", label = ("Split into rows by"), 
                                                            choices = list("None" = "none", "Gene symbol" = "gene_symbol", "ID" = "id", "Expt ID" = "expt_id", "Complex" = "complex",  "Replicate" = "replicate", "Condition" = "condition"), 
                                                            selected = "expt_id"),
                                                selectInput("split_by_col", label = ("Split into columns by"), 
                                                            choices = list("None" = "none", "Gene symbol" = "gene_symbol", "ID" = "id",  "Expt ID" = "expt_id", "Complex" = "complex",  "Replicate" = "replicate", "Condition" = "condition"), 
                                                            selected = "none"),
                                                selectInput("color_by_choices", label = ("Color by"), 
                                                            choices = list("Gene symbol" = "gene_symbol", "ID" = "id", "Expt ID" = "expt_id", "Complex" = "complex", "Replicate" = "replicate", "Condition" = "condition"), 
                                                            selected = "gene_symbol"),
                                                selectInput("pointShape", label = ("Point shape"), 
                                                            choices = list("dot" = 16, "square" = 15, "plus" = 3, "cross" = 4, "diamond" = 18,
                                                                           "Assign to gene symbol" = "gene_symbol", "Assign to ID" = "id", 
                                                                           "Assign to expt ID" = "expt_id", "Assign to complex" = "complex", "Assign to replicate" = "replicate", "Assign to condition" = "condition"), 
                                                            selected = 16),
                                                sliderInput("pointSize", "Point size", 0, 10, 1),
                                                
                                                selectInput("lineType", label = ("Line type"), 
                                                            choices = list("solid" = 1, "dashed" = 2, "dotted" = 3, 
                                                                           "Assign to gene symbol" = "gene_symbol", "Assign to ID" = "id", 
                                                                           "Assign to complex" = "complex", "Assign to replicate" = "replicate", "Assign to condition" = "condition"), 
                                                            # "Assign to expt ID" = "expt_id",
                                                            selected = 1),
                                                sliderInput("lineSize", "Line size", 0, 10, 1),
                                                
                                                textInput("plot_title", label = "Plot title", value = ""),
                                                actionButton("colorButton", "Change colors"),
                                                downloadButton('downloadPlot2', "Save plot")
                                                # shinySaveButton('save', 'Save plot', 'Save plot as...', filetype=list(pdf=c('pdf'), png=c('png')))
                                              ), # ends wellPanel
                                              
                                              # main panel contains the reactively rendered output, dependent on selected widget elements
                                              mainPanel(
                                                uiOutput('retrieve_protein_button'),
                                                div(dataTableOutput("prot_table"), style = "font-size:80%"),
                                                div(dataTableOutput("complex_table"), style = "font-size:80%"),
                                                tags$head(tags$style(type="text/css", "#prot_table table.dataTable tr {height:10%;}")),  
                                                tags$head(tags$style(type="text/css", "#prot_table table.dataTable tr.even { background-color: white; }")),  
                                                tags$head(tags$style(type="text/css", "#prot_table table.dataTable tr.odd { background-color: white; }")),  
                                                tags$head(tags$style(type="text/css", "#prot_table table.dataTable tr.odd { background-color: white; }")),  
                                                
                                                #  div(dataTableOutput("complex_table"), style = "font-size:80%"),
                                                
                                                #  hr(),
                                                plotOutput("plot_proteasome",
                                                           #                                        height = 800,
                                                           #                                        width=1000,
                                                           #height = "auto",
                                                           height = 600,
                                                           
                                                           width="100%",
                                                           dblclick = "proteasomePlot_dblclick",
                                                           brush = brushOpts(
                                                             id = "proteasomePlot_brush",
                                                             resetOnNew = TRUE)),
                                                hr()
                                              )  # ends mainPanel
                                            )  # ends sidebarLayout
                                   ),
                                   
                                   tabPanel("Spike-in quantification",
                                            sidebarLayout(
                                              # Sidebar with a slider input
                                              wellPanel(
                                                tags$style(type="text/css", '#leftPanel { width:250px; float:left;}'),
                                                id = "leftPanel",
                                                helpText("Select the relevant variable conditions"),
                                                selectizeInput('show_expt_id_std', 'Select expt ID', choices = NULL,  selected =  NULL, multiple=TRUE),  
                                                selectizeInput('show_trypsin_symbol', 'Select UniProt ID', choices = NULL,  selected =  NULL, multiple=TRUE),  
                                                selectInput("y_axis_choices_tryp", label = ("y variable"), 
                                                            choices = list("Raw intensity" = "raw.intensity",
                                                                           #"LFQ intensity" = "LFQ.intensity", 
                                                                           "MS count" = "MS.MS.count", 
                                                                           "Peptides count"= "peptides.count", 
                                                                           "Unique peptides only" ="unique.peptides.only",
                                                                           "Razor and unique peptides" = "razor.and.unique.peptides"#,
                                                                           #"Sequence coverage" ="sequence.coverage"
                                                            ), 
                                                            selected = "raw.intensity") ,
                                                checkboxGroupInput("y_axis_log_tryp",
                                                                   label = "y transformation",
                                                                   choices = c("log2", "normalize across fractions")),
                                                selectInput("split_by_std", label = ("Split into rows by"), 
                                                            choices = list("None" = "none", "UniProt ID" = "id", "Expt ID" = "expt_id"), 
                                                            selected = "expt_id"),
                                                selectInput("split_by_std_col", label = ("Split into columns by"), 
                                                            choices = list("None" = "none", "UniProt ID" = "id", "Expt ID" = "expt_id"), 
                                                            selected = "none"),
                                                selectInput("color_by_choices_std", label = ("Color by"), 
                                                            choices = list("UniProt ID" = "id", "Expt ID" = "expt_id"), 
                                                            selected = "id"),
                                                
                                                textInput("plot_tryp_title", label = "Plot title", value = ""),
                                                actionButton("colorButton_tryp", "Change colors"), 
                                                downloadButton('downloadPlot3', "Save plot")
                                                #  shinySaveButton('save', 'Save plot', 'Save plot as...', filetype=list(pdf=c('pdf'), png=c('png')))
                                                
                                              ),
                                              # Show a plot of the generated distribution
                                              mainPanel(
                                                uiOutput('retrieve_protein_std_button'),
                                                div(dataTableOutput("std_prot_table"), style = "font-size:80%"),
                                                tags$head(tags$style(type="text/css", "#std_prot_table table td {line-height:50%;}")),  
                                                tags$head(tags$style(type="text/css", "#std_prot_table table.dataTable tr.even { background-color: white; }")),  
                                                tags$head(tags$style(type="text/css", "#std_prot_table table.dataTable tr.odd { background-color: white; }")),
                                                #plotlyOutput("plot_trypsin"),
                                                plotOutput("plot_trypsin",
                                                           #  height = "auto",
                                                           height = 600,
                                                           width="100%",
                                                           dblclick = "proteasomeTrypsinPlot_dblclick",
                                                           brush = brushOpts(
                                                             id = "proteasomeTrypsinPlot_brush",
                                                             resetOnNew = TRUE)),
                                                hr()
                                                
                                                #   div(tableOutput("tryp_table"), style = "font-size:80%"),
                                                #    tags$head(tags$style(type="text/css", "#tryp_table table td {line-height:50%;}"))
                                              )
                                            ))
                                   
                                   
                       )), # ends plot tabPanel
              
              tabPanel("Tables", # shows complete metadata and data tables
                       tabsetPanel(
                         id = 'datasets',
                         tabPanel('Summary plots',downloadLink('save_sum_table', "Export"), h3('  '), dataTableOutput('summarydata')),
                         tabPanel('Individual proteins',downloadLink('save_indiv_table', "Export"), h3('  '), dataTableOutput('proteasomedata')),
                         tabPanel('Spike-in quantification',downloadLink('save_std_table', "Export"),  h3('  '), dataTableOutput('standardsedata'))
                       ) # ends tabsetPanel
              ), # ends tables tabPanel
              
              tabPanel("Correlations",
                       sidebarLayout(
                         # Sidebar with a slider input
                         wellPanel(
                           tags$style(type="text/css", '#leftPanelHeatmap { width:400px; float:left;}'),
                           id = "leftPanelHeatmap",
                           helpText("Select/indicate the relevant variable conditions.\nNote that the entries for Condition and Replicate will\ndetermine the labels printed in the heatmap."),
                           tags$style(type='text/css', "#heatmap.add { width:100%; margin-top: 0px;}"),
                           
                           fluidRow(column(2, offset = 0, style='padding:0px;',  textInput("heatmap.cond0", label = "Cond.*", value = "")),
                                    column(2, offset = 0, style='padding:0px;',  numericInput("heatmap.rep0", min = 1, label = "Repl.*", value = NULL)),
                                    column(6, offset = 0, style='padding:0px;', selectizeInput('show_expt_id_heatmap0', 'Select Expt. ID*', choices = NULL,  selected =  NULL, multiple=FALSE)),
                                    column(1, offset = 0, style='padding:0px; margin-top:30px;margin-left:25px;',   actionLink("heatmap.add", "add"))
                           ),
                           #                                                numericInput("num", label = h3("Numeric input"), value = 1)
                           tags$div(id = 'placeholderHeatmapButtons'),
                           fluidRow(column(1, offset = 0, style='padding:0px; margin-top:0px;margin-left:350px;',   actionLink("heatmap.delete", "delete"))),
                           
                           selectInput("heatmap_measurement", width='200px', label = ("Measurement"), 
                                       choices = list( "Raw" = "raw", "Smoothed" = "superSmu", "Fraction normalized (smoothed)" = "norm.value"), 
                                       selected = "superSmu"),
                           selectInput("heatmap_cor_method",  width='200px', label = ("Correlation method"), 
                                       choices = list("Pearson" = "pearson", "Kendall" = "kendall", "Spearman"  = "spearman"), 
                                       selected = "pearson"),
                           tags$head(tags$style(HTML('#plot_heatmap{background-color: #B8FFBC}')) ),
                           #actionButton("saveButton", "Save") )
                           actionButton("plot_heatmap", "Calculate correlations and plot heatmap"),
                           br(),br(),
                           helpText("Specify the range for the correlation values; this will limit the number of proteins shown in the heatmap."),
                           fluidRow(column(6, 
                                           style='padding:0px; margin-top:0px;margin-left:25px;', sliderInput("heatmap.pairwise.corr.range", "Corr. range for pairwise comp.:",min = -1, max = 1, step=0.01, value = c(-1,0.8)),
                                           sliderInput("heatmap.conditionwise.corr.range", "Corr. range for conditionwise comp.:",min = -1, max = 1, step=0.01, value = c(-1,0.8)))
                           ),
                           helpText("Select a subset of proteins from the top or bottom of the first set of heatmaps which will be shown in more detail below."),
                           fluidRow(
                             column(3, style='padding:0px; margin-top:0px;margin-left:25px;', radioButtons(inputId="heat.n.choice", label="", 
                                                                                                           choices=c("top","bottom"))),
                             column(2, style='padding:0px; margin-top:0px;margin-left:25px;', numericInput("heatmap.n.genes", min = 1, label = "", value = 10))
                           )
                         ),
                         # Show a plot of the generated distribution
                         absolutePanel(left = '450px', top = '100px', right = 0, bottom = 0,
                                       # all
                                       bsAlert("heatmapAlert"),
                                       fluidRow(id = "heatmap_fluid",
                                                splitLayout(cellWidths = c("50%", "50%"),  h4("Pairwise"),  h4("Conditionwise")),
                                                splitLayout(cellWidths = c("50%", "50%"),  h6("ALL"),  h6("ALL")),
                                                
                                                splitLayout(cellWidths = c("50%", "50%"), 
                                                            plotOutput("heatmap_pairwise_all",
                                                                       height = 600,
                                                                       width="100%"#,
                                                            ),                                                              
                                                            plotOutput("heatmap_conditionwise_all",
                                                                       height = 600,
                                                                       width="100%"#,
                                                            )
                                                )
                                       ),
                                       # top or bottom
                                       fluidRow(id = "heatmap_sub_fluid", 
                                                splitLayout(cellWidths = c("50%", "50%"),  h6("SUBSET"),  h6("SUBSET")),
                                                splitLayout(cellWidths = c("50%", "50%"), 
                                                            
                                                            plotOutput("heatmap_pairwise_subset",
                                                                       height = 600,
                                                                       width="100%"#,
                                                            ),                                                              
                                                            plotOutput("heatmap_conditionwise_subset",
                                                                       height = 600,
                                                                       width="100%"#,
                                                            )
                                                )
                                       ),
                                       
                                       fluidRow(
                                         splitLayout(cellWidths = c("50%", "50%"), 
                                                     tableOutput("top_n_pair"),                                                              
                                                     tableOutput("top_n_conditionwise")
                                         )
                                       )
                                       
                         )
                       )),
              
              
              tabPanel("Database", # shows complete metadata and data tables
                       tabsetPanel(
                         id = 'databses',
                         tabPanel('DB Viewer', 
                                  selectizeInput('show_expt_id_db_viewer', 'Select expt ID to view', choices = NULL,  selected =  NULL, multiple=FALSE), 
                                  br(), br(),
                                  htmlOutput("dv_viewer_summary")
                         ),
                         
                         tabPanel("DB Editor", # shows complete metadata and data tables
                                  selectizeInput('show_expt_id_db_browser', 'Select expt ID to edit', choices = NULL,  selected =  NULL, multiple=FALSE),  
                                  tags$head(
                                    tags$style(HTML('#deleteButton{background-color: #E77471}'))
                                  ),
                                  bsButton("EditButton", "Edit"), bsButton("deleteButton", "Delete"),
                                  bsAlert("deleteAlert"),
                                  br(), br(), 
                                  tags$div(id = 'placeholderEditButtons')
                         )           
                       ) # ends tabsetPanel
              )# 
              
              
  )
)

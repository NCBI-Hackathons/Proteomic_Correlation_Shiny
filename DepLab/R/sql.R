#' Reading MaxQuant data
#' 
#' @description 
#' This wrapper functions reads in the MaxQuant output,
#'  \emph{retrieves the proteins of interest}, and generates a long 
#'  \code{data.frame} for \emph{all types of measurement} that are 
#' available across all fractions (e.g., raw intensity, MS count).
#' 
#' @details This function is a wrapper for \code{\link{reading_MQ}},
#' \code{\link{cleaning_MQ}}, and  \code{\link{MQ_to_longFormat}}.
#' The function automatically focusses on a \emph{subset} of proteins,
#' e.g. 
#' 
#' The \emph{types of measurement} that are kept depend on the setting of 
#' \code{MQaccessor_to_plotname}.
#' 
#' @param filename path to MaxQuant output (\code{proteinGroups.txt})
#' @param expt.id unique identifier for the experiment to be read in
#' @param data.subset a string specifying what kind of data should be retrieved,
#' one of:
#' \itemize{
#' \item "poi": extract information for one or more proteins of interest -- 
#' if you choose this option, make sure to indicate the organism! If you want
#' to retrieve the bulk information from MaxQuant output, this is the right
#' parameter.
#' \item "trypsin": extract information for trypsin
#' \item a specific UniProt ID or yeast ORF, e.g. "P00761"
#' }
#' @param organism This determines how the proteins of interest should be
#' retrieved, i.e., indicate either 'human' or 'yeast', or NULL if data.subset 
#' is not 'poi'.
#' @param return.dt Boolean indicating whether a data.table shold be returned
#' instead of a data.frame. Default: FALSE
#' @param data.types List of value types that should be part of the final
#' data set that will be stored in the data base. Check \code{\link{MQaccessor_to_plotname}}
#' for possible options.
#' Default: list("peptides.count", "unique.peptides.only", 
#' "razor.and.unique.peptides","raw.intensity","MS.MS.count")
#' @param db.order.poi Vector indicating the \emph{order} of the columns for the
#' \code{data.frame} holding the values for the proteins of interest. This should 
#' match the expectations of \code{\link{initialize.database}} and
#' \code{\link{add.expt.to.database}}. Default:
#' c("gene_symbol", "fraction", "value", "measurement", "expt_id", "organism")
#' @param db.order.std Vector indicating the _order_ of the columns for the
#' data.frame holding the values for the spiked-in proteins. This should 
#' match the expectations of \code{\link{initialize.datbase}} and
#' \code{\link{add.expt.to.database}}. Default: c("id", "fraction",
#' "value", "measurement", "expt_id")
#'  
#' @return A long data.frame where each row corresponds to a single protein,
#' and \emph{one} value and type of measurement per fraction.
#' 
#' @seealso @seealso \code{\link{reading_MQ}}, \code{\link{cleaning_MQ}},
#' \code{\link{extract_proteinID}}
#' @export
read.MQ.data <- function(filename, expt.id, data.subset = "poi", organism = NULL, 
                         return.dt = FALSE, 
                         data.types = list("peptides.count", "unique.peptides.only", 
                                           "razor.and.unique.peptides", #"sequence.coverage",
                                           "raw.intensity", #"LFQ.intensity",
                                           "MS.MS.count"),
                         db.order.poi = c("gene_symbol","fraction", "value","measurement","expt_id","organism"),
                         db.order.std = c("id", "fraction", "value", "measurement", "expt_id")
) {
  
  #### 1. read MaxQuant output
  mq <- reading_MQ(filename)
  
  #### 2. subset the MaxQuant output and get info from protein ID column, 
  ####    either for proteins of interest for a given organism
  ####    or by specifying specific IDs, e.g., for the spike-in control(s)
  
  if (all(data.subset == "poi")) {# --> a) proteins of interest based on a defined organism
    
    if(is.null(organism)){
      stop("If you want to retrieve proteins of interest, 
           you must specify the organism, e.g., 'yeast'.")
    }
    if(is.null(db.order.poi)){
      stop("If you want to retrieve proteins of interest, 
           you must specify the order of the columns that is expected by the
           functions initializing and extending the data base. See ?read.MQ.data")
    }
    
    mq.clean <- cleaning_MQ(mq, remove.contaminants = TRUE, 
                            remove.decoys = TRUE, poi = organism)
    mq.ids <- extract_proteinID(mq.clean$Protein.IDs, routine = organism)
    
    }else{# --> b) lists of manually defined proteins (or "trypsin")
      if(!is.null(organism)){
        warning(paste("The organism will be ignored for retrieving data for", data.subset))
      }
      
      mq.clean <- cleaning_MQ(mq, remove.contaminants = FALSE,
                              remove.decoys = TRUE, spikeIn = data.subset)
      
      if(all(data.subset == "trypsin")){
        mq.ids <- extract_proteinID(mq.clean$Protein.IDs, routine = "trypsin")
      }else{
        if(length(data.subset > 1)){
          uniprot.regex <- c( paste(".*(", paste(data.subset, collapse = "|"), ").*", sep = ""), "\\1")
        }else{
          uniprot.regex <- c(paste(".*(", data.subset, ").*", sep = ""), "\\1")
        }
        mq.ids <- extract_proteinID(mq.clean$Protein.ID, regex = uniprot.regex,
                                    label = "id", routine = NULL)
      }
    }
  
  #### 3. generate skinny data.frames by iterating through the available
  ####    types of data
  types <- data_types
  
  mq.data <- lapply(types, function(x) MQ_to_longFormat(mq.clean, y = x, return.dt = TRUE, mq.ids))
  mq.data <- rbindlist(mq.data)
  setkey(mq.data, id)
  
  ##### 4. make format suitable for data base entry
  mq.data$expt_id <- as.factor(expt.id)
  mq.data$measurement <- as.factor(mq.data$measurement)
  
  if(all(data.subset == "poi")){
    mq.data$organism <- organism
    
    # this returns a data.table with 'id' (keyed) and 'gene_symbol'
    IDs <- read_IDinfo(organism)
    IDs$organism <- NULL
    
    mq.data <- IDs[mq.data]
    
    # check the validity of the merge
    if( dim( mq.data[is.na(gene_symbol)] )[1] > 0 ){
      diff <- length( unique( mq.data[is.na(gene_symbol)]$id ))
      ex <- head(unique(mq.data[is.na(gene_symbol)]$id))
      warning( paste( diff, "IDs of the MaxQuant output were not found in the general UniProt-based ID list. For example:", paste(ex, collapse = "; "),". Note that these will not be included in the data base.") )
    }
    
    mq.data <- mq.data[!is.na(gene_symbol)][, db.order.poi, with = FALSE]
    
    
  }else{
    colnames(mq.data)[colnames(mq.data) == "uniprot"] <- "id" # legacy code
    
    # make sure that the columns and their order match what is expected for the DB
    mq.data <- mq.data[, db.order.std, with = FALSE]
    #mq.data <- mq.data[, c("id", "fraction", "value", "measurement", "expt_id"), with = FALSE]
  }
  
  if( !return.dt ){
    mq.data <- data.frame(mq.data)
  }
  
  return(mq.data)
    }


#' Create the data base
#' 
#' @details This function initializes an \emph{empty} data base with the 
#' following tables:
#' \itemize{
#' \item id_info
#' \item expt_id: \code{ PRIMARY KEY (expt_id, organism)}
#' \item prot_data
#' \item frac_data: values for each fraction, 
#' \code{PRIMARY KEY (organism, expt_id, gene_symbol, measurement, fraction)}
#' \item std_data: values for the spike-in proteins, 
#' \code{PRIMARY KEY (expt_id, id, measurement, fraction)}
#' }
#' 
#' In addition to tables for the actual values, there are also tables for
#' the different types of possible meta-data that can be input via the shiny
#' app start page:
#' \itemize{
#' \item origin_data: info about the experiment, e.g. experimenter, genotype,
#' cell_type etc.
#' \item prefractionation_method_data
#' \item ms_method_data
#' \item data_processing_data
#' }
#' 
#' @param database.name path to data.base
#' @param organism.name name of the organism, e.g. "yeast" or "human"
#' @param force Boolean indicating whether a \emph{new} data base is going to
#' be made even if a file with the name \code{database.name} already exists.
#' Default: FALSE (= no overwriting)
#' 
#' @export
initialize.database <- function(database.name, organism.name, force = FALSE) {
  db <- dbConnect(SQLite(), dbname=database.name, cache_size = 5000)
  
  # If we are forcing initialization, we overwrite any previously existing data
  if (force) {
    dbGetQuery(conn = db, "drop table if exists id_info")
    dbGetQuery(conn = db, "drop table if exists expt_info")
    dbGetQuery(conn = db, "drop table if exists prot_data")
    dbGetQuery(conn = db, "drop table if exists frac_data")
    dbGetQuery(conn = db, "drop table if exists std_data")
    dbGetQuery(conn = db, "drop table if exists origin_data")
    dbGetQuery(conn = db, "drop table if exists prefractionation_method_data")
    dbGetQuery(conn = db, "drop table if exists ms_method_data")
    dbGetQuery(conn = db, "drop table if exists data_processing_data")
    
  }
  
  dbWriteTable(conn = db, name = "id_info",
               value = read_IDinfo(organism = organism.name),
               row.names = FALSE, append = TRUE, overwrite=FALSE)
  
  # Create table schemas for expt_info, frac_data, and std_data tables
  
  expt.info.sql <- "create table IF NOT EXISTS expt_info (expt_id TEXT NOT NULL, 
  organism TEXT NOT NULL,
  PRIMARY KEY (expt_id, organism))"
  dbGetQuery(conn = db, expt.info.sql)
  
  frac.data.sql <- "create table IF NOT EXISTS frac_data (gene_symbol TEXT NOT NULL,
  fraction INTEGER NOT NULL,
  value REAL NOT NULL,
  measurement TEXT NOT NULL,
  expt_id TEXT NOT NULL,
  organism TEXT NOT NULL,
  PRIMARY KEY (organism, expt_id, gene_symbol, measurement, fraction))"
  dbGetQuery(conn = db, gsub("  ", "", frac.data.sql))
  
  std.data.sql <- "create table IF NOT EXISTS std_data (id TEXT NOT NULL, 
  fraction INTEGER NOT NULL,
  value REAL NOT NULL,
  measurement TEXT NOT NULL,
  expt_id TEXT NOT NULL,
  PRIMARY KEY (expt_id, id, measurement, fraction))"
  dbGetQuery(conn = db, gsub("  ", "", std.data.sql))
  
  origin.data.sql <- "create table IF NOT EXISTS origin_data (expt_id TEXT NOT NULL,
  experimenter TEXT NOT NULL,
  genotype TEXT NOT NULL,               
  cell_type TEXT NOT NULL,
  harvest_date TEXT NOT NULL,
  buffer_composition TEXT NOT NULL,
  lysis_method TEXT NOT NULL,
  digestion_enzyme TEXT NOT NULL,
  notes TEXT,
  PRIMARY KEY (expt_id))"
  dbGetQuery(conn = db, gsub("  ", "", origin.data.sql))
  
  prefractionation.data.sql <- "create table IF NOT EXISTS prefractionation_method_data (expt_id TEXT NOT NULL,
  column_id TEXT NOT NULL,
  amount_protein_loaded REAL NOT NULL,               
  sample_vol_loaded REAL NOT NULL,
  lc_flow_rate REAL NOT NULL,
  lc_fraction_size REAL NOT NULL,
  time_per_fraction REAL NOT NULL,
  fractions_collected REAL NOT NULL,
  PRIMARY KEY (expt_id))"
  dbGetQuery(conn = db, gsub("  ", "",  prefractionation.data.sql))
  
  msmethod.data.sql <- "create table IF NOT EXISTS ms_method_data (expt_id TEXT NOT NULL,
  instrument_id TEXT NOT NULL,
  run_date TEXT NOT NULL,               
  method_length REAL NOT NULL,
  PRIMARY KEY (expt_id))"
  dbGetQuery(conn = db, gsub("  ", "",  msmethod.data.sql))
  
  dataproc.data.sql <- "create table IF NOT EXISTS data_processing_data (expt_id TEXT NOT NULL,
  processing_platform TEXT NOT NULL,
  search_algorithm TEXT NOT NULL,               
  filtering_algorithm TEXT NOT NULL,
  filtering_stringency TEXT NOT NULL,
  PRIMARY KEY (expt_id))"
  dbGetQuery(conn = db, gsub("  ", "",  dataproc.data.sql))
  
  dbDisconnect(db)
}

#' Add uploaded data to the existing data base
#' 
#' Adds data from a new MaxQuant output file to the data base.
#' 
#' @param database.name path to the existing data base, e.g., proteomics.db
#' @param expt.info metadata associated with the experiment
#' (one row per expt), e.g. \code{data.frame(expt_id = input$expt.id, organism = input$organism)}
#' @param prot.data data per protein per expt 
#' (incl all columns that are not frac.data)
#' @param frac.data data per protein per fraction per expt (here, mq.data)
#' @param std.data data per standard per fraction per expt (from mq.data)
#' @param origin.data required metadata about the specific experiment
#' @param prefractionation.data optional metatdata about the specific experiment, can be NULL
#' @param msmethod.data optional metatdata about the specific experiment, can be NULL
#' @param dataproc.data optional metatdata about the specific experiment, can be NULL
#' 
#' @examples
#' \dontrun{
#' ## create the database, and completely overwrite if it already exists (useful for debugging!)
#' initialize.database(database.name, organism.name = orga, force = TRUE)
#' 
#' ## read the data in and turn it into a format suitable for the data base
#' x <- read.MQ.data(filename, expt.id, data.subset = "poi", organism = orga)
#' head(x)
#' ## testing the manual input of protein IDs ( aka standards aka spike-ins)
#' y <- read.MQ.data(filename, expt.id, organism = NULL, data.subset = c("P07477", "Q0140","YAL003W"))
#' head(y)
#' 
#' # create some meta-data
#' origin_df <- data.frame(expt_id = expt.id, experimenter = "myself",
#'                      genotype = "unknown", cell_type = "yeast_cells",
#'                      harvest_date = "Nov 2016", buffer_composition = "TrisHCl", 
#'                      lysis_method = "standard", digestion_enzyme = "Trypsin", 
#'                      notes = NA)
#'  msmethods_df <- data.frame(expt_id = expt.id, instrument_id = "X0000",
#'                         run_date = "Oct 2016", method_length = 1)
#'  dataproc_df <- data.frame(expt_id = expt.id, processing_platform = "unknown", 
#'                       search_algorithm = "unknown", filtering_algorithm = "unknown",
#'                       filtering_stringency = "unknown")
#'  prefractionation_df <- data.frame(expt_id = expt.id, column_id = "x",
#'                                  amount_protein_loaded = 1, sample_vol_loaded = 1, 
#'                                  lc_flow_rate = 1, lc_fraction_size = 1, 
#'                                  time_per_fraction = 1, fractions_collected = 1)
#'
#' # add the data to the db
#' add.expt.to.database(database.name,
#'                      expt.info = data.frame(expt_id = expt.id, organism = "yeast"), 
#'                      prot.data = NULL, frac.data = x, std.data = y,
#'                       origin.data = origin_df, 
#'                      prefractionation.data = prefractionation_df, 
#'                      msmethod.data = msmethods_df, dataproc.data = dataproc_df)
#' }
#' @export
add.expt.to.database <- function(database.name, expt.info, prot.data, frac.data, std.data, origin.data, prefractionation.data, msmethod.data, dataproc.data) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  
  dbGetQuery(conn = db, "DROP INDEX IF EXISTS gene_symbol_idInfo")
  dbGetQuery(conn = db, "DROP INDEX IF EXISTS exptGeneMeas_fracData")
  dbGetQuery(conn = db, "DROP INDEX IF EXISTS organismGeneMeas_fracData")
  dbGetQuery(conn = db, "DROP INDEX IF EXISTS expt_id_exptID")
  
  dbWriteTable(conn = db, name = "expt_info", value = expt.info, row.names = FALSE, append = TRUE)
  dbWriteTable(conn = db, name = "frac_data", value = frac.data, row.names = FALSE, append = TRUE)
  dbWriteTable(conn = db, name = "std_data", value = std.data, row.names = FALSE, append = TRUE)
  dbWriteTable(conn = db, name = "origin_data", value = origin.data, row.names = FALSE, append = TRUE)
  
  if(!is.null(prefractionation.data)){
    dbWriteTable(conn = db, name = "prefractionation_method_data", value = prefractionation.data, row.names = FALSE, append = TRUE)
  }
  
  if(!is.null(msmethod.data)){
    dbWriteTable(conn = db, name = "ms_method_data", value = msmethod.data, row.names = FALSE, append = TRUE)
  }
  
  if(!is.null(dataproc.data)){
    dbWriteTable(conn = db, name = "data_processing_data", value = dataproc.data, row.names = FALSE, append = TRUE)
  }
  
  dbGetQuery(conn = db, "CREATE INDEX IF NOT EXISTS gene_symbol_idInfo ON id_info (gene_symbol)")
  dbGetQuery(conn = db, "CREATE INDEX IF NOT EXISTS exptGeneMeas_fracData ON frac_data (expt_id, gene_symbol, measurement)")
  dbGetQuery(conn = db, "CREATE INDEX IF NOT EXISTS organismGeneMeas_fracData ON frac_data (organism, gene_symbol, measurement)")
  dbGetQuery(conn = db, "CREATE INDEX IF NOT EXISTS expt_id_exptID ON expt_info (expt_id)")
  
  dbDisconnect(db)
}

#' Query the frac_data table I
#' 
#' Filtering by expt_id, gene_symbol and measurement.
query.measurements <- function(database.name, expt.ids, gene.symbols, measurement.type) {
  conn <- src_sqlite(path = database.name)
  dq <-  inner_join(tbl(conn, "frac_data"),
                    tbl(conn, "id_info"),
                    by = "gene_symbol") %>% 
    filter(expt_id %in% c('foo_dummy', expt.ids)) %>% # [bug in dplyr] hack to allow expt.ids vector to only contain one element
    filter(gene_symbol %in% c('foo_dummy', gene.symbols)) %>% 
    filter(measurement == measurement.type)
  # dq$query shows the SQL query that will be sent to the database
  dq %>% collect()
}


#' Query the frac_data table I (v2)
#' 
#' @description Filtering by expt_id, gene_symbol and measurement.
#' This version uses native sql instead of dplyr code.
#' 
#' @param database.name Path to data base.
#' @param expt.ids One or more experiment IDs for which the values should be retrieved.
#' @param gene.symbols The specific, unique gene symbols for which the values should be retrieved.
#' @param measurement.type One of c("peptides.count","unique.peptides.only","razor.and.unique.peptides", "sequence.coverage","raw.intensity").
#'
query.measurements.v2 <- function(database.name, expt.ids, gene.symbols, measurement.type) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  gene <- paste(paste("(", paste(paste("frac_data.gene_symbol=", paste(paste("'", gene.symbols, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  expt_id <- paste(paste("(", paste(paste("frac_data.expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  measurement <- paste(paste("(", paste(paste("frac_data.measurement=", paste(paste("'", measurement.type, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste("SELECT frac_data.gene_symbol, frac_data.fraction, frac_data.value, frac_data.measurement, frac_data.expt_id, id_info.organism, id_info.id FROM frac_data, id_info WHERE frac_data.gene_symbol = id_info.gene_symbol", paste( paste(gene, expt_id, sep=" AND "), measurement, sep = " AND "), sep=" AND ")
  dbGetQuery(conn = db, sqlcmd)
}


#' Convert ID to gene symbol
#' 
#' Filtering by id and organism
return.gene.symbols.from.id <- function(database.name, id.name, organism.name) {
  conn <- src_sqlite(path = database.name)
  dq <- conn %>% 
    tbl("id_info") %>%
    filter(organism == organism.name) %>%
    filter(id %in% c('foo_dummy', id.name))
  dq %>% collect()
}

#' Delete experiments from all tables
#' 
#' Delete by expt_id
#' 
delete.expt<- function(database.name, expt.ids) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  
  organism.name <- DepLab:::list.organism.by.expt.v2(database.name, expt.ids)$organism
  
  expt_id <- paste(paste("(", paste(paste("expt_info.expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste("DELETE FROM expt_info WHERE ", expt_id)
  dbGetQuery(conn = db, sqlcmd)
  
  expt_id <- paste(paste("(", paste(paste("frac_data.expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste("DELETE FROM frac_data WHERE ", expt_id)
  dbGetQuery(conn = db, sqlcmd)
  
  expt_id <- paste(paste("(", paste(paste("std_data.expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste("DELETE FROM std_data WHERE ", expt_id)
  dbGetQuery(conn = db, sqlcmd)
  
  #update gene list
  dbGetQuery(conn = db, paste("DROP TABLE IF EXISTS", organism.name))
  organism <- paste(paste("(", paste(paste("organism=", paste(paste("'", organism.name, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste(paste("CREATE TABLE IF NOT EXISTS", organism.name), "AS SELECT DISTINCT frac_data.gene_symbol FROM frac_data WHERE ", organism)
  dbGetQuery(conn = db, sqlcmd)
}



#' Edit experiment meta-data
#' 
#' Edit experiment meta-data, e.g. organism.
#' 
update.expt<- function(database.name, expt.ids, new.expt.id, new.organism.name) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  
  org <- new.organism.name
  new.expt.id <- paste(paste("'", new.expt.id, sep=""), "'", sep="")
  new.organism.name <- paste(paste("'", new.organism.name, sep=""), "'", sep="")
  
  expt_id <- paste(paste("(", paste(paste("expt_info.expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste(paste(paste( paste(paste("UPDATE expt_info SET expt_id = ", new.expt.id, sep=""), ", organism ="), new.organism.name), " WHERE"), expt_id)      
  dbGetQuery(conn = db, sqlcmd)
  
  expt_id <- paste(paste("(", paste(paste("frac_data.expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste(paste(paste( paste(paste("UPDATE frac_data SET expt_id = ", new.expt.id, sep=""), ", organism ="), new.organism.name), " WHERE"), expt_id)      
  dbGetQuery(conn = db, sqlcmd)
  
  expt_id <- paste(paste("(", paste(paste("std_data.expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste(paste(paste("UPDATE std_data SET expt_id = ", new.expt.id, sep=""), "WHERE"), expt_id)  
  dbGetQuery(conn = db, sqlcmd)
  
  expt_id <- paste(paste("(", paste(paste("origin_data.expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste(paste(paste("UPDATE origin_data SET expt_id = ", new.expt.id, sep=""), "WHERE"), expt_id)  
  dbGetQuery(conn = db, sqlcmd)
  
  expt_id <- paste(paste("(", paste(paste("data_processing_data.expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste(paste(paste("UPDATE data_processing_data SET expt_id = ", new.expt.id, sep=""), "WHERE"), expt_id)  
  dbGetQuery(conn = db, sqlcmd)
  
  expt_id <- paste(paste("(", paste(paste("ms_method_data.expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste(paste(paste("UPDATE ms_method_data SET expt_id = ", new.expt.id, sep=""), "WHERE"), expt_id)  
  dbGetQuery(conn = db, sqlcmd)
  
  expt_id <- paste(paste("(", paste(paste("prefractionation_method_data.expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste(paste(paste("UPDATE prefractionation_method_data SET expt_id = ", new.expt.id, sep=""), "WHERE"), expt_id)  
  dbGetQuery(conn = db, sqlcmd)
  
  #update gene list
  dbGetQuery(conn = db, paste("DROP TABLE IF EXISTS", org))
  organism <- paste(paste("(", paste(paste("organism=", paste(paste("'", org, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste(paste("CREATE TABLE IF NOT EXISTS", org), "AS SELECT DISTINCT frac_data.gene_symbol FROM frac_data WHERE ", organism)
  dbGetQuery(conn = db, sqlcmd)
}


#' Query the frac_data table
#' 
#' Filtering by expt_id, and measurement
query.measurements.by.expt.v2 <- function(database.name, expt.ids, measurement.type) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  expt_id <- paste(paste("(", paste(paste("frac_data.expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  measurement <- paste(paste("(", paste(paste("frac_data.measurement=", paste(paste("'", measurement.type, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste("SELECT fraction, value, expt_id FROM frac_data WHERE ", paste(expt_id, measurement, sep=" AND "))
  dbGetQuery(conn = db, sqlcmd)
}

#' Query the std_data table
#' 
#' Filtering by expt_id, gene_symbol and measurement.
#' 
query.std.measurements.v2 <- function(database.name, expt.ids, std.ids, measurement.type) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  expt_id <- paste("(std_data.expt_id=",  paste(paste("'", expt.ids, sep=""), "'", sep=""), sep="")
  std_id <- paste("std_data.id=",  paste(paste("'", std.ids, sep=""), "')", sep=""), sep="")
  measurement <- paste(paste("(", paste(paste("std_data.measurement=", paste(paste("'", measurement.type, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  expt_std_join <- paste(expt_id, std_id, sep=" AND ", collapse=" OR ")
  sqlcmd <- paste("SELECT * FROM std_data WHERE ", paste(measurement, paste("(", paste(expt_std_join, ")", sep=""), sep=""),  sep=" AND "))
  dbGetQuery(conn = db, sqlcmd)
}




#' Query the frac_data table
#' 
#' Query the data for fract data and provide gene symbols in output
#' 
query.measurements.by.expt.with.gene.symbol.v2 <- function(database.name, expt.ids, measurement.type) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  expt_id <- paste(paste("(", paste(paste("frac_data.expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  measurement <- paste(paste("(", paste(paste("frac_data.measurement=", paste(paste("'", measurement.type, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste("SELECT gene_symbol, fraction, value, expt_id FROM frac_data WHERE ", paste(expt_id, measurement, sep=" AND "))
  dbGetQuery(conn = db, sqlcmd)
}

#' Query the std_data table
#' 
#' Query the data for the spike-in while filtering by expt_id and measurement.
#' 
query.std.measurements.by.expt <- function(database.name, expt.ids, measurement.type) {
  conn <- src_sqlite(path = database.name)
  dq <- tbl(conn, "std_data") %>% 
    filter(expt_id %in% c('foo_dummy', expt.ids)) %>%
    filter(measurement == measurement.type)
  # dq$query shows the SQL query that will be sent to the database
  dq %>% collect()
}


list.expt.ids <- function(database.name) {
  conn <- src_sqlite(path = database.name)
  dq <- conn %>% 
    tbl("expt_info") %>% 
    select(expt_id)
  dq %>% collect() %>% distinct()
}

list.expt.ids.v2 <- function(database.name) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  dbGetQuery(conn = db, " SELECT expt_info.expt_id  FROM expt_info")
}


list.expt.ids.w.organism <- function(database.name) {
  conn <- src_sqlite(path = database.name)
  dq <- conn %>% 
    tbl("expt_info") %>% 
    select(expt_id, organism)
  dq %>% collect() %>% distinct()
}

list.expt.ids.w.organism.v2 <- function(database.name) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  dbGetQuery(conn = db, " SELECT * FROM expt_info")
}




list.organism.by.expt <- function(database.name, expt.ids) {
  conn <- src_sqlite(path = database.name)
  dq <- tbl(conn, "expt_info") %>% 
    filter(expt_id %in% c('foo_dummy', expt.ids)) # [bug in dplyr] hack to allow expt.ids vector to only contain one element
  dq %>% collect()
}


list.organism.by.expt.v2 <- function(database.name, expt.ids) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  expt_id <- paste(paste("(", paste(paste("expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste("SELECT organism FROM expt_info WHERE ", expt_id)
  dbGetQuery(conn = db, sqlcmd)
}


list.origin.data.by.expt.v2 <- function(database.name, expt.ids) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  expt_id <- paste(paste("(", paste(paste("expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste("SELECT * FROM origin_data WHERE ", expt_id)
  dbGetQuery(conn = db, sqlcmd)
}


list.prefractionation.data.by.expt.v2 <- function(database.name, expt.ids) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  expt_id <- paste(paste("(", paste(paste("expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste("SELECT * FROM prefractionation_method_data WHERE ", expt_id)
  dbGetQuery(conn = db, sqlcmd)
}

list.msmethod.data.by.expt.v2 <- function(database.name, expt.ids) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  expt_id <- paste(paste("(", paste(paste("expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste("SELECT * FROM ms_method_data WHERE ", expt_id)
  dbGetQuery(conn = db, sqlcmd)
}

list.data.processing.data.by.expt.v2 <- function(database.name, expt.ids) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  expt_id <- paste(paste("(", paste(paste("expt_id=", paste(paste("'", expt.ids, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste("SELECT * FROM data_processing_data WHERE ", expt_id)
  dbGetQuery(conn = db, sqlcmd)
}


list.expt.by.organism <- function(database.name, organism.name) {
  conn <- src_sqlite(path = database.name)
  dq <- tbl(conn, "expt_info") %>% 
    filter(organism %in% c('foo_dummy', organism.name)) # [bug in dplyr] hack to allow expt.ids vector to only contain one element
  dq %>% collect()
}


list.expt.by.organism.v2 <- function(database.name, organism.name) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  organism <- paste(paste("(", paste(paste("organism=", paste(paste("'", organism.name, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
  sqlcmd <- paste("SELECT expt_id FROM expt_info WHERE ", organism)
  dbGetQuery(conn = db, sqlcmd)
}


list.measurement.types <- function(database.name) {
  conn <- src_sqlite(path = database.name)
  dq <- conn %>% 
    tbl("frac_data") %>%
    select(measurement)
  dq %>% collect() %>% distinct()
}


list.all.gene.symbols <- function(database.name, organism.name) {
  conn <- src_sqlite(path = database.name)
  dq <- conn %>% 
    tbl("frac_data") %>%
    filter(organism == organism.name) %>%
    select(gene_symbol)
  dq %>% collect() %>% distinct()
}


list.all.gene.symbols.v2 <- function(database.name, organism.name, cur_len) {
  db <- dbConnect(SQLite(), dbname = database.name, cache_size = 5000)
  len_now <- as.numeric(dbGetQuery(conn = db, "SELECT MAX(_ROWID_) FROM frac_data LIMIT 1;")[1])
  if ( (len_now == cur_len) || (is.na(len_now)) ){
    organism <- paste(paste("(", paste(paste("organism=", paste(paste("'", organism.name, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
    sqlcmd <- paste(paste("CREATE TABLE IF NOT EXISTS", organism.name), "AS SELECT DISTINCT frac_data.gene_symbol FROM frac_data WHERE ", organism)
    dbGetQuery(conn = db, sqlcmd)
    sqlcmd <- paste("SELECT gene_symbol FROM", organism.name)
    dbGetQuery(conn = db, sqlcmd)
  } else {
    dbGetQuery(conn = db, paste("DROP TABLE IF EXISTS", organism.name))
    organism <- paste(paste("(", paste(paste("organism=", paste(paste("'", organism.name, sep=""), "'", sep=""), sep=""), collapse=" OR "), sep=""), ")", sep="")
    sqlcmd <- paste(paste("CREATE TABLE IF NOT EXISTS", organism.name), "AS SELECT DISTINCT frac_data.gene_symbol FROM frac_data WHERE ", organism)
    dbGetQuery(conn = db, sqlcmd)
    sqlcmd <- paste("SELECT gene_symbol FROM", organism.name)
    dbGetQuery(conn = db, sqlcmd)
  }
}


list.all.gene.symbols.in.id.table <- function(database.name, organism.name) {
  conn <- src_sqlite(path = database.name)
  dq <- conn %>% 
    tbl("id_info") %>%
    filter(organism == organism.name) %>%
    select(gene_symbol)
  dq %>% collect() %>% distinct()
}

list.std.gene.symbols <- function(database.name) {
  conn <- src_sqlite(path = database.name)
  dq <- conn %>% 
    tbl("std_data") %>%
    select(id)
  dq %>% collect() %>% distinct()
}


list.selected.std.gene.symbols <- function(database.name, expt.ids) {
  conn <- src_sqlite(path = database.name)
  dq <- conn %>% 
    tbl("std_data") %>%
    filter(expt_id %in% c('foo_dummy', expt.ids)) %>%
    select(id)
  dq %>% collect() %>% distinct()
}


list.selected.gene.symbols <- function(database.name, expt.ids) {
  conn <- src_sqlite(path = database.name)
  dq <- conn %>% 
    tbl("frac_data") %>% 
    filter(expt_id %in% c('foo_dummy', expt.ids)) %>%
    select(gene_symbol)
  dq %>% collect() %>% distinct()
}

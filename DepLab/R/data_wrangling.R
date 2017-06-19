#' Reading in MaxQuant output
#' 
#' @description 
#' This function expects the MaxQuant output "proteinGroups*.txt"
#' and returns a data frame.
#' 
#' @param filename path to MaxQuant output file
#' 
#' @export
reading_MQ <- function(filename){
  
  mq.in <- read.table(filename, header=TRUE, sep="\t",
                   strip.white = TRUE, fill = TRUE, stringsAsFactors = FALSE,
                   comment.char = "") # important to also capture cases with a # in the Fasta header
  
  return(mq.in)
  
}


#' Check whether a supplied ID conforms with the nomenclatures we're using
#' 
#' @param ID The UniProt or yeast gene ID to be checked.
#' @return Returns a Boolean that indicates whether the supplied ID meets the 
#' criteria of either the UniProt or the yeast gene nomenclature.
check_nomenclature <- function(ID = spikeIn){
  
  isUniProt <- grepl("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", ID)
  isYeast <- grepl("^[YQ]+", ID)
  correctID = isUniProt | isYeast
  
  return(correctID)
}


#' Cleaning MaxQuant output
#' 
#' @description 
#' Before generating a long data.frame, it is helpful to clean the MaxQuant 
#' output, i.e., to remove contaminants, decoy hits and other hits that do not
#' follow the nomenclature of annotated proteins.
#' There is a hierarchy to the parameters, i.e., contaminants are removed first,
#' followed by the removal of decoy hits etc.
#' 
#' @param mq.df Data.frame that contains the raw MaxQuant results.
#' (see \code{\link{reading_MQ}})
#' @param remove.contaminants Boolean to indicate whether rows with "CONT"
#' should be removed (should be set to FALSE if non-yeast proteins are of 
#' interest, including spike-ins)
#' @param remove.decoys Boolean to indicate whether rows with "REV" 
#' should be removed  (should always be set to TRUE)
#' @param poi Specify the organism for which the proteins of interest should be
#' retrieved; currently that can be:
#' \itemize{
#' \item "human": only those hits where the  protein ID conforms with the 
#' TrEMBL nomenclature
#' \item "yeast": only those hits where the protein ID conforms with the
#'  yeast standard nomenclature of either starting with Y or Q
#'  }
#'  Note that `poi` is mutually exclusive with `spikeIn`, i.e., either one
#'  should be set to NULL.
#' @param spikeIn Either "trypsin" or a string with a single UniProtID of
#' interest, e.g. "P00761". Choosing "trypsin" will retrieve values for trypsin 
#' from pig (UniProt ID P00761). This option is mutually exclusive with `poi`, 
#' which should be set to NULL.
#' 
#' @return A subset of the original MaxQuant data.frame.
#' 
#' @seealso \code{\link{MQ_to_longFormat}}, \code{\link{extract_proteinID}}
#'
#' @export
cleaning_MQ <- function(mq.df, remove.contaminants = TRUE,
                        remove.decoys = TRUE, 
                        poi = NULL, spikeIn = NULL){
  # Currently, this function will subset the data.frame more and more, thus
  # multiple filtering options may clash. E.g., if the data.frame is already 
  # filtered to only contain trypsin-related entries, it will most likely not
  # find anything related to a Uniprot search for non-trypsin proteins.
  if(dim(mq.df)[1] == 0)(warning("The input to cleaning_MQ is empty."))
  
  mq.out <- mq.df

  if(!remove.contaminants && !remove.decoys && is.null(poi) && is.null(spikeIn)){
    warning("Note that none of the offered filtering options is set. 
            The in-going data frame should be the same as the out-going one.")
    }
  
  if(remove.contaminants){
    mq.out <- subset(mq.out, !grepl("CON", mq.out$Protein.IDs))
  }
  
  if(remove.decoys){
    mq.out <- subset(mq.out, !grepl("REV", mq.out$Protein.IDs))
  }
  
  if(!is.null(poi)){
    
    if(poi == "yeast"){
      
      mq.out <- subset(mq.out, grepl("^[YQ]+", Protein.IDs))
      
    }else if(poi == "human"){
      
      # the massive regex in the middle is from TrEMBL (http://www.uniprot.org/help/accession_numbers)
      mq.out <- subset(mq.out, grepl("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})", Protein.IDs))
      
      }else{stop("If you would like to retrieve the proteins of interest for the
                 organism for which the experiment was done, specify one of the
                 available options: 'yeast' or 'human'.")
        }
  }
  
  # extracting spiked-in proteins
  if(!is.null(spikeIn)){
    
    if(!is.null(poi)){
      stop("The option to retrieve spike-in entries is not compatible with
           retrieving yeast or human entries. Set `poi = NULL`.")
    }
    
    if(remove.contaminants){
      warning("Extracting the results for spike-ins while at the same time 
              removing contaminants will probably not yield the desired results
              (if anything). Recommended settings: remove.contaminants = FALSE,
              remove.decoys = TRUE, poi = NULL")
    }
    
    if(all(spikeIn == "trypsin")){
      
      mq.out <- subset(mq.out, grepl("P00761$|P00761[^a-zA-Z0-9]", Protein.IDs) &
                           grepl("CON", Protein.IDs))
      }else{
        ID.check <- check_nomenclature(spikeIn)
        if(!all(ID.check)){stop("The ID(s) you supplied to `spikeIn =` do(es) 
            not meet the UniProt or yeast gene nomenclature criteria.")}
        reg.1 <- paste(paste(spikeIn, "$", sep = ""), collapse="|")
        reg.2 <- paste(paste(spikeIn, "[^a-zA-Z0-9]", sep = ""), collapse="|")
        reg.combi <- paste(reg.1, reg.2, sep = "|", collapse = "")
        mq.out <- subset(mq.out, grepl( reg.combi, Protein.IDs))
        }
    }
  
  # done cleaning
  if(dim(mq.out)[1] == 0){
    warning("None of the entries in the MaxQuant output survived the cleaning. 
            Check that you selected the correct organism for the data that you uploaded.")
    }
  
  return(mq.out)
}

#' Match the MaxQuant accessor to the plot name
#' 
#' @description
#' This is a helper function that matches the y-axis value name that will
#' be displayed in the fraction profile plot with the regex needed 
#' for extracting the data from the MaxQuant output.  
#' 
#' @details
#' This function is needed because the MaxQuant output has a separate column 
#' for each value per fraction, but we need a long data.frame where the fraction
#' is indicated in an extra column.
#' This function helps to a) retrieve all the columns from the MaxQuant output
#' related to a column and b) to keep a meaningful name for the value that can
#' be used for the axis labeling in the plot later on.
#'  
#' @param plotname type of value to be used for y-axis values in the fraction
#' profile plot
#' 
#' @return Regular expression for extracting all columns associated with the 
#' desired value of interest from the MaxQuant output.
#' 
#' @export
MQaccessor_to_plotname <- function(plotname){
  
  # the data frame contains all column headers from the MaxQuant output
  # that are present for each fraction (except Identification.type since its
  # values are strings instead of numbers)
  matching <- data.frame(mq_filter = c("Peptides\\.[0-9]+",
                                      "Unique\\.peptides\\.[0-9]+",
                                      "Razor.*\\.[0-9]+", # "Sequence\\.coverage\\.[0-9]+",
                                      "Intensity.*\\.[0-9]+",#  "LFQ\\.intensity.*\\.[0-9]+",
                                      "MS\\.MS.*\\.[0-9]+"),
                        ylabel = c("peptides.count",
                               "unique.peptides.only",
                               "razor.and.unique.peptides", # "sequence.coverage",
                               "raw.intensity",# "LFQ.intensity",
                               "MS.MS.count"),
                        stringsAsFactors = FALSE)
  
  mq2plot.regex <- subset(matching, ylabel == plotname)[,1]
  return(mq2plot.regex)
}

#' Pre-populate info table with all the known protein/gene names
#'
#' @description To initiate the id.info data base, this function retrieves the
#' information stored in the ID_lists folder of the DepLab package.
#' The lists should have two columns: 1) unique IDs, e.g. UniProt identifiers;
#' 2) easily recognizable gene names, e.g. the primary gene name from UniProt.
#' These will be pasted, so that the gene name column will also be unique.
#' 
#' @param organism string indicating the organism for which to retrieve the information.
#' 
#' @return Data.table with two columns, 'id' and 'gene_symbol'; 'id' is keyed.
#' Note that 'gene_symbol' is a bit of a misnomer as it's a mash up of the ID 
#' and the primary gene name.
#' 
#' @details
#' \bold{Yeast data} derived from
#' \url{http://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab}
#' and generated using an awk routine.
#' 
#' \bold{Human} info downloaded from uniprot.org with the following search:
#' organism:human AND organism:"Homo sapiens (Human) [9606]" 
#' AND proteome:up000005640
#' 
#' cut -f 1-2 *tab | sed "1d" > human_uniprot.entry.name
#' 
#' @seealso \code{\link{make_gene_name}}
#' @export
read_IDinfo <- function(organism) {
  id.info <- read.table( system.file("extdata", "ID_lists", paste(organism, "id.txt", sep="."), package = "DepLab"),
                         header = FALSE, stringsAsFactors = FALSE,
                         col.names = c("id", "gene_symbol"), sep = "\t")
  # get a unique "gene symbol"
  id.info <- make_gene_name(id.info)
  id.info <- data.table(id.info)
  setkey(id.info, id)
  id.info$organism <- organism
  
  return(id.info)
}

#' Generate unique, yet easily recognizable gene names
#'
#' @description To have a compromise for a unique protein/gene identifier and
#' an easily recognizable name, we're merging the unique UniProt ID with the
#' non-unique gene name (the one that is indicated as "primary" gene name by
#' UniProt).
#' If the gene_symbol string is empty, it will be replaced with the UniProt ID.
#' 
#' @param id.info data.frame with 2 columns: "id" and "gene_symbol"
#' @return data.frame with the exact same 2 columns, except that gene_symbol will
#' now contain: "UniProtID (gene name)"
#' @seealso \code{\link{read_IDinfo}}
make_gene_name <- function(id.info){
  
  if( length(unlist(stringr::str_extract_all(names(id.info), "^id$|^gene_symbol$"))) != 2){
    stop(paste("The data.frame supplied to make_gene_name should contain the 
              columns 'id' and 'gene_symbol'. Instead, it currently has the following ones:",
               names(id.info) ) )
  }
  
  id.info$gene_symbol <- with(id.info,
                            ifelse( grepl(".+", gene_symbol),
                                    paste0(id, " (", gene_symbol, ")" ), id) )
  
  if( length(unique(id.info$gene_symbol)) != length(id.info$gene_symbol) ){
    stop("The gene symbol entry modified by make_gene_names() is still not unique.
         Check the ID column of id.info.")
  }
  return(id.info)
}

#' Read in lists of known complexes
#' 
#' @description This function reads in a table where for each protein a possible
#' complex membership is noted.
#' 
#' @param filename Path to file. Should have one row per protein with at least:
#' 'id', 'complex'
#' @param organism Organism for which the complex list was generated, e.g. 'human'.
#' 
#' @return Data.frame with 'id', 'gene_symbol' and 'complex' where the gene symbol
#' corresponds to the gene symbol format used for all the raw data
#' ( < UniProtID (gene name) > ).
#' @export
read_complexes <- function(filename, organism){
  c.in <- data.table(read.table(filename, sep="\t", stringsAsFactors=FALSE, header=TRUE, quote="\"") )
  
  if( length(unlist(stringr::str_extract_all(names(c.in), "^id$|^gene_symbol$"))) != 2){
    stop(paste(filename, "should contain the headers 'id' and 'gene_symbol'.
               Instead, it currently has the following ones:", names(c.in) ) )
  }
  
  ID.check <- check_nomenclature(c.in$id)
  if( !all(ID.check) ){
    stop(paste("Not all ID(s) in", filename, "match the UniProt or yeast gene nomenclature."))
  }
  
  setkey(c.in, id)
  
  ids <- read_IDinfo(organism)
  
  c.merge <- ids[c.in][,c("id", "gene_symbol","complex"), with = FALSE] 
  
  return(data.frame(c.merge))
}

#' Extract ID information
#' 
#' @description
#' This function currently contains all the regex we've come up with to extract
#' IDs of interest as defined by the Dephoure lab.
#' 
#' @details
#' To parse the protein-centric information into the data base using a consistent
#' ID scheme, these information must be extracted from the sometimes erratic
#' MaxQuant output.
#' Older MaxQuant formats have the information in $Fasta.headers, newer ones will
#' have a $Gene.names column and an empty $Fasta.headers column. 
#' We will therefore rely on the $Protein.IDs only and retrieve the gene_symbol
#' info from our external data base.
#' 
#' @param prot.id Vector of protein.id entries or the name of the column 
#' with related information
#' @param routine Specify the type of IDs you would like to retrieve.
#' Currently, the following options are available:
#'  \itemize{
#'  \item "trypsin": retrieving the protein ID for pig trypsin
#'  \item "yeast": retrieves the protein IDs etc. that conform with the official 
#'  yeast nomenclature (either starting with Y or Q)
#'  \item "human": retrieves protein IDs that conform with the TrEMBL format 
#'  and mnemonic protein names of UniProt
#'  }
#' @param regex If you do not want to use one of the pre-established routines,
#' you should specify a **vector of two strings**, where the first string 
#' corresponds to the pattern and the second string to the replacement.
#' For example, regex =  c("F", "f") will replace
#' upper case F with lower case f within the protein ID column.
#' @param label If you don't want to use a routine, indicate a label for your 
#' customized filtering, e.g. "gene_IDs"
#' 
#' @return If a routine is chosen, this returns a list of ID vectors. If regex
#' and label are used, the list will contain only one vector of IDs.
#' Both results can directly be used with
#' \code{\link{MQaccessor_to_plotname()}}, which is used within
#' \code{\link{MQ_to_longFormat()}}.
#' 
#' @seealso \code{\link{MQ_to_longFormat}}
#' 
#' @export
extract_proteinID <- function(prot.id, routine = NULL, regex = NULL, label = NULL){

  if(!is.character(prot.id)){
    stop("Please supply a vector of characters for the protein IDs.")
  }
  out <- list()
  
  if( is.null(routine)){
    if( is.null(regex) | is.null(label)){
      stop("If you choose not to run a routine for extracting information from 
           the Protein.ID column, you must specify the regular expression of your 
           choice and a label for the type of information that you would like 
           to extract." )
    }
    
    out[[label]] <- gsub( regex[1], regex[2], prot.id)
    
  }else{
    if(routine == "trypsin"){
      out$id <- gsub(".*(P00761).*","\\1", prot.id)
      # out$id <- gsub(">([A-Z0-9.]+) .*", "\\1", prot.id)
      # out$id_type <- gsub(">([A-Z0-9a-z.]+) ([A-Za-z0-9_-]+):.+", "\\2", fasta.info)
    }
  
    if(routine == "yeast"){

      out$id <- gsub("^([YQ][0-9A-Z-]*)(;.*)*", "\\1", prot.id)
      #out$gene_symbol <- gsub("^>[YQ][0-9A-Z-]*[ ]([,A-Z0-9-]*)[ ].*", "\\1", fasta.info)
      #out$sgdid <- gsub(".*SGDID:([A-Z0-9a-z]*),.*", "\\1", fasta.info)
    }
    
    if(routine == "human"){
      out$id <- gsub("(^sp\\|)*([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(\\|.*)*(;.*)*","\\2", prot.id)
    }
    
    if(!unique(unlist( lapply(out,
                              function(x) grepl("^[A-Z0-9]", x) ) 
    ))) {
      warning("Double-check that the info in the supplied column matches the 
              expectation for the retrieval of human or yeast proteins.")
    }
  }
  
  return(out)
}


#' Generate a long data.frame for a specific MS value
#'
#' @description 
#' One problem with the MaxQuant output is that there is a separate column
#' for each value/fraction pair.
#' This function will generate a long data.frame where each row contains the
#' value of interest per fraction.
#'
#' @param mq.df data.frame based on MaxQuant output without 
#' decoy hits (REV) and contaminants (see also \cite{\link{cleaning_MQ}}
#' @param y type of value to be used for y-axis values
#' @param return.dt Boolean indicating whether a data.table or data.frame should be returned.
#' Default: FALSE
#' @param ... several vectors (or a list thereof), e.g. of protein IDs,
#' that will be merged with the actual values. These should therefore have the
#' same length as the data.frame!
#' 
#' @seealso \code{\link{reading_MQ}}, \code{\link{MQaccessor_to_plotname}}, 
#' \code{link{extract_proteinID}}, \code{\link{read.mq.data}}
#' 
#' @export
MQ_to_longFormat <- function(mq.df, y, return.dt = FALSE, ...){
  
  if( !is.list(...)){
    vars.cbind <- list(...)
  }else{ # is there a more elegant way to do this if we supply a list as the last argument?
    vars.cbind <- list(...) # makes a list of the list; list(...) is the only way I've found to gather the input for (...)
    vars.cbind <- vars.cbind[[1]]
  }
  
  if( !(dim(mq.df)[1] == unique( unlist( lapply(vars.cbind, length) )) ) ){
    stop("The MaxQuant data.frame and the additionally supplied IDs or other 
         information do not have the same number of observations." )
    }
  
  y_regex <- MQaccessor_to_plotname(y)
  mq.long <- mq.df[grepl(y_regex, names(mq.df))]
  
  # check the types of entries in mq.long
  cls <- unique(sapply(mq.long, class))
  
  mq.long <- data.table(cbind(..., mq.long))
  
  # set all values to numeric if some are not
  if( length(cls) > 1){
    start <- 1+length(vars.cbind)
    end <- dim(mq.long)[2]
    mq.long[, start:end := lapply(.SD, as.numeric), .SDcols = start:end ]
    }

  # making sure all columns are numeric
  #start <- 1+length(vars.cbind)
  #end <- dim(mq.long)[2]
  #mq.long[, start:end := lapply(.SD, as.numeric), .SDcols = start:end ]
  
  mq.long <- data.table::melt(mq.long,
                  id.vars = c(1:length(vars.cbind)), # specify the columns with the id variables
                  variable.name = "fraction")
  mq.long$fraction <- as.numeric(gsub(".*\\.([0-9]+)_.*","\\1", mq.long$fraction))
  names(mq.long)[length(names(mq.long))] <- "value"
  mq.long$measurement <- y

  if(!return.dt){
    mq.long <- data.frame(mq.long)
  }
  return(mq.long)
}

#' @title
#' Extract all information from the MaxQuant output that is
#' present only once per protein
#'
#' @description 
#' The MaxQuant table contains values per fraction, which can be extracted via
#' \code{\link{MQ_to_longFormat}}, and entries that are present only once per
#' protein ID. These entries include, but are not limited to,
#' the Protein ID itself, (sometimes) the FASTA headers, numbers of proteins and 
#' peptides detected as well as oxidation and phosphorylation sites.

#' @param mq.df data.frame based on MaxQuant output
#' @return data.frame of protein-specific entries
#'
#' @details
#' \bold{example usage:}
#' 
#' mq <- reading_MQ(filename = "../dephourelab/DepLab/inst/extdata/test_data/proteinGroups_100mM_new.txt")
#' 
#' prot.specific <- extract_protspec(mq)
#' 
#' @export
extract_protspec <- function(mq.df){
  
  ps.headers <- c("Protein.IDs", "Majority.protein.IDs", "Peptide.counts..all.",
                  "Peptide.counts..razor.unique.", "Peptide.counts..unique.",
                  "Fasta.headers", "Number.of.proteins", "Peptides",
                  "Razor...unique.peptides", "Unique.peptides", "Sequence.coverage....",
                  "Unique...razor.sequence.coverage....", "Unique.sequence.coverage....",
                  "Mol..weight..kDa.",  "Sequence.length", "Sequence.lengths", "Q.value", 
                  "Score", "Intensity", "MS.MS.Count", "Only.identified.by.site", "Reverse",
                  "Potential.contaminant", "id", "Peptide.IDs", "Peptide.is.razor", "Mod..peptide.IDs",
                  "Evidence.IDs", "MS.MS.IDs", "Best.MS.MS", "Oxidation..M..site.IDs", "Phospho..STY..site.IDs",
                  "Oxidation..M..site.positions", "Phospho..STY..site.positions")  
  
  ps.data <- mq.df[!grepl("_f", names(mq.df))]
  
  if( !identical( intersect(ps.headers, names(ps.data)),
                  ps.headers) )(warning(cat("The input file 
                                            does not have the following header(s):\n", 
                                            paste( setdiff( ps.headers, names(ps.data) ),
                                                   collapse = ","), 
                                            "\n Header(s) present in the input file 
                                            that we have not seen before:\n", 
                                            paste( setdiff( names(ps.data), ps.headers ),
                                                   collapse = ",") )
                                        ))
  
  return(ps.data)
}


#' Smooth measured values
#' 
#' @description 
#' This function allows for the smoothening of, e.g., raw intensity values
#' collected from the MaxQuant output using Friedman's superSmoother.
#' 
#' @param long.df data.frame with unique values per protein, expt. ID, 
#' measurement and fraction; preferably for only one measurement.
#' @param prot.identifier specify the name that identifies single proteins,
#' e.g. "gene_symbol"
#' 
#' @return data.frame with normalized intensity values
#' 
#' @seealso \code{\link{supsmu}}
#' 
#' @details
#' \bold{example usage:}
#'
#'                
#' @export
superSmooth_values <- function(long.df,  prot.identifier = NULL){
  
  if(length(unique(long.df$measurement)) > 1)(stop("Please supply only one type of measurement."))
  
  dt <- data.table(long.df)
  
  if(is.null(prot.identifier))(stop("You must specify the accessor for the unique 
                                      protein identifiers using 'prot.identifier'."))
    
  check_columns(c("expt_id", "measurement"), long.df, "long.df","superSmooth_values")
  
  setorderv(dt, c("expt_id", prot.identifier,"measurement", "fraction"))
  keys = c(prot.identifier, "expt_id")
  setkeyv(dt, keys) 
  
  dt.sm <- dt[, value := supsmu(x = seq(1, length(value)), y = value )$y, by = key(dt)]
  dt.sm <- dt.sm[, measurement:= paste(measurement, "superSmu", sep = "_")]
  
  return(as.data.frame(dt.sm))
}



#' Normalize measured values
#' 
#' @description 
#' This function allows for the normalization of, e.g., raw intensity values
#' collected from the MaxQuant output.
#' Currently, two normalization schemes are implemented.
#' 
#' @param long.df data.frame with unique values per protein, expt. ID, 
#' measurement and fraction; preferably for only one measurement.
#' @param norm.type choose the type of normalization to carried out; currently,
#' there are two possibilities:
#' \itemize{
#' \item 'fraction': normalize the measurements per protein across all fractions
#'  ("row-wise")
#' \item 'spike-in': normalize the measurements per fraction across all proteins
#'  ("column-wise")
#' }
#' @param prot.identifier specify the name that identifies single proteins,
#' e.g. "gene_symbol"; only necessary if norm.type == "fraction"
#' @param std.df long data.frame with values for the spike-in control with 'value',
#' 'measurement', and 'fraction'; this should only be the data for a
#'  \emph{single} protein, though.
#' 
#' @return data.table with normalized intensity values
#' 
#' @seealso \code{\link{norm_by_fraction}},
#' \code{\link{norm_by_spikeIn}}
#' 
#' @details
#' \bold{example usage:}
#'
#' ## fraction normalization 
#' mq.y.1 <- read.MQ.data(filename = system.file("extdata", "test_data", "proteinGroups_100mM_new.txt", package = "DepLab"), 
#'                        expt.id = "100mM", data.subset = "poi", organism = "yeast")
#' mq.y.3 <- read.MQ.data(filename = system.file("extdata", "test_data", "proteinGroups_300mM_new.txt", package = "DepLab"), 
#'                        expt.id = "300mM", data.subset = "poi", organism = "yeast")
#' mqcombi <- rbind(mq.y.1, mq.y.3)
#' fraction.norm <- normalize_values(long.df = subset(mqcombi, measurement == "raw.intensity"),
#'                                  norm.type = "fraction", prot.identifier = "gene_symbol")
#'
#' ## spike-in normalization
#' # first, read in the values for the spike-in (here: trypsin) 
#' y.std.1 <- read.MQ.data(filename = system.file("extdata", "test_data", "proteinGroups_100mM_new.txt", package = "DepLab"), 
#'                        expt.id = "100mM", data.subset = "trypsin", organism = NULL)
#' y.std.3 <- read.MQ.data(filename = system.file("extdata", "test_data", "proteinGroups_300mM_new.txt", package = "DepLab"),
#'                       expt.id = "300mM", data.subset = "trypsin", organism = NULL)
#' std.combi <- rbind(y.std.1, y.std.3)
#' 
#' # the subsetting for one type of measurement is absolutely required
#' std.norm <- normalize_values(long.df = subset(mqcombi, measurement == "raw.intensity"), 
#'                norm.type = "spike-in", prot.identifier = "gene_symbol", 
#'                std.df = subset(std.combi, measurement == "raw.intensity"))
#'                
#' @export
normalize_values <- function(long.df, norm.type = c("fraction","spike-in"),
                             prot.identifier = NULL,
                             std.df = NULL){

  if(length(unique(long.df$measurement)) > 1)(stop("Please supply only one type of measurement."))
  
  dt <- data.table(long.df)
    
  if(norm.type == "fraction" | norm.type == "superSmu"){
    if(is.null(prot.identifier))(stop("If you want to normalize across fractions,
                                      you must specify the accessor for the unique 
                                      protein identifiers using 'prot.identifier'."))
    
    keys = c(prot.identifier, "expt_id")
    setkeyv(dt, keys) #setkey*v* to allow for strings as keys
    dt.norm <- norm_by_fraction(dt)
    
    
  }else if(norm.type == "spike-in"){
    
    if(is.null(std.df))(stop("If you want to normalize using a spike-in control,
                             you must supply the long data.frame with the 
                             corresponding values for the spike-in using 'std.df'."))
    
    keys = c("fraction", "expt_id")
    setkeyv(dt, keys)
    
    std.dt <- data.table(std.df)
    setkeyv(std.dt, keys)
    dt.norm <- norm_by_spikeIn(dt, std.dt)
  
  }else{(stop("You did not specify an available normalization method."))
  }
  
  return(as.data.frame(dt.norm))
}



#' Normalize values across fractions
#' 
#' @description 
#' This normalization assumes that proteins present at equal levels in each
#' fraction should get a value of 1 in each fraction.
#' ("row-wise normalization")
#' 
#' @param long.dt data.table with one row per protein and fraction; the keys
#' should be set to identify unique proteins and experiment IDs.
#' Fractions should not be key'ed!
#' 
#' @return It returns virtually the same data.table except that the 'value'
#' column will contain normalized values and the 'measurement' column will have
#' been updated to reflect the normalization.
#' 
#' @seealso \code{\link{normalize_values}}
#' 
#' @export
#' @import data.table
norm_by_fraction <- function(long.dt){
  
  if(!("fraction" %in% names(long.dt) && 
         "value" %in% names(long.dt)))(stop("The data.table must have the accessors 'fraction' and 'value' in order to be normalized."))
  
  if(!haskey(long.dt))(stop("The data.table must have at least one key that is unique per protein, e.g. 'gene_name'."))
  
  if("fraction" %in% key(long.dt))(stop("The data.table has a key for fraction. This is incompatible with normalizing the measurements across fractions."))
  
  # in data.table queries, column names are used like variables
  # := is used to update an existing column or add a new column
  dt.fn <- long.dt[, value:= value * length(unique(fraction)) / sum(value), by = key(long.dt)]
  dt.fn <- dt.fn[, measurement:= paste(measurement, "fract.normalized", sep = "_")]
  
  # replace NaNs with 0
  dt.fn[is.na(value)]$value <- 0
  
  return(dt.fn)
}


#' Normalize values based on spike-in
#' 
#' @description 
#' This function will normalize the values of each fraction (for all proteins) to the
#' corresponding value from the spike-in control in that fraction
#' ("column-wise normalization").
#' 
#' @param long.dt data.table with one entry per protein and fraction with values
#' for one (!) type of measurement, e.g. raw intensity; it should have keys for
#' 'fraction', 'measurement', and 'expt.id'
#' @param spike.dt data.table of the same format as long.dt, but with values
#' from the standard protein that was spiked-in in equal amounts into each
#' fraction; it should have keys for 'fraction', 'measurement', and 'expt.id'
#' 
#' @return It returns a data.table with the same structure except that the
#' 'value' column will contain normalized values and the 'measurement' column 
#' will have been updated to reflect the normalization.
#' 
#' @seealso \code{\link{normalize_values}}
#' 
#' @export
#' @import data.table
norm_by_spikeIn <- function(long.dt, spike.dt){
  
  # make sure that there's only one spike-in per experimental id
  exid <- spike.dt[, .N, by = c("expt_id", "id") ]$expt_id
  
  if( any(duplicated( exid )) )(stop( paste("There should be only 1 spiked-in protein per experimental ID for the normalization. Check sample:", exid[ duplicated(exid) ] ))) 
  
  if(! ( haskey(long.dt) && haskey(spike.dt) &&
           grepl("fraction", key(long.dt)) && 
           grepl("fraction", key(spike.dt)) ))(stop("Check that both data.tables have at least a key for 'fraction'."))
  
  if(!("measurement" %in% names(long.dt) && 
         "value" %in% names(long.dt) &&
         "measurement" %in% names(spike.dt) &&
         "value" %in% names(spike.dt)))(stop("The data.tables must have the accessors 'fraction' and 'measurement' in order to be normalized."))
  
  if( length(unique(spike.dt$measurement)) > 1 ){
    stop(paste("The spike-in data table has more than one type of measurement:", paste(unique(spike.dt$measurement), collapse = ", "),  "This is incompatible with the normalization using norm_by_spikeIn()"))
  }
  
  ####################

  # join spike-in and individual proteins
  merged.dt <- long.dt[spike.dt]
  merged.dt <- merged.dt[!is.na(merged.dt$measurement)]
  if(dim(merged.dt)[1] < 1)(warning("After merging the spike-in table with the measurement table, no values are left."))
  
  # do the actual normalization
  merged.dt <- merged.dt[, value:= value/i.value]
  merged.dt <- merged.dt[, measurement:= paste(measurement, "spike-in.normalized", sep= "_")]
  
  # replace NaNs with 0
  merged.dt[is.na(value)]$value <- 0
  
  # return only those columns that the original long.dt has;
  # .SD allows for subsetting the data.table supplying strings
  merged.dt <- merged.dt[,.SD, .SDcols = names(long.dt)]
  setkeyv(merged.dt, key(long.dt))
  
  return(merged.dt)
}



#' Preparing data for plotting of "metadata"
#'
#' @description 
#' "Metadata" here simply refers to the summing up of all measurements
#' in one fraction across all available peptides/proteins.
#' 
#' @param mq.df data frame with MaxQuant data
#' @param value string indicating which sort of value to summarize; one of:
#' c("raw.intensity", "MS.MS.count", "peptides.count",
#' "unique.peptides.only", "razor.and.unique.peptides")
#' @return data.table suitable for plotting with ggplot2
#' @seealso \code{\link{plot_cumulativeValuesAcrossFractions}},
#' \code{\link{MQ_to_ggplot2}}
#' 
#' @export
MQ_to_cumulativePlot <- function(mq.df, value){
  mq.clean <- cleaning_MQ(mq.df, remove.contaminants = TRUE, remove.decoys = TRUE, poi = "yeast")
  df.plot <- MQ_to_longFormat(mq.clean, y = value, extract_proteinID(mq.clean$Protein.ID, routine = "yeast"))
  dt.plot <- data.table(df.plot)
  setkey(dt.plot, fraction)
 # dt.plot <- dt.plot[,.(peptide.count=sum(peptides.count, na.rm = TRUE)), by=fraction] # this must be generalized, not specific for peptides.count!
 dt.plot <- dt.plot[, .(sum( get("value"), na.rm = TRUE)), by = fraction]
 names(dt.plot)[2] <- value
  
  return(dt.plot)
}



#' Retrieve the UniProt information via the protein ID
#' 
#' @description 
#' This function will use the UniProt ID to access UniProt's information
#' about that protein via their API.
#' 
#' @param protein.ID String with the UniProt accessor, e.g. "P00761"
#' 
#' @return A string indicating the protein ID together with the first line of 
#' the information of the .txt file provided by UniProt
#' (e.g. http://www.uniprot.org/uniprot/P07761.txt)
#' plus the link to the corresponding html page, which provides more
#' information and is easier to look at.
#' 
#' @details
#' \bold{example usage:}
#'
#' get_UniProt_info("P00761")
#'
#' @export
get_UniProt_info <- function(protein.ID){
  
  # getURL does not care about the format, it will dump everything it finds into
  # on object that then needs to be turned into a table of sorts
  url.dump <- RCurl::getURL( paste("http://www.uniprot.org/uniprot/", 
                                   protein.ID, ".txt", sep = "") )
  
  dump.df <- read.table(textConnection(url.dump),
                        header = F, sep = "\t", strip.white=TRUE, fill= TRUE,
                        quote="", stringsAsFactors = FALSE)
  
  out.string <- paste(protein.ID, dump.df[1,])
  return(out.string)
}


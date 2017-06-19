#' Wrapper for generating matrices of correlation coefficients
#' 
#' 
#' @param in.dt data.table with one row per protein fraction and different
#' columns for different types of measurements
#' @param measurement the type of measurement that should be used for the correlation
#' calculation, e.g. "raw" or "superSmu"
#' @param corr.method string indicating which method should be used by the base R cor() function.
#' One of "pearson" (default), "kendall", or "spearman"
#' @param uniquifiers vector of strings that define a unique sample,
#' e.g.: c("gene_symbol","condition","replicate")
#' @param calc.mode indicate whether you want to calculate the correlation for each pair of samples ["pairwise"], 
#' or between conditions ["condition"], or both ["both"].
#'
#' @return List with one or two matrices for heatmap visualization.
#' 
#' @seealso \code{\link{cor}}
corrHM_wrap <- function(in.dt, measurement="superSmu", corr.method = "pearson", uniquifiers=uniquifiers, calc.mode = c("both","pairwise","condition")){
  
  check_columns(c("gene_symbol", measurement, uniquifiers), in.dt, "in.dt", "corrHM_wrap")
  
  out <- list()
  
  if(calc.mode == "both" | calc.mode == "pairwise"){
    pw.cors.dt <- lapply(unique(in.dt$gene_symbol), 
                         function(x) pw_corr_wrap(in.dt, x,
                                                  measurement = measurement,
                                                  uniq.factors = c("fraction", uniquifiers),
                                                  cast.form = "fraction ~ condition + replicate",
                                                  method="pearson")) %>% rbindlist
    
    if(sum(dim(pw.cors.dt)) == 0){
      warning(str(in.dt), "\nCorrelations could not be computed for any of the proteins.
           This is most likely due to lack of conditions and/or replicates, i.e., if 
           there's just one value per protein fraction. Check the structure of the input
           data.table above.")
      }else{
        #print("pw.cors.dt")
       #print(head(pw.cors.dt))
        pw.cors.mat <- corrHM_prep(pw.cors.dt)
        #print("pw.cors.mat")
        #print(head(pw.cors.mat))
        out$pairwise <- pw.cors.mat 
      }
  }
  
  if(calc.mode == "both" | calc.mode == "condition"){
    
    # calculate median values for each fraction per gene and condition
    agg.dt <- aggregate_dt(in.dt, measurement = measurement,
                           agg_by = c("fraction", uniquifiers[!grepl("replicate",uniquifiers)]) )
   
    # calculate correlations for the median values
    # again, this will take a while to compute
    # TO ADD: support for mean in addition to median
    agg.cors.dt <- lapply(unique(agg.dt$gene_symbol),
                          function(x) pw_corr_wrap(agg.dt, x,
                                                   measurement = "median.norm.intensity",
                                                   uniq.factors = c("fraction","gene_symbol", "condition"),
                                                   cast.form = "fraction  ~  condition" )) %>% rbindlist
    
    if(sum(dim(agg.cors.dt)) == 0){
      warning(print(str(agg.dt)),"Condition-wise correlations could not be computed for any of the proteins.
           This is most likely due to lack of conditions, i.e., there's just one value per
           protein fraction. Check the structure of the aggregated data.table shown above
           this error message and/or perhaps use calc.mode == 'pairwise'.")
      }else{
        agg.cors.mat <- corrHM_prep(agg.cors.dt)
        out$condition_comp <- agg.cors.mat
      }
  }
  
  return(out)
}


#' Aggregating replicate values
#' 
#' @param in.dt tidy data.table
#' @param measurement string indicating the type of measurment that will be used
#' for the aggregation. must be a column header in in.dt, e.g. "superSmu"
#' @param agg_by strings indicating the combination of column headers that identify
#' unique proteins, e.g., c("condition","replicate","gene_symbol")
#' @return data.table with aggregated values
#' @details 
#' # summed.dt <- aggregate_dt(hmdat, measurement = "superSmu",agg_by = uniquifiers)
  # head(summed.dt)
    # gene_symbol           condition replicate median.norm.intensity
    # 1: A0AV02 (SLC12A8)        EV        02          2.476554e-09
    # 2: A0AV02 (SLC12A8)        WT        01          4.171135e+04
    # 3:  A0AVI2 (FER1L5)     D109N        01          0.000000e+00
    # 4: A0AVI4 (TMEM129)        EV        02          1.524699e-09
    # 5:    A0AVT1 (UBA6)     D109N        01          7.792054e+06
    # 6:    A0AVT1 (UBA6)        EV        01          8.808870e+06
aggregate_dt <- function(in.dt, measurement, agg_by){
  
  check_columns(c(agg_by, measurement), in.dt, "in.dt", "aggregate_dt")
  
  # to do: add switch function here to allow for mean and sd, too
  agg.dt <- in.dt[, median(get(measurement), na.rm = TRUE), by = agg_by]
  
  names(agg.dt)[which(names(agg.dt) == "V1")] <-"median.norm.intensity"
  
  return(agg.dt)
}


#' Preparing tidy correlation data.table for heatmap visualization
#' 
#' @param in.dt data.table in tidy format with:, e.g. output from pw_corr_map()
#' @return matrix fit for most heatmap functions c("gene_symbol","comparison","correlation"),
#' i.e.: proteins per row, samples per columns, 
#' gene symbols as row names. If there's only one column, the matrix will be sorted in descending order.
corrHM_prep <- function(in.dt){
  
  check_columns(c("gene_symbol","comparison","correlation"), in.dt, "in.dt", "corrHM_prep")
  
  # get one column per pw comparison
  dtw <- dcast(in.dt, gene_symbol ~ comparison, value.var = "correlation")
  n.cols <- dim(dtw)[2] - 1
  n.rows <- dim(dtw)[1]
  
  if( min(c(n.cols, n.rows)) == 0){
    stop("The data.table submitted to corrHM_prep does not have enough entries.")
  }
    
  cdt.mat <- as.matrix(dtw[, -c("gene_symbol"), with = FALSE], nrow = n.rows, ncol = n.cols)
  colnames(cdt.mat) <- names(dtw[, -c("gene_symbol"), with = FALSE])
  rownames(cdt.mat) <- dtw$gene_symbol
  
  # remove proteins where more than 1 comparison is missing
  # and sort of there's only one column
  if(n.cols > 1){
    cdt.mat <- cdt.mat[apply(cdt.mat, 1, function(x) length(which(is.na(x)))) <= 1,]
  }else{
    col.nme <- dimnames(cdt.mat)[[2]]
    cdt.mat <- cdt.mat[!is.na(cdt.mat), ] %>% sort(., decreasing = TRUE) %>% as.matrix(., ncol = 1) 
    colnames(cdt.mat) <- col.nme
  }
  
  return(cdt.mat)
  
}


#' Wrapper for calculating pairwise correlations
#' 
#' This function will calcule the pairwise correlation of multiple profiles for 
#' one specific protein or gene.
#' 
#' @param in.dt data.table in tidy format
#' @param gene must match an entry in the gene_symbol column of the in.dt
#' @param mesasurement string indicating the type of value that should be used for
#' the correlation calculation, e.g. "superSmu" or "raw"
#' @param uniq.factors vector of strings that determines how the pairwise
#' correlations will be calculated, e.g. c("gene_symbol","fraction","condition","replicate")
#' @param cast.form formula for casting the pairwise comparisons, in string format,
#'  e.g., "fraction  ~  condition + replicate"
#' @param ... additional parameters for \code{\link{cor}}
#' 
#' @return tidy data.table of pairwise correlations; if no correlation can be
#' determined (e.g. because there are not enough values), nothing will be returned
#' \tabular{llll}{
#' \t comparison \t correlation  \t gene_symbol \cr
#' [,1] \t D109N_01_vs_EV_01   \t 0.7964631    \t A0AVT1 (UBA6) \cr
#' [,2] \t D109N_01_vs_EV_02   \t 0.7659069    \t A0AVT1 (UBA6) \cr
#' [,3] \t D109N_01_vs_WT_01   \t 0.8525967    \t A0AVT1 (UBA6) \cr
#' [,4] \t D109N_01_vs_WT_02   \t 0.7876166    \t A0AVT1 (UBA6) \cr
#' [,5] \t EV_01_vs_EV_02      \t 0.9859125    \t A0AVT1 (UBA6) \cr
#' [,6] \t EV_01_vs_WT_01      \t 0.9897667    \t A0AVT1 (UBA6) \cr
#' [,7] \t EV_01_vs_WT_02   \t 0.9984263    \t A0AVT1 (UBA6) \cr
#' [,8] \t EV_02_vs_WT_01   \t 0.9587254    \t A0AVT1 (UBA6) \cr
#' [,9] \t EV_02_vs_WT_02   \t 0.9866967    \t A0AVT1 (UBA6) \cr
#' [,10] \t WT_01_vs_WT_02   \t 0.9865275    \t A0AVT1 (UBA6) 
#' }
#' 
#' @examples
#'  \dontrun{
#'  # retrieve data
#'  hmdat <- as.data.table(query.measurements.by.expt.with.gene.symbol.v2(DB, exps, "raw.intensity"))
#'  setorder(hmdat, gene_symbol, expt_id,fraction)
#'  setnames(hmdat, "value","raw")
#'  
#'  # calculate pairwise correlation; this takes a while to compute
#'  cors.dt <- lapply(unique(hmdat$gene_symbol), 
#'                    function(x) pw_corr_wrap(hmdat, x, measurement = "superSmu",
#'                                          uniq.factors = c("fraction", uniquifiers),
#'                                           cast.form = "fraction ~ condition + replicate",
#'                                           method="pearson")) %>% rbindlist
#'
#' # obtain matrix for heatmap plotting
#' cors.mat <- corrHM_prep(cors.dt)
#'  }
pw_corr_wrap <- function(in.dt, gene, measurement, 
                         uniq.factors = c("gene_symbol","fraction","condition","replicate"), 
                         cast.form, ...){
  
    check_columns(c("fraction",uniq.factors, measurement), in.dt, "in.dt", "pw_corr_wrap")
    
    v <- in.dt[gene_symbol == gene, c(uniq.factors, measurement), with = FALSE]
    w <- dcast(v, as.formula(cast.form), value.var = measurement)
    
    w <- w[ , which(unlist(lapply(w, function(col) sum(col) > 0 ) ) ), with=F]
    w <- w[, -"fraction", with =FALSE]
    
    if(dim(w)[1] > 1 & dim(w)[2] > 1){
      # generate correlation matrix
      cor.out <- cor(w, ...)
      
      # turn correlation matrix into data table
      cor.dt <- get_cor_dt(cor.out)
      cor.dt$gene_symbol <- gene
      return(cor.dt)
    }

}

#' @title Turn the 2D-correlation matrix into a 1D-data table
#' @param cor.mat matrix of pairwise correlations, obtained via cor()
#' @return a data.table with one column for the pairwise correlations,
#' one column describing the comparison
#' 
#' @examples
#' \dontrun{
#' cors.dt <- lapply(unique(dat.test$id),
#'                   function(x) {
#'                      w <- dcast(dat.test[id == x,
#'                                          c("fraction","expt_id","replicate","superSmu"), 
#'                                          with = FALSE], 
#'                                 fraction  ~  expt_id + replicate, value.var = "superSmu")
#'                      w <- w[,which(unlist(lapply(w, function(col) sum(col) > 0 ) ) ), with=F]
#'                      w <- w[, -"fraction", with =FALSE]
#'                      
#'                      if(dim(w)[1] > 1 & dim(w)[2] > 1){
#'                         cor.out <- cor(w, method  = "pearson")
#'                         cor.dt <- get_cor_dt(cor.out)
#'                         cor.dt$id <- x
#'                         return(cor.dt)
#'                       }
#'          }) %>% rbindlist
#'}
get_cor_dt <- function(cor.mat){
  # get the names of the comparisons
  seen = NULL
  comp = NULL
  one = dimnames(cor.mat)[[1]]
  for(N_i in c(1:(length(one)-1) )){
    N_j <- N_i + 1
    comp <- c(comp, paste(one[N_i], one[N_j:length(one)], sep = "_vs_"))
  }
  
  # get the lower part of the correlation matrix
  corr = cor.mat[lower.tri(cor.mat)]
  
  return(data.table(comparison = comp, correlation = corr))
}

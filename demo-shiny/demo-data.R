#library("tidyr")
library("DepLab")
library("DepLabData")
library("NMF")
library("heatmaply")
library("shinyHeatmaply")
# install.packages(c("heatmaply", "shinyHeatmaply"))
# devtools::install_github("NCBI-Hackathons/Proteomic_Correlation_Shiny/DepLab")
# devtools::install_github("julia-wrobel/DepLabData")

# bundle the logic for setting up the data and the plots in functions
# export global objects for the data and plots at the end of this script
# source this script from the Shiny app.R
# that way, we can tweak the data and plot setup functions and logic while maintaining the same exported objects for the Shiny app. 


# ~~~~~ HEATMAP ~~~~~ #
make_heatmap_data <- function(){
    # data from the data package
    data("WT_trial1")
    
    wt1 <- MQ_to_longFormat(WT_trial1, y = "raw.intensity",
                            return.dt = TRUE,
                            extract_proteinID(WT_trial1$Protein.IDs,
                                              routine = "human"))
    
    # remove entries with all zero values
    wt1_sum <- wt1[, sum(value), id]
    keep <- wt1_sum[V1 > 0]$id
    wt1 <- wt1[id %in% keep]
    
    # make matrix
    wt1_wide <-  dcast(wt1, id ~fraction)
    wt1_mat <- as.matrix(wt1_wide[,-"id", with=FALSE])
    rownames(wt1_mat) <- wt1_wide$id
    
    wt1_mat_subset <- log1p(wt1_mat[1:500,])

    return(wt1_mat_subset) 
}

#make_heatmap <- function(){
    # interactive heatmaply Plotly heatmap
#    heatmaply_data <- make_heatmap_data()
    
#    my_heatmap <- heatmaply(x = heatmaply_data, distfun = function(x) dist((1-cor(t(x), method = "pearson"))), scale = "none", Colv = FALSE, plot_method = "ggplot") 
#    return(my_heatmap)
#}

make_ranked_clustered_data <- function(mat = make_heatmap_data()){
    # calculate the dendrogram clustering distances and add it to the object
    # note: mat = matrix returned by `make_heatmap_data`
    
    # get the distances between entries
    mat_dist <- dist((1-cor(t(mat), method = "pearson")))
    # get the clustering
    mat_clust <- hclust(mat_dist)
    
    # make df with the labels and the rank order
    rank_labels <- data.frame(label = as.character(mat_clust[["labels"]]), order = mat_clust[["order"]], stringsAsFactors = FALSE)
    
    # clean up the original matrix, convert to df, add rownames as column
    mat_df <- as.data.frame(mat)
    mat_df[["label"]] <- rownames(mat_df)
    rank_labels[["label"]]
    mat_merged <- merge(x = mat_df, y = rank_labels, by = "label")
    return(mat_merged)
}

# ~~~~~ QC PLOTS ~~~~~ #
data("WT_trial1")
data("WT_trial2")
data("EV_trial1")
data("EV_trial2")

ids <- extract_proteinID(WT_trial1$Protein.IDs, routine = "human")
WT1 <- lapply(c("raw.intensity", "peptides.count", "sequence.coverage"), function(val_of_interest){
  MQ_to_longFormat(WT_trial1, y= val_of_interest, return.dt = TRUE, 
                   list(id = ids$id, expt_id = WT_trial1$expt_id) 
  )}) %>% rbindlist

ids <- extract_proteinID(WT_trial2$Protein.IDs, routine = "human")
WT2 <- lapply(c("raw.intensity", "peptides.count"), function(val_of_interest){
  MQ_to_longFormat(WT_trial2, y= val_of_interest, return.dt = TRUE, 
                   list(id = ids$id, expt_id = WT_trial2$expt_id) 
  )}) %>% rbindlist

ids <- extract_proteinID(EV_trial1$Protein.IDs, routine = "human")
EV1 <- lapply(c("raw.intensity", "peptides.count"), function(val_of_interest){
  MQ_to_longFormat(EV_trial1, y= val_of_interest, return.dt = TRUE, 
                   list(id = ids$id, expt_id = EV_trial1$expt_id) 
  )}) %>% rbindlist

ids <- extract_proteinID(EV_trial2$Protein.IDs, routine = "human")
EV2 <- lapply(c("raw.intensity", "peptides.count"), function(val_of_interest){
  MQ_to_longFormat(EV_trial2, y= val_of_interest, return.dt = TRUE, 
                   list(id = ids$id, expt_id = EV_trial2$expt_id) 
  )}) %>% rbindlist
test_data <- rbindlist(list(WT1, WT2, EV1, EV2))
rm(WT1)
rm(WT2)
rm(EV1)
rm(EV2)

#' Plot sum of MQ variables per protein
#'
#' @description This function calculates the sum of values across all fractions, whose
#' distributions are then shown via violin plots.
#'
#' @return ggplot object with a violin and boxplot for every \code{expt_id}
#' @param in.dt data.table in long format with one type of measurement
#' in a column named "value"
#' @param uniq.factors vector of column accessors that define unique protein entries
#' (Default: c("gene_symbol","expt_id"))
#' @param value.name name that will be used for display, should be descriptive of the
#' type of measurement that is being summed up 
plot_sums_per_protein <- function(in.dt, uniq.factors = c("gene_symbol", "expt_id"), value.name = NULL){
  
  # summing up whatever is stored in value for
  # a group defined by uniq.factors
  summing <- in.dt[ , sum(value), uniq.factors]
  
  setnames(summing, "V1", "sum")
  
  # make plot
  P <- ggplot(summing, aes(x = expt_id, y = log2(sum))) + geom_violin() + geom_boxplot(notch = TRUE)
  
  # add descriptive name to the column that contains the sum
  if(!is.null(value.name)){
    P <- P + ylab(paste0("log2(",value.name,")"))
  }else{
    P <- P + ylab("log(sum)")
  }
  
  return(P)
}




# ~~~~~ PROFILE PLOT ~~~~~ #
make_profile_plot_data <- function(){
    # copy/paste from other function in case other function changes
    # data from the data package
    data("WT_trial1")

    wt1 <- MQ_to_longFormat(WT_trial1, y = "raw.intensity",
                            return.dt = TRUE,
                            list(id = extract_proteinID(WT_trial1$Protein.IDs, routine = "human")$id,
                                 expt_id = WT_trial1$expt_id)  
                            )
    
    # remove entries with all zero values
    wt1_sum <- wt1[, sum(value), id]
    keep <- wt1_sum[V1 > 0]$id
    wt1 <- wt1[id %in% keep]
    
    return(wt1)
}

make_profile_plot <- function(df, selected_IDs = c("A0FGR8")){
   dt <- data.table(df)
   df.sub <- as.data.frame(dt[id %in% selected_IDs])
    my_profile_plot <- plot_profile(df.sub,
                                    what = c("id", "expt_id"),
                                    color.by = "id", line.smooth = FALSE)
    return(my_profile_plot) 
}




# ~~~~~ OBJECTS TO USE IN THE SHINY ~~~~~ #
# qc plots
p_rawInt <- plot_sums_per_protein(test_data[measurement == "raw.intensity"], 
                      uniq.factors = c("id","expt_id"), 
                      value.name = "raw intensity") +
  theme_bw() + geom_hline(yintercept = 20, linetype = "dashed", color = "gray")

p_nPep <- plot_sums_per_protein(test_data[measurement == "peptides.count"], 
                               uniq.factors = c("id","expt_id"), 
                               value.name = "raw intensity") +
  theme_bw() + geom_hline(yintercept = 1, linetype = "dashed", color = "gray")

# heatmap
#my_heatmap <- make_heatmap()
ranked_clustered_data <- make_ranked_clustered_data()

# xy plots
profile_plot_data <- make_profile_plot_data()
my_profile_plot <- make_profile_plot(df = profile_plot_data)


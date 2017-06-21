library("tidyr")
library("DepLab")
library("DepLabData")
library("NMF")
library("heatmaply")
library("shinyHeatmaply")
# install.packages(c("heatmaply", "shinyHeatmaply"))


# bundle the logic for setting up the data and the plots in functions
# export global objects for the data and plots at the end of this script
# source this script from the Shiny app.R
# that way, we can tweak the data and plot setup functions and logic while maintaining the same exported objects for the Shiny app. 


# ~~~~~ HEATMAP ~~~~~ #
make_heatmap_data <- function(){
    # data from the data package
    data("WT_trial1")
    colnames(WT_trial1)
    
    wt1 <- MQ_to_longFormat(WT_trial1, y = "raw.intensity",
                            return.dt = TRUE,
                            extract_proteinID(WT_trial1$Protein.IDs,
                                              routine = "human"))
    wt1 <- wt1[, -4]
    wt1_wide <- spread(wt1, fraction, value)
    wt1_mat <- as.matrix(wt1_wide[,-1])
    rownames(wt1_mat) <- wt1_wide$id
    wt1_mat <- wt1_mat[rowSums(wt1_mat)>0,]
    
    wt1_mat_subset <- log1p(wt1_mat[1:500,])
    wt1_mat_subset <- as.data.frame(wt1_mat_subset)
    return(wt1_mat_subset) # is a matrix
}

make_heatmap <- function(){
    # interactive heatmaply Plotly heatmap
    heatmaply_data <- make_heatmap_data()
    
    my_heatmap <- heatmaply(x = heatmaply_data, distfun = function(x) dist((1-cor(t(x), method = "pearson"))), scale = "none", Colv = FALSE, plot_method = "ggplot") 
    return(my_heatmap)
}

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




# ~~~~~ PROFILE PLOT ~~~~~ #
make_profile_plot_data <- function(){
    # copy/paste from other function in case other function changes
    # data from the data package
    data("WT_trial1")
    colnames(WT_trial1)
    
    wt1 <- MQ_to_longFormat(WT_trial1, y = "raw.intensity",
                            return.dt = TRUE,
                            extract_proteinID(WT_trial1$Protein.IDs,
                                              routine = "human"))
    wt1 <- wt1[, -4]
    wt1_wide <- spread(wt1, fraction, value)
    wt1_mat <- as.matrix(wt1_wide[,-1])
    rownames(wt1_mat) <- wt1_wide$id
    wt1_mat <- wt1_mat[rowSums(wt1_mat)>0,]
    
    wt1_mat_subset <- log1p(wt1_mat[1:500,])
    wt1_mat_subset <- as.data.frame(wt1_mat_subset)
    
    # need the Experiment ID column for the dataset
    wt1$expt_id <- "WT1"
    return(wt1)
}

make_profile_plot <- function(df, selected_IDs = c("A0FGR8")){
    my_profile_plot <- plot_profile(df[which(as.character(df[["id"]]) %in% selected_IDs) ,], what = c("id", "expt_id"), color.by = "id", line.smooth = FALSE)
    return(my_profile_plot) # is a matrix
}


# ~~~~~ OBJECTS TO USE IN THE SHINY ~~~~~ #
my_heatmap <- make_heatmap()
ranked_clustered_data <- make_ranked_clustered_data()

profile_plot_data <- make_profile_plot_data()
my_profile_plot <- make_profile_plot(df = profile_plot_data)


library(shiny)
library("DepLab")
library("DepLabData")
library("NMF")
library("heatmaply")
library("shinyHeatmaply")
# install.packages(c("heatmaply", "shinyHeatmaply"))


# ~~~~~ HEATMAP ~~~~~ #
# data from the data package
# code from daviderisso
library("tidyr")
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

# old method ofr heatmap
# pdf(file = "./heatmpa.pdf")
# aheatmap(log1p(wt1_mat[1:500,]), Colv = NA, distfun = "pearson",
#          scale = "none")
# dev.off()

# new interactive heatmap
wt1_mat_subset <- log1p(wt1_mat[1:500,])
wt1_mat_subset <- as.data.frame(wt1_mat_subset)
#my_heatmap <- heatmaply(x = wt1_mat_subset, distfun = function(x) dist((1-cor(t(x), method = "pearson"))), scale = "none", Colv = FALSE, plot_method = "ggplot") 

order.outcome = as.matrix(wt1_mat_subset) 
J = dim(order.outcome)[2] ## J is number of observations per subject
N = dim(order.outcome)[1] ## N is number of subjects

rownames(order.outcome) = seq(1:N)
colnames(order.outcome) = seq(1:J)
names(dimnames(order.outcome))= c("subj", "grid")


outcome.df = melt(order.outcome) 

#outcome.df %>% filter(subj == 1) %>%
#ggplot(., aes(grid, value, group = subj)) + geom_point() + geom_path()




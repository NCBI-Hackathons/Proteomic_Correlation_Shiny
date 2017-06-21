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
my_heatmap <- heatmaply(x = wt1_mat_subset, distfun = function(x) dist((1-cor(t(x), method = "pearson"))), scale = "none", Colv = FALSE, plot_method = "ggplot") 


# ~~~~~ PROFILE PLOT ~~~~~ #
# need the Experiment ID column for the dataset
wt1$expt_id <- "WT1"
selected_IDs <- c("A0FGR8")
my_profile_plot <- plot_profile(wt1[which(as.character(wt1[["id"]]) %in% selected_IDs) ,], what = c("id", "expt_id"), color.by = "id", line.smooth = FALSE)


#####Nick's additions#####

# ~~~~~ PCA ~~~~~ #
#####PCA; DN_Trial2 proteins#####
data("DN_trial2")
df<-DN_trial2
protids<-df$Protein.IDs
rownames(df)<-protids
pattern<-"intensity"
colinds<-grep(pattern,colnames(df))
mat<-df[,colinds]
rownames(mat)<-protids
pr_mat<-prcomp(t(mat),
               center = T,
               scale. = F)
#plot(pr_mat, type = "l")
#summary(pr_mat)
library(ggfortify)
id<-rownames(mat)
my_pca<-ggplot(pr_mat$rotation,aes(label=id),hjust=0, vjust=0)+
  geom_point(aes(x=PC1,y=PC2))+
  ggtitle("PCA; DN_Trial2 intensities")
library(plotly)
my_pca<-ggplotly(my_pca,source="my_pca")
#####3D PCA#####
m<-data.frame(pr_mat$rotation)
my_3dpca<-plot_ly(m,x= ~PC1,y= ~PC2,z= ~PC3,hovertext=id) %>% 
  add_markers() %>%
  layout(title = "PCA; DN_Trial2 proteins")
my_3dpca
# #####MDS; DN_Trial2 proteins#####
# data("DN_trial2")
# df<-DN_trial2
# protids<-df$Protein.IDs
# rownames(df)<-protids
# pattern<-"intensity"
# colinds<-grep(pattern,colnames(df))
# mat<-df[,colinds]
# rownames(mat)<-protids
# id<-rownames(mat)
# d<-dist(mat,method="euclidean")
# mds<-cmdscale(d, eig = TRUE)
# m<-mds$points
# colnames(m)<-c("X1","X2")
# my_mds<-ggplot(m,aes(label=id),hjust=0, vjust=0)+
#   geom_point(aes(x=X1,y=X2),color="orangered")+
#   ggtitle("MDS; DN_Trial2 proteins")
# library(plotly)
# ggplotly(my_mds)
# #####3D MDS#####
# mds<-cmdscale(d, k=3,eig = TRUE)
# m<-data.frame(mds$points)
# colnames(m)<-c("X1","X2","X3")
# library(plotly)
# my_3dmds<-plot_ly(m,x= ~X1,y= ~X2,z= ~X3,hovertext=id,color="YlOrRd") %>% 
#   add_markers() %>%
#   layout(title = "MDS; DN_Trial2 proteins")
# my_3dmds
# #####Tsne#####
# library(Rtsne)
# data("DN_trial2")
# df<-DN_trial2
# protids<-df$Protein.IDs
# rownames(df)<-protids
# pattern<-"intensity"
# colinds<-grep(pattern,colnames(df))
# mat<-df[,colinds]
# rownames(mat)<-protids
# id<-rownames(mat)
# library(plotly)
# tsne_out<-Rtsne(mat,dim=2,check_duplicates=F)
# m<-data.frame(tsne_out$Y)
# colnames(m)<-c("tsne_1","tsne_2")
# my_tsne<-ggplot(data=m,aes(label=id))+
#   geom_point(aes(x=tsne_1,y=tsne_2),color="orangered")+
#   ggtitle("tSNE; DN_Trial2 proteins")
# my_tsne
# #####3D Tsne#####
# tsne_out<-Rtsne(mat,dim=3,check_duplicates=F)
# m<-data.frame(tsne_out$Y)
# colnames(m)<-c("tsne_1","tsne_2","tsne_3")
# my_3dtsne<-plot_ly(m,x= ~tsne_1,y= ~tsne_2,z= ~tsne_3,
#                    marker = list(size = 5),hovertext=id,color="YlOrRd") %>% 
#   add_markers() %>%
#   layout(title = "tSNE; DN_Trial2 proteins")
# my_3dtsne
##########


#####load Data#####
load("/Users/nickgiangreco/Downloads/matrices_for_pca.rda")
library(ggfortify)
#####PCA; mean matrices#####
mat<-mean_npeaks
pr_mat<-prcomp(mat,
               center = F,
               scale. = T)
plot(pr_mat, type = "l")
summary(pr_mat)
autoplot(pr_mat,label=F,label.size=5,size=4)+ggtitle("PCA; mean_npeaks")
mat<-t(mean_bestpeaks)
pr_mat<-prcomp(mat,
               center = F,
               scale. = F)
plot(pr_mat, type = "l")
summary(pr_mat)
autoplot(pr_mat,label=T,label.size=5,size=0)+ggtitle("PCA; mean_bestpeaks")
#####MDS; mean matrices#####
mat<-mean_npeaks
d<-dist(mat,method="euclidean")
mds<-cmdscale(d, eig = TRUE)
autoplot(mds,label=F,label.size=5,size=3)+ggtitle("MDS; mean_npeaks")
mat<-t(mean_bestpeaks)
d<-dist(mat,method="euclidean")
mds<-cmdscale(d, eig = TRUE)
autoplot(mds,label=T,label.size=5,size=0)+ggtitle("MDS; mean_bestpeaks")
#####PCA; WT_Trial1#####
data("WT_trial1")
df<-WT_trial1
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
g<-ggplot(pr_mat$rotation,aes(label=id),hjust=0, vjust=0)+
  geom_point(aes(x=PC1,y=PC2))+
  ggtitle("PCA; WT_Trial1 proteins")
library(plotly)
ggplotly(g)
#####PCA; WT_Trial2#####
data("WT_trial2")
df<-WT_trial2
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
g<-ggplot(pr_mat$rotation,aes(label=id),hjust=0, vjust=0)+
  geom_point(aes(x=PC1,y=PC2))+
  ggtitle("PCA; WT_Trial2 proteins")
library(plotly)
ggplotly(g)
#####PCA; EV_Trial1#####
data("EV_trial1")
df<-EV_trial1
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
g<-ggplot(pr_mat$rotation,aes(label=id),hjust=0, vjust=0)+
  geom_point(aes(x=PC1,y=PC2))+
  ggtitle("PCA; EV_Trial1 proteins")
library(plotly)
ggplotly(g)
#####PCA; EV_Trial2#####
data("EV_trial2")
df<-EV_trial2
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
g<-ggplot(pr_mat$rotation,aes(label=id),hjust=0, vjust=0)+
  geom_point(aes(x=PC1,y=PC2))+
  ggtitle("PCA; EV_Trial2 proteins")
library(plotly)
ggplotly(g)
#####PCA; DN_Trial2 intensities#####
data("DN_trial2")
df<-DN_trial2
protids<-df$Protein.IDs
rownames(df)<-protids
pattern<-"intensity"
colinds<-grep(pattern,colnames(df))
mat<-df[,colinds]
rownames(mat)<-protids
pr_mat<-prcomp(mat,
               center = T,
               scale. = F)
#plot(pr_mat, type = "l")
#summary(pr_mat)
library(ggfortify)
int<-colnames(mat)
g<-ggplot(pr_mat$rotation,aes(label=int),hjust=0, vjust=0)+
  geom_point(aes(x=PC1,y=PC2))+
  ggtitle("PCA; DN_Trial2 intensities")
library(plotly)
ggplotly(g)
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
g<-ggplot(pr_mat$rotation,aes(label=id),hjust=0, vjust=0)+
  geom_point(aes(x=PC1,y=PC2))+
  ggtitle("PCA; DN_Trial2 intensities")
library(plotly)
ggplotly(g)
#####3D PCA#####
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
library(plotly)
m<-data.frame(pr_mat$rotation)
p<-plot_ly(m,x= ~PC1,y= ~PC2,z= ~PC3,hovertext=id) %>% 
  add_markers() %>%
  layout(title = "PCA; DN_Trial2 proteins")
p
#####Tsne#####
library(Rtsne)
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
library(plotly)
m<-data.frame(pr_mat$rotation)
tsne_out<-Rtsne(m,check_duplicates=F)
plot(tsne_out$Y)
d<-dist(m,method="euclidean")
tsne_out<-Rtsne(d)
plot(tsne_out$Y)

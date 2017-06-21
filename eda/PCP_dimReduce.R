#####load Data#####
load("/Users/nickgiangreco/Downloads/matrices_for_pca.rda")
library(ggfortify)
#####PCA; mean matrices#####
mat<-mean_bestpeaks
pr_mat<-prcomp(t(mat),
               center = F,
               scale. = T)
plot(pr_mat, type = "l")
summary(pr_mat)
autoplot(pr_mat,label=F,label.size=5,size=4)+ggtitle("PCA; mean_npeaks")
mat<-mean_bestpeaks
pr_mat<-prcomp(mat,
               center = F,
               scale. = F)
plot(pr_mat, type = "l")
summary(pr_mat)
autoplot(pr_mat,label=F,label.size=5,size=0)+ggtitle("PCA; mean_bestpeaks")
#####MDS; mean matrices#####
mat<-mean_bestpeaks
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
#####MDS; DN_Trial2 proteins#####
data("DN_trial2")
df<-DN_trial2
protids<-df$Protein.IDs
rownames(df)<-protids
pattern<-"intensity"
colinds<-grep(pattern,colnames(df))
mat<-df[,colinds]
rownames(mat)<-protids
id<-rownames(mat)
d<-dist(mat,method="euclidean")
mds<-cmdscale(d, eig = TRUE)
m<-mds$points
colnames(m)<-c("X1","X2")
my_mds<-ggplot(m,aes(label=id),hjust=0, vjust=0)+
  geom_point(aes(x=X1,y=X2),color="orangered")+
  ggtitle("MDS; DN_Trial2 proteins")
library(plotly)
ggplotly(my_mds)
#####3D MDS#####
mds<-cmdscale(d, k=3,eig = TRUE)
m<-data.frame(mds$points)
colnames(m)<-c("X1","X2","X3")
library(plotly)
my_3dmds<-plot_ly(m,x= ~X1,y= ~X2,z= ~X3,hovertext=id,color="YlOrRd") %>% 
  add_markers() %>%
  layout(title = "MDS; DN_Trial2 proteins")
my_3dmds
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
my_3dpca<-plot_ly(m,x= ~PC1,y= ~PC2,z= ~PC3,hovertext=id) %>% 
  add_markers() %>%
  layout(title = "PCA; DN_Trial2 proteins")
my_3dpca
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
colnames(tsne_out)<-c("tsne_1","tsne2")
my_tsne<-ggplot(data=tsne_out$Y)+geom_point(aes(x=tsne_out$Y[,1],y=tsne_out$Y[,2],color="blue"))
my_tsne
tsne_out<-Rtsne(m,dim=3,check_duplicates=F)
colnames(tsne_out)<-c("tsne_1","tsne_2","tsne_3")
m<-data.frame(tsne_out$Y)
my_3dtsne<-plot_ly(m,x= ~tsne_1,y= ~tsne_2,z= ~tsne_3,
                   marker = list(size = 5),hovertext=id,color="YlOrRd") %>% 
  add_markers() %>%
  layout(title = "tSNE; DN_Trial2 proteins")
my_3dtsne

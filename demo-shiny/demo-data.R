mq.y.1 <- read.MQ.data(filename = system.file("extdata", "test_data", "proteinGroups_100mM_new.txt", package = "DepLab"), expt.id = "100mM", data.subset = "poi", organism = "yeast")
mq.y.3 <- read.MQ.data(filename = system.file("extdata", "test_data", "proteinGroups_300mM_new.txt", package = "DepLab"), expt.id = "300mM", data.subset = "poi", organism = "yeast")
mqcombi <- rbind(mq.y.1, mq.y.3)

smu <- superSmooth_values(long.df = subset(mqcombi, measurement == "raw.intensity"), prot.identifier = "gene_symbol")
fraction.norm <- normalize_values(long.df = smu,
                                  norm.type = "fraction", prot.identifier = "gene_symbol")


y.std.1 <- read.MQ.data(filename = system.file("extdata", "test_data", "proteinGroups_100mM_new.txt", package = "DepLab"), expt.id = "100mM", data.subset = "trypsin", organism = NULL)
y.std.3 <- read.MQ.data(filename = system.file("extdata", "test_data", "proteinGroups_300mM_new.txt", package = "DepLab"), expt.id = "300mM", data.subset = "trypsin", organism = NULL)
std.combi <- rbind(y.std.1, y.std.3)

std.norm <- normalize_values(long.df = subset(mqcombi, measurement == "raw.intensity"), norm.type = "spike-in", prot.identifier = "gene_symbol", std.df = subset(std.combi, measurement == "raw.intensity"))


mqcombi.plot <- subset(mqcombi, measurement == "raw.intensity" & gene_symbol %in% c("YAL003W (EFB1)", "YAL005C (SSA1)"))
P <- plot_profile(mqcombi.plot, what = c("gene_symbol","expt_id"), color.by = "gene_symbol", split.by = "gene_symbol", line.smooth = TRUE)
# print(P)

# plot normalized values
# mqcombi.plot.norm.frac <- subset(fraction.norm, gene_symbol %in% c("YAL016W (TPD3)", "YAL005C (SSA1)", "YAL003W (EFB1)"))
# plot_profile(mqcombi.plot.norm.frac, what = c("gene_symbol","expt_id"), color.by = "gene_symbol", split.by = "gene_symbol", line.smooth = FALSE)

mqcombi.plot.norm.std <- subset(std.norm, gene_symbol %in% c("YAL016W (TPD3)", "YAL005C (SSA1)", "YAL003W (EFB1)") )

my_profile_plot <- plot_profile(mqcombi.plot.norm.std, what = c("gene_symbol","expt_id"), color.by = "expt_id", split.by = "gene_symbol", line.smooth = FALSE)


ggplotly(my_profile_plot)

# adding marker for molecular weight
# mwmark <- data.frame(expt_id = c("100mM","300mM"), MWmarker = c(15,25))
# P + geom_vline(data = mwmark, aes(xintercept = MWmarker), linetype="dashed")

# smoothing line instead of geom_line
# ggplot(subset(fraction.norm, grepl("YAL005C", gene_symbol)), aes(x=fraction, y=value, colour=expt_id)) + 
#     geom_smooth(span = 0.3, se = TRUE) +
#     geom_point()+ theme_bw() +
#     facet_grid(.~expt_id)




# make a heatmap the old way
# data from the data package
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
library("NMF")
# pdf(file = "./heatmpa.pdf")
# aheatmap(log1p(wt1_mat[1:500,]), Colv = NA, distfun = "pearson",
#          scale = "none")
# dev.off()

# install.packages(c("heatmaply", "shinyHeatmaply"))
library("heatmaply")
library("shinyHeatmaply")
my_heatmap <- heatmaply(x = log1p(wt1_mat[1:500,]), distfun = function(x) dist((1-cor(t(x), method = "pearson"))), scale = "none", Colv = FALSE)

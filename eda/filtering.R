library(magrittr)
library(ggplot2)
library(DepLab)
library(tidyverse)
library(matrixStats)

DB <- "day9.db"
exps <- c(paste0("DN_9d_0", c(1:3)), paste0("EV_9d_0", c(1:3)), paste0("WT_9d_0",c(1:3)))
raw <- as.data.table(DepLab:::query.measurements.by.expt.with.gene.symbol.v2(DB, exps, "raw.intensity"))

## take the log
## raw$log_value <- log(raw$value + 1)
## raw[,log_value:= log(value)]

# filtered <- raw[raw[, .N, gene_symbol][N == 270]$gene_symbol, on = "gene_symbol"]

## normalizing and smoothing

# ## tidyversion
# filtered %>%
#   group_by(gene_symbol, expt_id) %>%
#   do(cbind(., smoothed=supsmu(x = seq(1, length(.$log_value)), y = .$log_value)$y)) -> filtered_smoothed
#
# filtered_smoothed %>%
#   filter(gene_symbol == "A0FGR8 (ESYT2)") -> tmp
#
# plot_profile(tmp, y = "log_value", color.by = "expt_id")
# plot_profile(tmp, y = "smoothed", color.by = "expt_id")

filtered2 <- raw[raw[, .N, gene_symbol][N == 270]$gene_symbol, on = "gene_symbol"]
filtered2$measurement <- "filtered"

filt.norm <- normalize_values(long.df = filtered2, norm.type = "fraction", prot.identifier = "gene_symbol")
filt.norm <- superSmooth_values(filt.norm, "gene_symbol") # this is a data.frame!

## old (log)
# filtered <- filtered[, value:= log(value + 1)]
# filtered$measurement <- "log_filtered"
# filtered <- superSmooth_values(filtered, "gene_symbol") # this is a data.frame!

# ## check the values
# filtered %>%
#   filter(gene_symbol == "A0FGR8 (ESYT2)") -> tmp2
#
# plot_profile(tmp2, y = "log_value", color.by = "expt_id")
# plot_profile(tmp2, y = "value", color.by = "expt_id")

## smoothing will create negative values --> putting them to 0
filt.norm %>%
  mutate(value = pmax(0, value)) -> filt.norm
head(filt.norm)

## wide data
wide_smu <- dcast(filt.norm, "gene_symbol + fraction ~ expt_id")
head(wide_smu)

## add condition column
filt.norm$condition <- gsub("_[0-9]+$","", filt.norm$expt_id)

pca <- prcomp(t(as.matrix(wide_smu[,3:11])))

condition_factor <- as.factor(rep(unique(filt.norm$condition), each=3))
plot(pca$x, pch=19, col=(2:4)[condition_factor])

## compute average pairwise correlation
filtered.dt <- data.table(filt.norm)
filtered.dt$condition <- gsub("_[0-9]+$","", filtered.dt$expt_id)
filtered.dt$replicate <- gsub(".*_([0-9]+$)", "\\1", filtered.dt$expt_id)

cors <- DepLab:::corrHM_wrap(filtered.dt[,-"measurement", with = FALSE], measurement="value", corr.method = "spearman", uniquifiers=c("gene_symbol", "condition", "replicate"), calc.mode = "both")
colMeans(cors$condition_comp)

mean_cor <- colMeans(cors$pairwise)
cor_mat <- matrix(data = NA, ncol=9, nrow=9)
diag(cor_mat) <- 1
cor_mat[2:9,1] <- mean_cor[1:8]
cor_mat[3:9, 2] <- mean_cor[9:15]
cor_mat[4:9, 3] <- mean_cor[16:21]
cor_mat[5:9, 4] <- mean_cor[22:26]
cor_mat[6:9, 5] <- mean_cor[27:30]
cor_mat[7:9, 6] <- mean_cor[31:33]
cor_mat[8:9, 7] <- mean_cor[34:35]
cor_mat[9, 8] <- mean_cor[36]
rownames(cor_mat) <- colnames(cor_mat) <- unique(filtered$expt_id)
pheatmap::pheatmap(cor_mat, cluster_rows = FALSE,
                   cluster_cols = FALSE, scale = "none")

## Filtering WT_9d_03
filt.norm %>%
  filter(expt_id != "WT_9d_03") -> filt.norm

## filtering based on correlations
correlations <- lapply(unique(filt.norm$condition),
                       function(cond) {
    filt.norm %>%
    filter(condition == cond) %>%
    dcast("gene_symbol + fraction ~ expt_id") %>%
    group_by(gene_symbol) %>%
    do(data.frame(min_cor =  min(cor(.[,-(1:2)])))) %>%
    mutate(condition = cond)
})

correlations <- do.call(rbind, correlations)

ggplot(correlations, aes(x = condition, y = min_cor)) +
  geom_boxplot()

correlations %>%
  group_by(condition) %>%
  summarize(q25 = quantile(min_cor, probs = .25, na.rm=TRUE))

correlations %>%
  group_by(gene_symbol) %>%
  summarize(cor = min(min_cor)) %>%
  filter(cor >= .9) %>%
  pull(gene_symbol) %>% as.character -> good_proteins

filt.norm %>%
  filter(gene_symbol %in% good_proteins) -> day9_filtered
save(day9_filtered, file="day9_filtered_norm.rda")


### peaks
find_peaks <- function(x, thresh = 0) {
  pks <- which(diff(sign(diff(x)))==-2) + 1
##  which_pks <- (x[pks] - x[max(0, pks - 3)] > thresh) &
##    (x[pks] - x[min(length(x), pks + 3)] > thresh)
  return(pks) #[which_pks])
}

filt.norm %>%
     group_by(gene_symbol, expt_id) %>%
     do(data.frame(peak = find_peaks(.$value))) -> peaks

num_peaks <- na.omit(dcast(peaks, formula = "gene_symbol ~ expt_id",
                           fun.aggregate = length,
                           value.var = "peak"))

npeaks <- as.matrix(num_peaks[,-1])
rownames(npeaks) <- num_peaks[,1]
npeaks <- npeaks[rowVars(npeaks)>0,]

library(pheatmap)
pheatmap(npeaks, scale = "none", distfun = "manhattan", labels_row = FALSE)

mean_npeaks <- t(apply(npeaks, 1, tapply, condition_factor[-9], mean))
colnames(mean_npeaks) <- c("dn", "ev", "wt")
head(mean_npeaks)

pheatmap(mean_npeaks, scale = "none", distfun = "manhattan")

filt.norm %>%
  group_by(gene_symbol, expt_id) %>%
  summarize(peak = which.max(value)) -> peaks

best_peaks <- na.omit(dcast(peaks, formula = "gene_symbol ~ expt_id",
                            value.var = "peak"))

bestpeaks <- as.matrix(best_peaks[,-1])
rownames(bestpeaks) <- best_peaks[,1]
bestpeaks <- bestpeaks[rowVars(bestpeaks)>0,]

pheatmap(bestpeaks, scale = "none", labels_row = FALSE)

head(mean_npeaks)

filt.norm %>%
  filter(gene_symbol == "A3KMH1 (VWA8)") %>%
  ggplot(aes(x = fraction, y = value, color = expt_id,
             group=expt_id)) +
  geom_line() + facet_grid(~condition)
head(npeaks)

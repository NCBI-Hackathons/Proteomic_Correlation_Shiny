library(magrittr)
library(ggplot2)
library(DepLab)
library(tidyverse)
library(matrixStats)

DB <- "day9.db"
exps <- c(paste0("DN_9d_0", c(1:3)), paste0("EV_9d_0", c(1:3)), paste0("WT_9d_0",c(1:3)))
raw <- as.data.table(DepLab:::query.measurements.by.expt.with.gene.symbol.v2(DB, exps, "raw.intensity"))

## take the log
raw$log_value <- log(raw$value + 1)
## raw[,log_value:= log(value)]

filtered <- raw[raw[, .N, gene_symbol][N == 270]$gene_symbol, on = "gene_symbol"]

## smoothing

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

filtered <- raw[raw[, .N, gene_symbol][N == 270]$gene_symbol, on = "gene_symbol"]
filtered <- filtered[, value:= log(value + 1)]
filtered$measurement <- "log_filtered"
filtered <- superSmooth_values(filtered, "gene_symbol") # this is a data.frame!

# ## check the values
# filtered %>%
#   filter(gene_symbol == "A0FGR8 (ESYT2)") -> tmp2
#
# plot_profile(tmp2, y = "log_value", color.by = "expt_id")
# plot_profile(tmp2, y = "value", color.by = "expt_id")

## smoothing will create negative values --> putting them to 0
filtered %>%
  mutate(value = pmax(0, value)) -> filtered
head(filtered)

## wide data
wide_smu <- na.omit(dcast(filtered, "gene_symbol + fraction ~ expt_id"))
head(wide_smu)

## add condition column
filtered$condition <- gsub("_[0-9]+$","", filtered$expt_id)

pca <- prcomp(t(as.matrix(wide_smu[,3:11])))

condition_factor <- as.factor(rep(unique(filtered$condition), each=3))
plot(pca$x, pch=19, col=(2:4)[condition_factor])

## compute average pairwise correlation
filtered.dt <- data.table(filtered)
filtered.dt$condition <- gsub("_[0-9]+$","", filtered.dt$expt_id)
filtered.dt$replicate <- gsub(".*_([0-9]+$)", "\\1", filtered.dt$expt_id)

cors <- DepLab:::corrHM_wrap(filtered.dt[,-"measurement", with = FALSE], measurement="value", corr.method = "pearson", uniquifiers=c("gene_symbol", "condition", "replicate"), calc.mode = "both")
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
filtered %>%
  filter(expt_id != "WT_9d_03") -> filtered

## filtering based on correlations
correlations <- lapply(unique(filtered$condition),
                       function(cond) {
    filtered %>%
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
  filter(cor >= .5) %>%
  pull(gene_symbol) %>% as.character -> good_proteins

filtered %>%
  filter(gene_symbol %in% good_proteins) -> day9_filtered
save(day9_filtered, file="day9_filtered.rda")

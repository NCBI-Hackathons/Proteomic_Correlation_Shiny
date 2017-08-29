library(magrittr)
library(ggplot2)
library(DepLab)
library(tidyverse)
library(matrixStats)

load("day9_filtered_norm.rda")

day9_filtered %>%
  group_by(gene_symbol, condition, fraction) %>%
  summarize(mean_value = mean(value)) -> mean_profiles

comparisons <- list("EVvDN" = c("EV_9d", "DN_9d"),
                    "EVvWT" = c("EV_9d", "WT_9d"),
                    "WTvDN" = c("WT_9d", "DN_9d"))

bw_group_cors <- lapply(comparisons, function(comp) {
  mean_profiles %>%
    filter(condition %in% comp) %>%
    dcast("gene_symbol + fraction ~ condition") %>%
    group_by(gene_symbol) %>%
    set_colnames(c("gene_symbol", "fraction", "cond1", "cond2")) %>%
    do(data.frame(cor = cor(.$cond1, .$cond2, method = "spearman"))) %>%
    arrange(cor)
})
bw_group_cors

mean_profiles %>%
  filter(gene_symbol == "Q08378 (GOLGA3)") %>%
  ggplot(aes(x = fraction, y = mean_value, color = condition)) +
  geom_line()


day9_filtered %>%
  filter(gene_symbol == "Q08378 (GOLGA3)") %>%
  ggplot(aes(x = fraction, y = value, color = expt_id,
             group=expt_id)) +
  geom_line() + facet_grid(~condition)

## Here we compute correlations to look at proteins that change across conditions. 

## The input is the filtered and smoothed data and the output is a ranking of proteins from the most interesting to the least interesting

library(magrittr)
library(ggplot2)
library(DepLabData)
library(DepLab)
library(tidyverse)
library(matrixStats)

data("WT_trial1")
data("WT_trial2")
data("DN_trial1")
data("DN_trial2")
data("EV_trial1")
data("EV_trial2")

WT_trial1 %>%
  MQ_to_longFormat(., y = "raw.intensity", return.dt = TRUE,
                   extract_proteinID(.$Protein.IDs, routine = "human")) -> wild_type_1
WT_trial2 %>%
  MQ_to_longFormat(., y = "raw.intensity", return.dt = TRUE,
                   extract_proteinID(.$Protein.IDs, routine = "human")) -> wild_type_2

DN_trial1 %>%
  MQ_to_longFormat(., y = "raw.intensity", return.dt = TRUE,
                   extract_proteinID(.$Protein.IDs, routine = "human")) -> dn_1
DN_trial2 %>%
  MQ_to_longFormat(., y = "raw.intensity", return.dt = TRUE,
                   extract_proteinID(.$Protein.IDs, routine = "human")) -> dn_2

EV_trial1 %>%
  MQ_to_longFormat(., y = "raw.intensity", return.dt = TRUE,
                   extract_proteinID(.$Protein.IDs, routine = "human")) -> ev_1
EV_trial2 %>%
  MQ_to_longFormat(., y = "raw.intensity", return.dt = TRUE,
                   extract_proteinID(.$Protein.IDs, routine = "human")) -> ev_2

## combine samples
wild_type_1 %>%
  mutate(expt_id = "wt1", exp_cond = "wt") -> wt1
wild_type_2 %>%
  mutate(expt_id = "wt2", exp_cond = "wt") -> wt2
dn_1 %>%
  mutate(expt_id = "dn1", exp_cond = "dn") -> dn1
dn_2 %>%
  mutate(expt_id = "dn2", exp_cond = "dn") -> dn2
ev_1 %>%
  mutate(expt_id = "ev1", exp_cond = "ev") -> ev1
ev_2 %>%
  mutate(expt_id = "ev2", exp_cond = "ev") -> ev2

long_df <- rbind(wt1, wt2, dn1, dn2, ev1, ev2)

## smoothing
smu <- superSmooth_values(long.df = long_df, prot.identifier = "id")

## smoothing will create negative values --> putting them to 0
smu %>%
  mutate(value = pmax(0, value)) -> smu

## wide data
wide_smu <- na.omit(dcast(smu, "id + fraction ~ expt_id"))

# Filtering

protein_list <- unique(as.character(wide_smu$id))
length(protein_list)

smu %>%
  filter(id %in% protein_list) -> smu_filtered
smu_filtered %>%
  mutate(value = log1p(value)) -> smu_log

conditions <- c("ev", "wt", "dn")
correlations <- list()
correlations <- lapply(conditions, function(cond) {
  smu_log %>%
    filter(exp_cond == cond) %>%
    dcast("id + fraction ~ expt_id") %>%
    group_by(id) %>%
    do(data.frame(min_cor =  min(cor(.[,-(1:2)])))) %>%
    mutate(exp_cond = cond)
})

correlations <- do.call(rbind, correlations)

correlations %>%
  group_by(id) %>%
  summarize(cor = min(min_cor)) %>%
  filter(cor >= .5) %>%
  pull(id) %>% as.character -> good_proteins

smu_filtered %>%
  filter(id %in% good_proteins) -> smu_good

# Correlation

smu_good %>%
  group_by(id, exp_cond, fraction) %>%
  summarize(mean_value = mean(value)) -> mean_profiles

mean_profiles %>%
  mutate(mean_value = log1p(mean_value)) -> mean_profiles_log

comparisons <- list("EVvDN" = c("ev", "dn"),
                    "EVvWT" = c("ev", "wt"),
                    "WTvDN" = c("wt", "dn"))

bw_group_cors <- lapply(comparisons, function(comp) {
  mean_profiles_log %>%
    filter(exp_cond %in% comp) %>%
    dcast("id + fraction ~ exp_cond") %>%
    group_by(id) %>%
    set_colnames(c("id", "fraction", "cond1", "cond2")) %>%
    do(data.frame(cor = cor(.$cond1, .$cond2))) %>%
    arrange(cor)
})

save(bw_group_cors, file="protein_list.rda")

find_peaks <- function(x, thresh = 1e3) {
  pks <- which(diff(sign(diff(x)))==-2) + 1
  which_pks <- (x[pks] - x[max(0, pks - 3)] > thresh) &
    (x[pks] - x[min(length(x), pks + 3)] > thresh)
  return(pks[which_pks])
}

smu_good %>%
  group_by(id, expt_id) %>%
  do(data.frame(peak = find_peaks(.$value))) -> peaks

num_peaks <- na.omit(dcast(peaks, formula = "id ~ expt_id",
                           fun.aggregate = length,
                           value.var = "peak"))

npeaks <- as.matrix(num_peaks[,-1])
rownames(npeaks) <- num_peaks[,1]
npeaks <- npeaks[rowVars(npeaks)>0,]

x <- gl(3, 2)
mean_npeaks <- t(apply(npeaks, 1, tapply, x, mean))
colnames(mean_npeaks) <- c("dn", "ev", "wt")
head(mean_npeaks)

smu_good %>%
  group_by(id, expt_id) %>%
  summarize(peak = which.max(value)) -> peaks

best_peaks <- na.omit(dcast(peaks, formula = "id ~ expt_id",
                            value.var = "peak"))

bestpeaks <- as.matrix(best_peaks[,-1])
rownames(bestpeaks) <- best_peaks[,1]
bestpeaks <- bestpeaks[rowVars(bestpeaks)>0,]

x <- gl(3, 2)
mean_bestpeaks <- t(apply(bestpeaks, 1, tapply, x, mean))
colnames(mean_bestpeaks) <- c("dn", "ev", "wt")
head(mean_bestpeaks)

## NMF::aheatmap(mean_npeaks, scale = "none", distfun = "manhattan")
## NMF::aheatmap(mean_bestpeaks, scale = "none", distfun = "manhattan")

save(mean_npeaks, mean_bestpeaks, file="peaks.rda")

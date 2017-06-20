library(magrittr)
library(ggplot2)
library(DepLabData)
library(DepLab)
library(tidyverse)

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

data_list = list(wt1 = rename(wild_type_1, wt1 = value),
                 wt2 = rename(wild_type_2, wt2 = value),
                 dn1 = rename(dn_1, dn1 = value),
                 dn2 = rename(dn_2, dn2 = value),
                 dn1 = rename(ev_1, ev1 = value),
                 dn2 = rename(ev_2, ev2 = value))

data_all = reduce(data_list, left_join) %>%
  select(id, fraction, measurement, everything()) %>%
  arrange(id, fraction)
filtered <- na.omit(data_all)

filtered %>%
  group_by(id) %>%
  summarize(cor_wt = cor(wt1, wt2),
            cor_dn = cor(dn1, dn2),
            cor_wt1dn1 = cor(wt1, dn1),
            cor_wt1dn2 = cor(wt1, dn2),
            cor_wt2dn2 = cor(wt2, dn2)) %>%
  na.omit -> cors

plot(density(cors$cor_dn), main="DN")
plot(density(cors$cor_wt), main="WT")
plot(density(cors$cor_wt1dn1), main="WT1 vs DN1")
plot(density(cors$cor_wt1dn2), main="WT1 vs DN2")

head(cors)
cors_all <- cor(filtered[,4:9])

pdf("global_cor.pdf")
NMF::aheatmap(cors_all)
dev.off()

filtered %>%
  group_by(id) %>%
  summarize(max_wt1 = which.max(wt1),
            max_wt2 = which.max(wt2),
            max_dn1 = which.max(dn1),
            max_dn2 = which.max(dn2),
            max_ev1 = which.max(ev1),
            max_ev2 = which.max(ev2)
            ) -> Modes
cors_mode <- cor(Modes[,2:7])
heatmap(cors_mode)

modmat <- as.matrix(Modes[,-1])
rownames(modmat) <- Modes$id

library(limma)
x <- gl(3, 2)
x
design <- model.matrix(~x)
fit <- lmFit(modmat, design)
fit <- eBayes(fit)
topTable(fit)

library(cowplot)

filtered %>%
  filter(id == "Q9Y6M7") %>%
  ggplot(aes(x = fraction)) +
  geom_line(aes(y = wt1, color = "wt")) +
  geom_line(aes(y = wt2, color = "wt")) +
  geom_line(aes(y = dn1, color = "dn")) +
  geom_line(aes(y = dn2, color = "dn")) +
  geom_line(aes(y = ev1, color = "ev")) +
  geom_line(aes(y = ev2, color = "ev")) -> Q9Y6M7

save_plot("Q9Y6M7_profiles.pdf", Q9Y6M7)

filtered %>%
  filter(id == "P60953") %>%
  ggplot(aes(x = fraction)) +
  geom_line(aes(y = wt1, color = "wt")) +
  geom_line(aes(y = wt2, color = "wt")) +
  geom_line(aes(y = dn1, color = "dn")) +
  geom_line(aes(y = dn2, color = "dn")) +
  geom_line(aes(y = ev1, color = "ev")) +
  geom_line(aes(y = ev2, color = "ev")) -> P60953

save_plot("P60953_profiles.pdf", P60953)

filtered %>%
  filter(id == "P85037") %>%
  ggplot(aes(x = fraction)) +
  geom_line(aes(y = wt1, color = "wt")) +
  geom_line(aes(y = wt2, color = "wt")) +
  geom_line(aes(y = dn1, color = "dn")) +
  geom_line(aes(y = dn2, color = "dn")) +
  geom_line(aes(y = ev1, color = "ev")) +
  geom_line(aes(y = ev2, color = "ev")) -> P85037

save_plot("P85037_profiles.pdf", P85037)


### smoothing before the analysis

## add experimental info

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

smu <- superSmooth_values(long.df = long_df, prot.identifier = "id")

smu %>%
  filter(id == "P85037") %>%
  ggplot(aes(x = fraction, y = value, group = expt_id,
             color = exp_cond)) +
  geom_line()

long_df %>%
  filter(id == "P85037") %>%
  ggplot(aes(x = fraction, y = value, group = expt_id,
             color = exp_cond)) +
  geom_line()

wide_df <- na.omit(dcast(long_df, "id + fraction ~ expt_id"))
head(wide_df)
cors_df <- cor(wide_df[,-(1:2)])
cors_df

wide_smu <- na.omit(dcast(smu, "id + fraction ~ expt_id"))
cors_smu <- cor(wide_smu[,-(1:2)])
cors_smu

heatmap(cors_df)
heatmap(cors_smu)

pdf("global_cor_smoothed.pdf")
NMF::aheatmap(cors_smu)
dev.off()

wide_smu %>%
  group_by(id) %>%
  summarize(max_wt1 = which.max(wt1),
            max_wt2 = which.max(wt2),
            max_dn1 = which.max(dn1),
            max_dn2 = which.max(dn2),
            max_ev1 = which.max(ev1),
            max_ev2 = which.max(ev2)
  ) -> Modes
cors_mode <- cor(Modes[,2:7])
heatmap(cors_mode)

modmat <- as.matrix(Modes[,-1])
rownames(modmat) <- Modes$id

wide_smu %>%
  group_by(id)

### plan:
# 1. filter out proteins not correlated between replicates
# 2. identify proteins that are not correlated between conditions (already in shiny app)
# 3. identify peaks (pdfCluster or mclust)
# 4. compare peak locations (what if there are multiple peaks?)

# the below function was inspired by the quantmod package (findPeak)
find_peaks <- function(x, thresh = 1e3) {
  pks <- which(diff(sign(diff(x)))==-2) + 1
  which_pks <- (x[pks] - x[pks - 3] > thresh) &
    (x[pks] - x[pks + 3] > thresh)
  return(pks[which_pks])
}

smu %>%
  filter(id == "Q7Z417") %>%
  ggplot(aes(x = fraction, y = value, group = expt_id,
             color = exp_cond)) +
  geom_line()


smu %>%
  filter(id == "P85037") %>%
  ggplot(aes(x = fraction, y = value, group = expt_id,
             color = exp_cond)) +
  geom_line()

  smu %>%
  filter(id %in% c("Q7Z417", "P85037")) %>%
  group_by(id, expt_id) %>%
  do(as.data.frame(find_peaks(.$value)))


library(limma)
x <- gl(3, 2)
x
design <- model.matrix(~x)
fit <- lmFit(modmat, design)
fit <- eBayes(fit)
round(topTable(fit)[,1:2])

smu %>%
  filter(id == "Q7Z417") %>%
  ggplot(aes(x = fraction, y = value, group = expt_id,
             color = exp_cond)) +
  geom_line()

smu %>%
  filter(id == "Q7Z417" & expt_id == "wt1") -> test_data

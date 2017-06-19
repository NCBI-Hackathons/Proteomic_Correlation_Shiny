library(magrittr)
library(ggplot2)
library(DepLabData)
library(DepLab)
library(tidyverse)

data("WT_trial1")
data("WT_trial2")
data("DN_trial1")
data("DN_trial2")

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

data_list = list(wt1 = rename(wild_type_1, wt1 = value),
                 wt2 = rename(wild_type_2, wt2 = value),
                 dn1 = rename(dn_1, dn1 = value),
                 dn2 = rename(dn_2, dn2 = value))

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
cors_all <- cor(filtered[,4:7])
heatmap(cors_all)

filtered %>%
  group_by(id) %>%
  summarize(max_wt1 = which.max(wt1),
            max_wt2 = which.max(wt2),
            max_dn1 = which.max(dn1),
            max_dn2 = which.max(dn2)
            ) -> Modes
cors_mode <- cor(Modes[,2:5])
heatmap(cors_mode)

modmat <- as.matrix(Modes[,-1])
rownames(modmat) <- Modes$id

library(limma)
x <- gl(2, 2)
x
design <- model.matrix(~x)
fit <- lmFit(modmat, design)
fit <- eBayes(fit)
topTable(fit)

filtered %>%
  filter(id == "Q8TB61") %>%
  ggplot(aes(x = fraction)) +
  geom_line(aes(y = wt1, color = "wt")) +
  geom_line(aes(y = wt2, color = "wt")) +
  geom_line(aes(y = dn1, color = "dn")) +
  geom_line(aes(y = dn2, color = "dn"))

filtered %>%
  filter(id == "P31943") %>%
  ggplot(aes(x = fraction)) +
  geom_line(aes(y = wt1, color = "wt")) +
  geom_line(aes(y = wt2, color = "wt")) +
  geom_line(aes(y = dn1, color = "dn")) +
  geom_line(aes(y = dn2, color = "dn"))

diff <- apply(modmat, 1, function(x) mean(x[1:2]) - mean(x[3:4]))
head(sort(abs(diff), decreasing = TRUE))

filtered %>%
  filter(id == "Q9UPN9") %>%
  ggplot(aes(x = fraction)) +
  geom_line(aes(y = wt1, color = "wt")) +
  geom_line(aes(y = wt2, color = "wt")) +
  geom_line(aes(y = dn1, color = "dn")) +
  geom_line(aes(y = dn2, color = "dn"))






#### Admixture plot, tidyverse version ####
library(tidyverse)
library(here)

#### Sid all data
#maya sid all
#CV error (K=1): 0.51004
#CV error (K=2): 0.43361
#CV error (K=3): 0.43744
#CV error (K=4): 0.45047
#CV error (K=5): 0.47139
#so k = 2,3,4 all good opts

#k = 2:
#         Pop0
# Pop0
# Pop1    0.205

#k = 3:
#         Pop0    Pop1
# Pop0
# Pop1    0.078
# Pop2    0.183   0.217

#k = 4
#         Pop0    Pop1    Pop2
# Pop0
# Pop1    0.148
# Pop2    0.213   0.185
# Pop3    0.207   0.227   0.077

#### Inputs ####
qopt_path_sid_all <- here("2bRAD/admix_sid_all/mydata_noclones_k4.qopt")   # ngsAdmix output; swap in the k2 file for K=2
meta_path_sid_all <- here("2bRAD/admix_sid_all/sid_all_metadata_noclones.csv")

# colors: one per cluster, so the vector length must equal K
# cols_lineage_k2 <- c("#bcbddc", "#807dba", "#3f007d")
# cols_lineage_k2 <- c("#3f007d", "#807dba")

#### Read admixture proportions and attach sample metadata ####
# The qopt rows are in the same order as the bam list, so the metadata must be too.
# Metadata columns 1 and 2 are taken as individual ID and population/group.
admix_sid_all <- read_table(qopt_path_sid_all, col_names = FALSE, show_col_types = FALSE) %>%
  select(where(~ any(!is.na(.x)))) %>%   # ngsAdmix lines end in a space, which readr reads as an extra all-NA column
  set_names(paste0("cluster", seq_len(ncol(.)))) %>%
  bind_cols(
    read_csv(meta_path_sid_all, show_col_types = FALSE) %>%
      select(ind = sample_id, pop = lineage) %>%
      mutate(ind = as.character(ind))
  )

# Optional: set population order and rename to short codes, e.g.
# admix <- admix %>%
#   mutate(pop = fct_relevel(pop, "A", "B", "C")) %>%
#   mutate(pop = fct_recode(pop, "AA" = "A", "BB" = "B", "CC" = "C"))

#### Long format: one row per individual x cluster ####
admix_long_sid_all <- admix_sid_all %>%
  pivot_longer(starts_with("cluster"), names_to = "cluster", values_to = "proportion") %>%
  group_by(ind) %>%
  mutate(dominant = cluster[which.max(proportion)],
         max_prop = max(proportion)) %>%
  ungroup() %>%
  # order individuals within each population by dominant cluster, then by how strongly assigned
  arrange(pop, dominant, desc(max_prop)) %>%
  mutate(ind = fct_inorder(ind))

#### Admixture plot ####
plot_sid_all <- ggplot(admix_long_sid_all, aes(x = ind, y = proportion, fill = cluster)) +
  geom_col(width = 1) +
  facet_grid(~pop, scales = "free_x", space = "free_x") +
  #scale_fill_manual(values = cols_lineage_k2) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Ancestry proportion", fill = "Cluster") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.spacing = unit(0.1, "lines"),
    strip.background = element_blank()
  )
plot_sid_all

ggsave(plot_sid_all, filename = here("2bRAD/admixture_sid_all_noclones_k4.pdf"), h = 4, w = 15)

#### Cluster affiliations (clusters with proportion > 0.25, as in the original) ####
cluster_admix_sid_all <- admix_long_sid_all %>%
  group_by(ind) %>%
  summarize(cluster = paste(str_remove(cluster[proportion > 0.25], "cluster"),
                            collapse = "."))

#cluster_admix



#### Sid combo data
#### Inputs ####
qopt_path_sid_combo <- here("2bRAD/admix_sid_combo/mydata_noclones_k2.qopt")   # ngsAdmix output; swap in the k2 file for K=2
meta_path_sid_combo <- here("2bRAD/admix_sid_combo/sid_combo_metadata_noclones.csv")

# colors: one per cluster, so the vector length must equal K
# cols_lineage_k2 <- c("#bcbddc", "#807dba", "#3f007d")
# cols_lineage_k2 <- c("#3f007d", "#807dba")

#### Read admixture proportions and attach sample metadata ####
# The qopt rows are in the same order as the bam list, so the metadata must be too.
# Metadata columns 1 and 2 are taken as individual ID and population/group.
admix_sid_combo <- read_table(qopt_path_sid_combo, col_names = FALSE, show_col_types = FALSE) %>%
  select(where(~ any(!is.na(.x)))) %>%   # ngsAdmix lines end in a space, which readr reads as an extra all-NA column
  set_names(paste0("cluster", seq_len(ncol(.)))) %>%
  bind_cols(
    read_csv(meta_path_sid_combo, show_col_types = FALSE) %>%
      select(ind = sample_id, pop = lineage) %>%
      mutate(ind = as.character(ind))
  )

# Optional: set population order and rename to short codes, e.g.
# admix <- admix %>%
#   mutate(pop = fct_relevel(pop, "A", "B", "C")) %>%
#   mutate(pop = fct_recode(pop, "AA" = "A", "BB" = "B", "CC" = "C"))

#### Long format: one row per individual x cluster ####
admix_long_sid_combo <- admix_sid_combo %>%
  pivot_longer(starts_with("cluster"), names_to = "cluster", values_to = "proportion") %>%
  group_by(ind) %>%
  mutate(dominant = cluster[which.max(proportion)],
         max_prop = max(proportion)) %>%
  ungroup() %>%
  # order individuals within each population by dominant cluster, then by how strongly assigned
  arrange(pop, dominant, desc(max_prop)) %>%
  mutate(ind = fct_inorder(ind))

#### Admixture plot ####
plot_sid_combo <- ggplot(admix_long_sid_combo, aes(x = ind, y = proportion, fill = cluster)) +
  geom_col(width = 1) +
  facet_grid(~pop, scales = "free_x", space = "free_x") +
  #scale_fill_manual(values = cols_lineage_k2) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Ancestry proportion", fill = "Cluster") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.spacing = unit(0.1, "lines"),
    strip.background = element_blank()
  )
plot_sid_combo

ggsave(plot_sid_combo, filename = here("2bRAD/admixture_sid_combo_noclones_k2.pdf"), h = 4, w = 10)

#### Cluster affiliations (clusters with proportion > 0.25, as in the original) ####
cluster_admix_sid_combo <- admix_long_sid_combo %>%
  group_by(ind) %>%
  summarize(cluster = paste(str_remove(cluster[proportion > 0.25], "cluster"),
                                        collapse = "."))

#cluster_admix


### Siderastrea siderea #####
#### Inputs ####
qopt_path_sid <- here("2bRAD/admix_sid/mydata_k2.qopt")   # ngsAdmix output; swap in the k2 file for K=2
meta_path_sid <- here("2bRAD/admix_sid/sid_metadata.csv")

# colors: one per cluster, so the vector length must equal K
# cols_lineage_k2 <- c("#bcbddc", "#807dba", "#3f007d")
# cols_lineage_k2 <- c("#3f007d", "#807dba")

#### Read admixture proportions and attach sample metadata ####
# The qopt rows are in the same order as the bam list, so the metadata must be too.
# Metadata columns 1 and 2 are taken as individual ID and population/group.
admix_sid <- read_table(qopt_path_sid, col_names = FALSE, show_col_types = FALSE) %>%
  select(where(~ any(!is.na(.x)))) %>%   # ngsAdmix lines end in a space, which readr reads as an extra all-NA column
  set_names(paste0("cluster", seq_len(ncol(.)))) %>%
  bind_cols(
    read_csv(meta_path_sid, show_col_types = FALSE) %>%
      select(ind = sample_id, pop = reef_bay) %>%
      mutate(ind = as.character(ind))
  )

# Optional: set population order and rename to short codes, e.g.
# admix <- admix %>%
#   mutate(pop = fct_relevel(pop, "A", "B", "C")) %>%
#   mutate(pop = fct_recode(pop, "AA" = "A", "BB" = "B", "CC" = "C"))

#### Long format: one row per individual x cluster ####
admix_long_sid <- admix_sid %>%
  pivot_longer(starts_with("cluster"), names_to = "cluster", values_to = "proportion") %>%
  group_by(ind) %>%
  mutate(dominant = cluster[which.max(proportion)],
         max_prop = max(proportion)) %>%
  ungroup() %>%
  # order individuals within each population by dominant cluster, then by how strongly assigned
  arrange(pop, dominant, desc(max_prop)) %>%
  mutate(ind = fct_inorder(ind))

#### Admixture plot ####
plot_sid <- ggplot(admix_long_sid, aes(x = ind, y = proportion, fill = cluster)) +
  geom_col(width = 1) +
  facet_grid(~pop, scales = "free_x", space = "free_x") +
  #scale_fill_manual(values = cols_lineage_k2) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Ancestry proportion", fill = "Cluster") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.spacing = unit(0.1, "lines"),
    strip.background = element_blank()
  )

plot_sid

# ggsave(here("figures", "admixture_k2.pdf"), width = 8, height = 3)

#### Cluster affiliations (clusters with proportion > 0.25, as in the original) ####
cluster_admix_sid <- admix_long_sid %>%
  group_by(ind) %>%
  summarize(cluster = paste(str_remove(cluster[proportion > 0.25], "cluster"),
                                        collapse = "."))

#cluster_admix

#### Branching porites species #####

#CV error (K=1): 0.60727
#CV error (K=2): 0.43157
#CV error (K=3): 0.34114
#CV error (K=4): 0.36052
#CV error (K=5): 0.46157

#k = 3 is best or maybe k = 4

#### Inputs ####
qopt_path_por <- here("2bRAD/admix_por/mydata_k2.qopt")   # ngsAdmix output; swap in the k2 file for K=2
meta_path_por <- here("2bRAD/admix_por/por_metadata.csv")

# colors: one per cluster, so the vector length must equal K
# cols_lineage_k2 <- c("#bcbddc", "#807dba", "#3f007d")
# cols_lineage_k2 <- c("#3f007d", "#807dba")

#### Read admixture proportions and attach sample metadata ####
# The qopt rows are in the same order as the bam list, so the metadata must be too.
# Metadata columns 1 and 2 are taken as individual ID and population/group.
admix_por <- read_table(qopt_path_por, col_names = FALSE, show_col_types = FALSE) %>%
  select(where(~ any(!is.na(.x)))) %>%   # ngsAdmix lines end in a space, which readr reads as an extra all-NA column
  set_names(paste0("cluster", seq_len(ncol(.)))) %>%
  bind_cols(
    read_csv(meta_path_por, show_col_types = FALSE) %>%
      select(ind = sample_id, pop = reef_bay) %>%
      mutate(ind = as.character(ind))
  )

# Optional: set population order and rename to short codes, e.g.
# admix <- admix %>%
#   mutate(pop = fct_relevel(pop, "A", "B", "C")) %>%
#   mutate(pop = fct_recode(pop, "AA" = "A", "BB" = "B", "CC" = "C"))

#### Long format: one row per individual x cluster ####
admix_long_por <- admix_por %>%
  pivot_longer(starts_with("cluster"), names_to = "cluster", values_to = "proportion") %>%
  group_by(ind) %>%
  mutate(dominant = cluster[which.max(proportion)],
         max_prop = max(proportion)) %>%
  ungroup() %>%
  # order individuals within each population by dominant cluster, then by how strongly assigned
  arrange(pop, dominant, desc(max_prop)) %>%
  mutate(ind = fct_inorder(ind))

#### Admixture plot ####
plot_por <- ggplot(admix_long_por, aes(x = ind, y = proportion, fill = cluster)) +
  geom_col(width = 1) +
  facet_grid(~pop, scales = "free_x", space = "free_x") +
  #scale_fill_manual(values = cols_lineage_k2) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Ancestry proportion", fill = "Cluster") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.spacing = unit(0.1, "lines"),
    strip.background = element_blank()
  )

plot_por

ggsave(plot_por, filename = here("2bRAD/admixture_por_k3.pdf"), h = 4, w = 8)

#### Cluster affiliations (clusters with proportion > 0.25, as in the original) ####
cluster_admix_por <- admix_long_por %>%
  group_by(ind) %>%
  summarize(cluster = paste(str_remove(cluster[proportion > 0.25], "cluster"),
                            collapse = ".")) |> 
  rename(sample_id = ind) |> 
  mutate(sample_id = as.character(sample_id))
  

# por_meta <- read_csv(here("2bRAD/por/por_metadata.csv")) %>%
#   mutate(sample_id = as.character(sample_id))
# 
# lineage_por <- cluster_admix_por |> left_join(por_meta) 

write.csv(cluster_admix_por, here("2bRAD/por/lineage_por_k2.csv"), row.names=FALSE)
#edit to manually add lineage with "L1, L2, etc)


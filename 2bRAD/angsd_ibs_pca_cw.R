#Maya Powell 
###Process rewritten from Hannah Aichelman TVE IBS script###

#### Setup ####
library(tidyverse)
library(vegan)
library(here)
library(ggdendro)

### Branching Porites species

#read in meta
por_meta <- read_csv(here("2bRAD/por/por_metadata.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))
# id_col used below: adjust "sample_id" to whatever column matches your bam/ngsrelate order
lin_k2 <- read_csv(here("2bRAD/por/lineage_por_k2.csv"), show_col_types = FALSE) %>% 
  mutate(sample_id = as.character(sample_id), lineage = as.character(lineage)) %>%
  rename(lin_k2 = lineage) %>% select(sample_id, lin_k2)
lin_k3 <-read_csv(here("2bRAD/por/lineage_por_k3.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id), lineage = as.character(lineage)) %>%
  rename(lin_k3 = lineage) %>% select(sample_id, lin_k3)
lin_k4 <-read_csv(here("2bRAD/por/lineage_por_k4.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id), lineage = as.character(lineage)) %>%
  rename(lin_k4 = lineage) %>% select(sample_id, lin_k4)

por_lin <- por_meta |> left_join(lin_k2) |> 
  left_join(lin_k3) |> left_join(lin_k4)

#read in ibs
por_ibs <- as.matrix(read.table(here("2bRAD/por/myresult.ibsMat")))
dimnames(por_ibs) <- list(por_meta$sample_id, por_meta$sample_id)

#dendrogram
por_den <- ggdendrogram(hclust(as.dist(por_ibs), method = "average")) +
  labs(x = "Sample ID", y = "IBS distance") +
  #geom_hline(yintercept = 0.265) +
  theme_bw()
por_den

#pcoa + cap
por_pcoa <- capscale(por_ibs ~ 1)
por_cap  <- capscale(por_ibs ~ reef_bay, data = por_meta)
adonis2(por_ibs ~ reef_bay, data = por_meta)
#          Df SumOfSqs   R2      F Pr(>F)    
# Model     1  0.50217 0.44 17.286  0.001 ***
# Residual 22  0.63913 0.56                  
# Total    23  1.14130 1.00      

por_eig <- tibble(
  axis = seq_along(por_pcoa$CA$eig),
  eigenvalue = por_pcoa$CA$eig) %>%
  mutate(pct_var = eigenvalue / sum(eigenvalue[eigenvalue > 0]) * 100)

por_scores <- scores(por_pcoa, display = "sites") %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(por_lin, by = "sample_id")

#pcoa plot
por_pcoa_plot <- ggplot(por_scores, aes(x = MDS1, y = MDS2, color = reef_bay, shape = lin_k4)) +
  geom_point(size = 2, stroke = 1) +
  stat_ellipse(type = "t", linewidth = 1) +
  labs(x = paste0("MDS1 (", round(por_eig$pct_var[1], 1), "% variance)"),
       y = paste0("MDS2 (", round(por_eig$pct_var[2], 1), "% variance)"),
    color = "Site", title = "Branching Porites sp.") +
  theme_bw()
por_pcoa_plot

#### Por: relatedness ####
por_rel <- read.table(here("2bRAD/por/ngsrelate_vcf.res"), header = TRUE)
por_lineages <- read_csv(here("2bRAD/por/lineage_por_k2.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))

# por_meta must be in the same row order as the bam list fed into ngsrelate,
# since ngsrelate's a/b columns are 0-indexed positions in that list.
por_meta <- read_csv(here("2bRAD/por/por_metadata.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))

por_lookup <- por_meta %>%
  mutate(merge_id = row_number() - 1) %>%
  select(merge_id, sample_id) %>%
  mutate(sample_id = as.character(sample_id)) %>%
  left_join(por_lineages, by = c("sample_id"))

por_rel_annotated <- por_rel %>%
  left_join(por_lookup, by = c("a" = "merge_id")) %>%
  rename(sample_id_a = sample_id, lineage_a = lineage) %>%
  left_join(por_lookup, by = c("b" = "merge_id")) %>%
  rename(sample_id_b = sample_id, lineage_b = lineage) %>%
  mutate(lineage_comparison = paste(lineage_a, lineage_b, sep = "_")) %>%
  filter(!is.na(lineage_a), !is.na(lineage_b))

# heatmap
ggplot(por_rel_annotated, aes(x = sample_id_a, y = sample_id_b, fill = rab)) +
  geom_tile() +
  scale_fill_viridis_c(name = "rab") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# boxplot by lineage comparison
por_rab_boxplot <- por_rel_annotated %>%
  filter(lineage_comparison %in% c("L1_L1", "L2_L2", "L3_L3", "L4_L4")) %>%
  ggplot(aes(x = lineage_comparison, y = rab, fill = lineage_comparison)) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter() +
  labs(x = NULL, y = "Pairwise relatedness (rab)") +
  theme_bw() +
  theme(legend.position = "none")

por_rab_boxplot

# summary table
por_related_comp <- por_rel_annotated %>%
  group_by(lineage_comparison) %>%
  summarize(mean_rab = mean(rab, na.rm = TRUE), sd_rab = sd(rab, na.rm = TRUE), n = n())
por_related_comp

#### Siderastrea siderea ####

#read in meta
sid_meta <- read_csv(here("2bRAD/sid/sid_metadata.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))
# id_col used below: adjust "sample_id" to whatever column matches your bam/ngsrelate order

#read in ibs
sid_ibs <- as.matrix(read.table(here("2bRAD/sid/myresult.ibsMat")))
dimnames(sid_ibs) <- list(sid_meta$sample_id, sid_meta$sample_id)

#dendrogram
sid_den <- ggdendrogram(hclust(as.dist(sid_ibs), method = "average")) +
  labs(x = "Sample ID", y = "IBS distance") +
  #geom_hline(yintercept = 0.265) +
  theme_bw()
sid_den

#pcoa + cap
sid_pcoa <- capscale(sid_ibs ~ 1)
sid_cap  <- capscale(sid_ibs ~ reef_bay, data = sid_meta)
adonis2(sid_ibs ~ reef_bay, data = sid_meta)
#          Df SumOfSqs      R2      F Pr(>F)
# Model     1  0.03342 0.04963 1.1489  0.165
# Residual 22  0.64000 0.95037              
# Total    23  0.67342 1.00000  

sid_eig <- tibble(
  axis = seq_along(sid_pcoa$CA$eig),
  eigenvalue = sid_pcoa$CA$eig) %>%
  mutate(pct_var = eigenvalue / sum(eigenvalue[eigenvalue > 0]) * 100)

sid_scores <- scores(sid_pcoa, display = "sites") %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(sid_meta, by = "sample_id")

#pcoa plot
sid_pcoa_plot <- ggplot(sid_scores, aes(x = MDS1, y = MDS2, color = reef_bay)) +
  geom_point(size = 2, stroke = 1) +
  stat_ellipse(type = "t", linewidth = 1) +
  labs(x = paste0("MDS1 (", round(sid_eig$pct_var[1], 1), "% variance)"),
       y = paste0("MDS2 (", round(sid_eig$pct_var[2], 1), "% variance)"),
       color = "Site", title = "Siderastrea siderea") +
  theme_bw()
sid_pcoa_plot

# #### sid: relatedness ####
# sid_rel <- read.table(here("2bRAD/sid/ngsrelate.res"), header = TRUE)
# sid_lineages <- read_csv(here("2bRAD/sid/sid_lineages.csv"), show_col_types = FALSE)
# 
# # sid_meta must be in the same row order as the bam list fed into ngsrelate,
# # since ngsrelate's a/b columns are 0-indexed positions in that list.
# sid_lookup <- sid_meta %>%
#   mutate(merge_id = row_number() - 1) %>%
#   select(merge_id, sample_id) %>%
#   left_join(sid_lineages, by = c("sample_id" = "gen_site"))
# 
# sid_rel_annotated <- sid_rel %>%
#   left_join(sid_lookup, by = c("a" = "merge_id")) %>%
#   rename(sample_id_a = sample_id, lineage_a = lineage) %>%
#   left_join(sid_lookup, by = c("b" = "merge_id")) %>%
#   rename(sample_id_b = sample_id, lineage_b = lineage) %>%
#   mutate(
#     sample_id_a = str_remove(sample_id_a, "_CLONE$"),
#     sample_id_b = str_remove(sample_id_b, "_CLONE$"),
#     lineage_comparison = paste(lineage_a, lineage_b, sep = "_")
#   ) %>%
#   filter(!is.na(lineage_a), !is.na(lineage_b))
# 
# # heatmap
# ggplot(sid_rel_annotated, aes(x = sample_id_a, y = sample_id_b, fill = rab)) +
#   geom_tile() +
#   scale_fill_viridis_c(name = "rab") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # boxplot by lineage comparison
# sid_rel_annotated %>%
#   filter(lineage_comparison %in% c("L1_L1", "L2_L2", "L3_L3")) %>%
#   ggplot(aes(x = lineage_comparison, y = rab, fill = lineage_comparison)) +
#   geom_boxplot() +
#   labs(x = NULL, y = "Pairwise relatedness (rab)") +
#   theme_bw() +
#   theme(legend.position = "none")
# 
# # summary table
# sid_rel_annotated %>%
#   group_by(lineage_comparison) %>%
#   summarize(mean_rab = mean(rab, na.rm = TRUE), sd_rab = sd(rab, na.rm = TRUE), n = n())

#### Siderastrea siderea - combo with Belize Data and Hannah Data ####

#combine metadata to get Hannah data
#phys metadata with lineage info
#sra metadata with sample info
# phys_meta <- read_csv(here("2bRAD/sid_combo/phys_metadata.csv"), show_col_types = FALSE) |> 
#   rename(sample_id = gen_site) |> 
#   group_by(sample_id) |> 
#   distinct(sample_id, .keep_all = TRUE)
# sra_meta <- read_csv(here("2bRAD/sid_combo/SraRunTable.csv"), show_col_types = FALSE) |> 
#   select(bam,`Sample Name`) |> 
#   rename(sample_id = `Sample Name` ) 
# all_meta <- sra_meta |> left_join(phys_meta, by = "sample_id") |> 
#   select(bam, sample_id, sitename,reef,lineage)
# write.csv(all_meta, here("2bRAD/sid_combo/hannah_metadata.csv"))
#NEED TO DO IT THIS WAY TO KEEP ALL SAMPLES IN SAME ORDER BC OTHERWISE ISSUES ARISE WITH MATCHING METADATA!!!
#order is order or samples

#read in meta
sid_combo_meta <- read_csv(here("2bRAD/sid_combo/sid_combo_metadata.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))
# id_col used below: adjust "sample_id" to whatever column matches your bam/ngsrelate order

# sid_combo_noclones <- sid_combo_meta |> filter(!remove %in% "y")
# sid_clones <- sid_combo_meta |> filter(remove == "y")
# 
# write.csv(sid_combo_noclones, here("2bRAD/sid_combo/sid_combo_metadata_noclones.csv"))
# write.csv(sid_clones, here("2bRAD/sid_combo/sid_clones.csv"))

#read in ibs
sid_combo_ibs <- as.matrix(read.table(here("2bRAD/sid_combo/myresult.ibsMat")))
dimnames(sid_combo_ibs) <- list(sid_combo_meta$sample_id, sid_combo_meta$sample_id)

#dendrogram
sid_combo_den <- ggdendrogram(hclust(as.dist(sid_combo_ibs), method = "average")) +
  labs(x = "Sample ID", y = "IBS distance") +
  #geom_hline(yintercept = 0.265) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
sid_combo_den

#pcoa + cap
sid_combo_pcoa <- capscale(sid_combo_ibs ~ 1)
sid_combo_cap  <- capscale(sid_combo_ibs ~ reef_bay, data = sid_combo_meta)
adonis2(sid_combo_ibs ~ reef_bay, data = sid_combo_meta)
#          Df SumOfSqs      R2      F Pr(>F)
# Model     1  0.03342 0.04963 1.1489  0.165
# Resid_comboual 22  0.64000 0.95037              
# Total    23  0.67342 1.00000  

sid_combo_eig <- tibble(
  axis = seq_along(sid_combo_pcoa$CA$eig),
  eigenvalue = sid_combo_pcoa$CA$eig) %>%
  mutate(pct_var = eigenvalue / sum(eigenvalue[eigenvalue > 0]) * 100)

sid_combo_scores <- scores(sid_combo_pcoa, display = "sites") %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(sid_combo_meta, by = "sample_id")

#pcoa plot
sid_combo_pcoa_plot <- ggplot(sid_combo_scores, aes(x = MDS1, y = MDS2, color = reef_bay, shape = lineage)) +
  geom_point(size = 2, stroke = 1) +
  stat_ellipse(type = "t", linewidth = 1) +
  labs(x = paste0("MDS1 (", round(sid_combo_eig$pct_var[1], 1), "% variance)"),
       y = paste0("MDS2 (", round(sid_combo_eig$pct_var[2], 1), "% variance)"),
       title = "Siderastrea siderea combo") +
  theme_bw()
sid_combo_pcoa_plot


#### Siderastrea siderea - combo with Belize Data and Hannah Data ####

#combine metadata to get Hannah data
#phys metadata with lineage info
#sra metadata with sample info
# phys_meta <- read_csv(here("2bRAD/sid_combo/phys_metadata.csv"), show_col_types = FALSE) |> 
#   rename(sample_id = gen_site) |> 
#   group_by(sample_id) |> 
#   distinct(sample_id, .keep_all = TRUE)
# sra_meta <- read_csv(here("2bRAD/sid_combo/SraRunTable.csv"), show_col_types = FALSE) |> 
#   select(bam,`Sample Name`) |> 
#   rename(sample_id = `Sample Name` ) 
# all_meta <- sra_meta |> left_join(phys_meta, by = "sample_id") |> 
#   select(bam, sample_id, sitename,reef,lineage)
# write.csv(all_meta, here("2bRAD/sid_combo/hannah_metadata.csv"))
#NEED TO DO IT THIS WAY TO KEEP ALL SAMPLES IN SAME ORDER BC OTHERWISE ISSUES ARISE WITH MATCHING METADATA!!!
#order is order or samples

#read in meta
sid_all_meta <- read_csv(here("2bRAD/sid_all/sid_all_metadata.csv"), show_col_types = FALSE) %>%
  mutate(sample_id = as.character(sample_id))
# id_col used below: adjust "sample_id" to whatever column matches your bam/ngsrelate order

# sid_all_noclones <- sid_all_meta |> filter(!remove %in% "y")
# sid_clones <- sid_all_meta |> filter(remove == "y")
# 
# write.csv(sid_all_noclones, here("2bRAD/sid_all/sid_all_metadata_noclones.csv"))
# write.csv(sid_clones, here("2bRAD/sid_all/sid_clones.csv"))

#read in ibs
sid_all_ibs <- as.matrix(read.table(here("2bRAD/sid_all/myresult.ibsMat")))
dimnames(sid_all_ibs) <- list(sid_all_meta$sample_id, sid_all_meta$sample_id)

#dendrogram
sid_all_den <- ggdendrogram(hclust(as.dist(sid_all_ibs), method = "average")) +
  labs(x = "Sample ID", y = "IBS distance") +
  #geom_hline(yintercept = 0.265) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
sid_all_den

#pcoa + cap
sid_all_pcoa <- capscale(sid_all_ibs ~ 1)
sid_all_cap  <- capscale(sid_all_ibs ~ reef_bay, data = sid_all_meta)
adonis2(sid_all_ibs ~ reef_bay, data = sid_all_meta)
#          Df SumOfSqs      R2      F Pr(>F)
# Model     1  0.03342 0.04963 1.1489  0.165
# Resid_allual 22  0.64000 0.95037              
# Total    23  0.67342 1.00000  

sid_all_eig <- tibble(
  axis = seq_along(sid_all_pcoa$CA$eig),
  eigenvalue = sid_all_pcoa$CA$eig) %>%
  mutate(pct_var = eigenvalue / sum(eigenvalue[eigenvalue > 0]) * 100)

sid_all_scores <- scores(sid_all_pcoa, display = "sites") %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(sid_all_meta, by = "sample_id")

#pcoa plot
sid_all_pcoa_plot <- ggplot(sid_all_scores, aes(x = MDS1, y = MDS2, color = reef_bay, shape = lineage)) +
  geom_point(size = 2, stroke = 1) +
  stat_ellipse(type = "t", linewidth = 1) +
  labs(x = paste0("MDS1 (", round(sid_all_eig$pct_var[1], 1), "% variance)"),
       y = paste0("MDS2 (", round(sid_all_eig$pct_var[2], 1), "% variance)"),
       title = "Siderastrea siderea all") +
  theme_bw()
sid_all_pcoa_plot

### Saving figures
ggsave(sid_den, filename = here("2bRAD/sid_den_plot.pdf"), h = 4, w = 6)
ggsave(sid_all_den, filename = here("2bRAD/sid_all_den_plot.pdf"), h = 4, w = 10)
ggsave(sid_combo_den, filename = here("2bRAD/sid_combo_den_plot.pdf"), h = 4, w = 10)
ggsave(por_den, filename = here("2bRAD/por_den_plot.pdf"), h = 4, w = 6)

ggsave(sid_pcoa_plot, filename = here("2bRAD/sid_pcoa_plot.pdf"), h = 4, w = 6)
ggsave(sid_all_pcoa_plot, filename = here("2bRAD/sid_all_pcoa_plot.pdf"), h = 4, w = 6)
ggsave(sid_combo_pcoa_plot, filename = here("2bRAD/sid_combo_pcoa_plot.pdf"), h = 4, w = 6)
ggsave(por_pcoa_plot, filename = here("2bRAD/por_pcoa_plot_lin_k4.pdf"), h = 4, w = 6)

ggsave(por_rab_boxplot, filename = here("2bRAD/por_rab_k2_plot.pdf"), h = 4, w = 6)
ggsave(sid_rab_boxplot, filename = here("2bRAD/sid_rab_k4_plot.pdf"), h = 4, w = 6)

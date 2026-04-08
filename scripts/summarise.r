# Bradley April 2026
# Compare the replication of QTLs and novelty of effector genes nominated by Tenk10k vs IBDverse

# Libraries
library(tidyverse)
library(ggplot2)
library(patchwork)
library(data.table)
library(ggrepel)

# Hard code paths
repo.dir <- 'IBDVerse-sc-eQTL-code/'
data.dir <- paste0(repo.dir,'data/')
coloc.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/results/2025_06_11_IBDverse_coloc_all_gwas/collapsed/'
henry_mr <- 'data/tenk10k_crohns_mr_genes.tsv'
metaf = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/IBDVerse-sc-eQTL-code/data/eqtl_metadata.txt.gz" # Didn't want to back up this big file
source(paste0(repo.dir,'qtl_plot/helper_functions.R'))
pi1f = "pi1_pairs-Tenk10K/results/pi1_all.txt"
tenk_annotation = "data/tenk10k_cell_type_annotation.tsv"
out.dir = "results"
if(!dir.exists(out.dir)){
  dir.create(out.dir)
}
ibdverse_immune = c("B", "Plasma", "Myeloid", "T")
annotation_file = paste0(repo.dir, "/data/all_IBDverse_annotation_mastersheet.csv")
annotation_map = read.csv(annotation_file) %>% 
              rename(label_new = "JAMBOREE_ANNOTATION") %>% 
              mutate(
                  label_new = gsub("_", " ", label_new)
              ) %>% 
              select(leiden, Category, label_new)

#################################
# Assess novelty of Tenk10K colocs vs IBDverse, grouped by major population
#################################
# Load IBDverse colocs and subset (Cuomo )
known.coloc.df <- get_colocs(coloc.dir, known_ibd_only = TRUE) %>%
    filter(PP.H4.abf > 0.75)
ibdverse_sub <- known.coloc.df %>%
  filter(gwas_trait == "CD")

# Load the sets of genes nominated by Tenk10K (coloc and MR)
mr = read.delim(henry_mr)
mr_genes = unique(mr$ensembl_gene_id)

# Calculate the novelty of genes in TenK10K vs IBDverse
ibdverse_genes = ibdverse_sub %>%
  pull(phenotype_id) %>%
  unique()

all_sub_de = length(unique(mr$ensembl_gene_id))
all_sub_novel_de = length(setdiff(unique(mr$ensembl_gene_id), ibdverse_genes))

cats = unique(known.coloc.df$category_new)
cats = cats[-grep("All Cells", cats)]
tenk_not_novel_res = do.call(rbind, lapply(cats, function(x){
  print(paste0("Processing category: ", x))
  
  # Get n effector genes
  de = ibdverse_sub %>% 
    filter(category_new == !!x) %>% 
    pull(phenotype_id) %>% 
    unique()
  nde = length(de)

  # Get novel effector genes 
  notnovelde = intersect(de, mr_genes)

  # Make contingency
  a <- length(notnovelde)
  b <- length(setdiff(all_sub_novel_de, notnovelde))
  c <- length(setdiff(de, notnovelde))
  d <- length(setdiff(all_sub_de, union(de, all_sub_novel_de)))
  
  contingency <- matrix(c(a, b, c, d), nrow = 2,
                        dimnames = list(
                          disease_effector = c("in_dataset", "not_in_dataset"),
                          novel_disease_effector   = c("in_dataset", "not_in_dataset")
                        ))

  fishres = fisher.test(contingency)
  or = fishres$estimate
  res = data.frame(category = x, nde = nde, n_notnovelde = length(notnovelde), OR=or)
  return(res)

}))

# Append the proportion of cells from from each major population which come from blood samples
cellmeta = fread(metaf)
prop_blood = cellmeta %>%
  group_by(predicted_category) %>%
  summarise(
    total_cells = n(),
    blood_cells = sum(tissue == "blood")
  ) %>%
  mutate(
    prop_blood = blood_cells / total_cells
  ) %>%
  select(predicted_category, prop_blood) %>%
  mutate(predicted_category = gsub("_ct", "", predicted_category))

tenk_not_novel_res = tenk_not_novel_res %>%
  left_join(prop_blood, by = c("category" = "predicted_category")) 

write.table(tenk_not_novel_res, paste0(out.dir, "/tenk_not_novel_res.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

tenk_not_novel_res = tenk_not_novel_res %>% 
    filter(is.finite(OR)) # MESENCHYMAL has OR of Inf

lm_fit <- lm(OR ~ log10(prop_blood), data = tenk_not_novel_res)
slope <- signif(coef(lm_fit)[[2]], 3)
pval <- signif(summary(lm_fit)$coefficients["log10(prop_blood)", "Pr(>|t|)"], 3)
# Print
print(paste0("Slope: ", slope, ". p=", pval))

# Plot this
tenk_not_novel_res$annotation_type = "Major population"
ggplot(tenk_not_novel_res, aes(x = log10(prop_blood), y = OR)) + 
  geom_point(aes(fill = category, size = annotation_type), shape = 21, stroke = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
  geom_text_repel(aes(label = category), size = 4, max.overlaps = Inf, point.padding = unit(5, "lines")) +
  scale_size_manual(values = annot.class.sizes, name = 'Annotation granularity') +
  scale_fill_manual(values = umap.category.palette, name = 'Major population') +
  labs(
    x = "log10 fraction of cells derived from blood samples",
    y ="Odds ratio for IBDverse disease effector gene\nalready found in TenK10K"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(paste0(out.dir,"/Blood_proportion_vs_notnovelDE_in_Henry-OR.png"), width = 6.5, height = 5)


###############
# Quantifying replication rates of eQTLs identified in Tenk10K vs IBDverse
###############
pi1 = read.delim(pi1f) %>%
    filter(!is.na(pi1)) %>%
    left_join(
        read.delim(tenk_annotation), by = c("discovery" = "cell_type")) %>%
    mutate(
        label_machine = gsub("dMean__", "", dataset_id),
        label_machine = gsub("_all", "", label_machine)
    ) %>% 
    left_join(
        annot.mapping %>% 
            select(label_machine, label_new, category_new)
    ) 

# Build labelled subsets for each comparison, then combine
pi1_max_all = pi1 %>%
    group_by(discovery) %>%
    slice_max(pi1) %>%
    ungroup() %>%
    mutate(comparison = "IBDverse - All")

pi1_max_gut = pi1 %>%
    filter(
        !grepl("_blood", dataset_id),
        !grepl("_ct", dataset_id),
        !grepl("unannotated", dataset_id)
    ) %>%
    group_by(discovery) %>%
    slice_max(pi1) %>%
    ungroup() %>%
    mutate(comparison = "IBDverse - TI/Rectum only")

pi1_max_gut_nonimmune = pi1 %>%
    filter(
        !grepl("_blood", dataset_id),
        !grepl("_ct", dataset_id),
        !(category_new %in% ibdverse_immune),
        !grepl("unannotated", dataset_id)
    ) %>%
    group_by(discovery) %>%
    slice_max(pi1) %>%
    ungroup() %>%
    mutate(comparison = "IBDverse - TI/Rectum, non-immune only")

pi1_combined = bind_rows(pi1_max_all, pi1_max_gut, pi1_max_gut_nonimmune) %>%
    mutate(comparison = factor(comparison, levels = c("IBDverse - All", "IBDverse - TI/Rectum only", "IBDverse - TI/Rectum, non-immune only")))

median_combined = pi1_combined %>%
    group_by(major_cell_type, comparison) %>%
    summarise(median_pi1 = median(pi1), .groups = "drop")

ggplot(pi1_combined, aes(x = pi1, fill = major_cell_type)) +
    geom_histogram(bins = 50, color = "black") +
    facet_grid(major_cell_type ~ comparison, scales = "free_y") +
    geom_vline(data = median_combined,
               aes(xintercept = median_pi1),
               lty = "dashed", color = "black", size = 1) +
    geom_text(
        data = median_combined,
        aes(x = 0, y = Inf, label = paste0("median: ", signif(median_pi1, 2))),
        inherit.aes = FALSE, hjust = 0, vjust = 1.5, size = 3.5
    ) +
    theme_classic() +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold"),
        strip.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)
    ) +
    labs(x = "Maximum π1 per cell-type", y = "Frequency") +
    xlim(c(0, 1))

ggsave(paste0(out.dir, "/pi1_TenK10K_vs_IBDverse.png"), width = 12, height = 12)

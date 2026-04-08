# Bradley April 2026
# Compare the replication of QTLs and novelty of effector genes nominated by Tenk10k vs IBDverse

# Libraries
library(tidyverse)
library(ggplot2)
library(patchwork)

# Hard code paths
repo.dir <- '../IBDVerse-sc-eQTL-code/'
data.dir <- paste0(repo.dir,'data/')
coloc.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/results/2025_06_11_IBDverse_coloc_all_gwas/collapsed/'
henry_mr <- 'data/tenk10k_crohns_mr_genes.tsv'
metaf = paste0(data.dir, "eqtl_metadata.txt.gz")
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

# Calculate the novelty of genes in TenK10K vs IBDverse
ibdverse_genes = ibdverse_sub %>%
  pull(phenotype_id) %>%
  unique()

all_sub_de = length(unique(mr$ensembl_gene_id))
all_sub_novel_de = length(setdiff(unique(mr$ensembl_gene_id), ibdverse_genes))

cats = unique(known.coloc.df$category_new)
tenk_novel_all = do.call(rbind, lapply(cats, function(x){
    print(paste0("Testing against: ", x))

    # Get n effector genes (TenK10K vs IBDverse)
    de = mr %>%
        pull(ensembl_gene_id) %>%
        unique()

    nde = length(de)

    # Get novel effector genes for this category
    ibdv_sub = ibdverse_sub %>%
        filter(category_new == x) %>%
        pull(phenotype_id) %>%
        unique()

    novelde = setdiff(de, ibdv_sub)

    # Make contingency
    a <- length(novelde)
    b <- length(setdiff(all_sub_novel_de, novelde))
    c <- length(setdiff(de, novelde))
    d <- length(setdiff(all_sub_de, union(de, all_sub_novel_de)))

    contingency <- matrix(c(a, b, c, d), nrow = 2,
                        dimnames = list(
                            disease_effector = c("in_dataset", "not_in_dataset"),
                            novel_disease_effector   = c("in_dataset", "not_in_dataset")
                        ))

    fishres = fisher.test(contingency)
    or = fishres$estimate
    res = data.frame(category = x, nde = nde, n_novelde = length(novelde), OR=or)
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

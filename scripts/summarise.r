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
pi1f = "results/pi1_all.txt"
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
# Load IBDverse colocs and subset
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

all_sub_de = unique(mr$ensembl_gene_id)
all_sub_novel_de = setdiff(ibdverse_genes, mr_genes)

cats = unique(known.coloc.df$category_new)
cats = cats[-grep("All Cells", cats)]

# Function to compute OR for enrichment of novelty (or non-novelty) of IBDverse disease effector genes vs an external gene set
# ibdverse_df: ibdverse coloc dataframe (filtered to trait of interest)
# external_genes: vector of gene IDs from the external dataset
# ibdverse_all_genes: all IBDverse disease effector genes (universe)
# cat_col: column name in ibdverse_df to group by
# novel: if TRUE, tests enrichment for novelty (not in external); if FALSE, tests enrichment for non-novelty (in external)
# relative: if FALSE (default), anchors the not-novel background to the full external gene set (comparison vs external dataset);
#           if TRUE, restricts background to IBDverse only (comparison relative to other IBDverse categories)
compute_novelty_or <- function(ibdverse_df, external_genes, ibdverse_all_genes, cats, cat_col, novel = TRUE, relative = FALSE, path = NULL, prefix = NULL) {
  do.call(rbind, lapply(cats, function(x) {
    de = ibdverse_df %>%
      filter(!!sym(cat_col) == x) %>%
      pull(phenotype_id) %>%
      unique()
    nde = length(de)

    if (novel) {
      # Question: are cat genes enriched for novelty vs external dataset?
      # a: cat genes not in external (novel, in cat)
      # b: non-cat IBDverse genes not in external (novel, not in cat)
      # c: cat genes in external (not novel, in cat)
      # d: all external genes not in cat (not novel background, anchored to external)
      focal <- setdiff(de, external_genes)
      a <- length(focal)
      b <- length(setdiff(setdiff(ibdverse_all_genes, de), external_genes))
      c <- length(intersect(de, external_genes))
      d <- length(setdiff(external_genes, de))
      n_focal <- length(focal)
      col_name <- "n_novelde"
    } else if (!relative) {
      # Question: are cat genes enriched for non-novelty vs external dataset?
      # a: cat genes in external (not novel, in cat)
      # b: all external genes not in cat (not novel background, anchored to external)
      # c: cat genes not in external (novel, in cat)
      # d: non-cat IBDverse genes not in external (novel background)
      focal <- intersect(de, external_genes)
      a <- length(focal)
      b <- length(setdiff(external_genes, de))
      c <- length(setdiff(de, external_genes))
      d <- length(setdiff(setdiff(ibdverse_all_genes, de), external_genes))
      n_focal <- length(focal)
      col_name <- "n_notnovelde"
    } else {
      # Question: are cat genes relatively enriched for non-novelty vs other IBDverse categories?
      # Universe is IBDverse only — no external-only genes enter the background.
      # a: cat genes in external (not novel, in cat)
      # b: non-cat IBDverse genes in external (not novel background, within IBDverse)
      # c: cat genes not in external (novel, in cat)
      # d: non-cat IBDverse genes not in external (novel background)
      focal <- intersect(de, external_genes)
      a <- length(focal)
      b <- length(intersect(setdiff(ibdverse_all_genes, de), external_genes))
      c <- length(setdiff(de, external_genes))
      d <- length(setdiff(setdiff(ibdverse_all_genes, de), external_genes))
      n_focal <- length(focal)
      col_name <- "n_notnovelde"
    }

    contingency <- matrix(c(a, b, c, d), nrow = 2,
                          dimnames = list(
                            category = c(paste0("in_", x), paste0("not_in_", x)),
                            external_gene_set = c("in_external", "not_in_external")
                          ))

    if (!is.null(path) && !is.null(prefix)) {
      out_f = file.path(path, paste0(prefix, "_", gsub("[^A-Za-z0-9_]", "_", x), "_contingency.tsv"))
      write.table(as.data.frame(contingency), out_f, sep = "\t", quote = FALSE, col.names = NA)
    }

    fishres = fisher.test(contingency)
    res = data.frame(category = x, nde = nde, OR = fishres$estimate,
                     CI_low = fishres$conf.int[1], CI_high = fishres$conf.int[2],
                     pval = fishres$p.value)
    res[[col_name]] <- n_focal
    return(res)
  }))
}

# Are IBDverse cat genes enriched for novelty vs mr_genes?
cont_dir = paste0(out.dir, "/novelty_contingency_tables")
if(!dir.exists(cont_dir)){
    dir.create(cont_dir)
}
mr_novel_res = compute_novelty_or(ibdverse_sub, mr_genes, ibdverse_genes, cats, "category_new", novel = TRUE, relative = FALSE, path = cont_dir, prefix = "mr_novelty_IBDv_vs_TenK10K")

# Append proportion of cells from blood
prop_blood_f = paste0(out.dir, "/prop_blood_ibdverse.csv")
if(!file.exists(prop_blood_f)){
    print("..Loading the big raw meta")
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

    write.csv(prop_blood, paste0(out.dir, "/prop_blood_ibdverse.csv"), row.names = FALSE)

} else {
    print("Loading the pre-computed proportion of blood cells per category")
    prop_blood = read.csv(prop_blood_f)
}

mr_novel_res = mr_novel_res %>%
  left_join(prop_blood, by = c("category" = "predicted_category")) 

mr_novel_res = mr_novel_res %>% 
    filter(is.finite(OR), OR > 0) # MESENCHYMAL has OR of 0

lm_fit <- lm(OR ~ log10(prop_blood), data = mr_novel_res)
slope <- signif(coef(lm_fit)[[2]], 3)
pval <- signif(summary(lm_fit)$coefficients["log10(prop_blood)", "Pr(>|t|)"], 3)
# Print
print(paste0("Slope: ", slope, ". p=", pval))

# Plot this
mr_novel_res$annotation_type = "Major population"
ggplot(mr_novel_res, aes(x = log10(prop_blood), y = OR)) + 
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

ggsave(paste0(out.dir,"/Blood_proportion_vs_novelDE_in_Henry-OR.png"), width = 6.5, height = 5)

# Are IBDverse cat genes enriched for non-novelty vs mr_genes?
tenk_not_novel_res = compute_novelty_or(ibdverse_sub, mr_genes, ibdverse_genes, cats, "category_new", novel = FALSE, relative = FALSE, path = cont_dir, prefix = "mr_non_novelty_IBDv_vs_TenK10K")

# Append the proportion of cells from from each major population which come from blood samples
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

# Are IBDverse cat genes RELATIVELY enriched for non-novelty vs other IBDverse categories?
tenk_not_novel_rel_res = compute_novelty_or(ibdverse_sub, mr_genes, ibdverse_genes, cats, "category_new", novel = FALSE, relative = TRUE, path = cont_dir, prefix = "mr_non_novelty_relative_IBDv_vs_TenK10K")

write.table(tenk_not_novel_rel_res, paste0(out.dir, "/tenk_not_novel_res_relative.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

tenk_not_novel_rel_res = tenk_not_novel_rel_res %>%
  left_join(prop_blood, by = c("category" = "predicted_category")) %>%
  filter(is.finite(OR))

lm_fit <- lm(OR ~ log10(prop_blood), data = tenk_not_novel_rel_res)
slope <- signif(coef(lm_fit)[[2]], 3)
pval <- signif(summary(lm_fit)$coefficients["log10(prop_blood)", "Pr(>|t|)"], 3)
print(paste0("Slope: ", slope, ". p=", pval))

tenk_not_novel_rel_res$annotation_type = "Major population"
ggplot(tenk_not_novel_rel_res, aes(x = log10(prop_blood), y = OR)) +
  geom_point(aes(fill = category, size = annotation_type), shape = 21, stroke = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
  geom_text_repel(aes(label = category), size = 4, max.overlaps = Inf, point.padding = unit(5, "lines")) +
  scale_size_manual(values = annot.class.sizes, name = 'Annotation granularity') +
  scale_fill_manual(values = umap.category.palette, name = 'Major population') +
  labs(
    x = "log10 fraction of cells derived from blood samples",
    y = "Odds ratio for IBDverse disease effector gene already found\nin TenK10K (relative to other categories)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(paste0(out.dir, "/Blood_proportion_vs_notnovelDE_in_Henry-OR_relative.png"), width = 6.5, height = 5)

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

# Save annot mapping
write.table(annot.mapping, paste0(out.dir, "/IBDverse_annotation_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

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

write.table(pi1_combined, paste0(out.dir, "/pi1_combined.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(median_combined, paste0(out.dir, "/pi1_median_combined.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Summarise the max, min and max cell-type per comparison
pi1_combined %>% 
    group_by(comparison) %>% 
    summarise(
        min = min(pi1), 
        max = max(pi1), 
        max_celltype = major_cell_type[which.max(pi1)],
        .groups = "drop"
    )

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

# Plot this as a boxplot instead
ggplot(pi1_combined, aes(x = major_cell_type, y = pi1, fill = comparison)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.7) +
    theme_classic() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)
    ) +
    labs(x = "Major cell-type", y = "Maximum π1 per cell-type", fill = NULL) +
    ylim(c(0,1)) 
  
ggsave(paste0(out.dir, "/pi1_TenK10K_vs_IBDverse_boxplot.png"), width = 10, height = 8)

# ==============================================================================
# Statistical & Ecological Analysis Pipeline (R script)
# Project: Humanisation of Cattle Gut Microbiome
# Author: Angsuman Das
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(tidyverse)
  library(ANCOMBC)
  library(ComplexHeatmap)
  library(circlize)
  library(randomForest)
  library(pROC)
})

# 1. Load Metadata & Abundance Matrices
metadata <- read.csv("metadata_manifest_template.csv", row.names = 1)
# otu_table <- read.table("03_taxonomy/combined_species_counts.tsv", header = TRUE, row.names = 1, sep = "\t")
# ps <- phyloseq(otu_table(as.matrix(otu_table), taxa_are_rows = TRUE), sample_data(metadata))

# 2. Alpha Diversity Analysis
calculate_alpha_diversity <- function(physeq_obj) {
  alpha_df <- estimate_richness(physeq_obj, measures = c("Observed", "Shannon", "Chao1", "Simpson"))
  alpha_df$cohort <- sample_data(physeq_obj)$cohort
  alpha_df$district <- sample_data(physeq_obj)$district
  
  # Wilcoxon rank-sum test between Urban Stray and Rural Domestic
  shannon_test <- wilcox.test(Shannon ~ cohort, data = alpha_df)
  message(sprintf("Shannon Diversity Wilcoxon p-value: %.5e", shannon_test$p.value))
  
  p <- ggplot(alpha_df, aes(x = cohort, y = Shannon, fill = cohort)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2) +
    theme_classic(base_size = 14) +
    labs(title = "Gut Microbiome Alpha Diversity (Shannon Index)",
         x = "Cohort", y = "Shannon Diversity Index") +
    scale_fill_manual(values = c("Urban_Stray" = "#D95F02", "Rural_Domestic" = "#1B9E77"))
  
  return(list(data = alpha_df, plot = p, test = shannon_test))
}

# 3. Beta Diversity & PERMANOVA (Bray-Curtis & Aitchison Distance)
run_beta_diversity <- function(physeq_obj) {
  # Normalize to relative abundance
  ps_rel <- transform_sample_counts(physeq_obj, function(x) x / sum(x))
  bray_dist <- phyloseq::distance(ps_rel, method = "bray")
  
  meta_df <- as(sample_data(ps_rel), "data.frame")
  
  # PERMANOVA controlling for district
  permanova_res <- adonis2(bray_dist ~ cohort + district, data = meta_df, permutations = 999)
  print(permanova_res)
  
  # Ordination
  ord_pcoa <- ordinate(ps_rel, method = "PCoA", distance = bray_dist)
  p_ord <- plot_ordination(ps_rel, ord_pcoa, color = "cohort", shape = "district") +
    stat_ellipse(type = "t", linetype = 2) +
    theme_bw(base_size = 14) +
    labs(title = "PCoA of Gut Microbiome (Bray-Curtis Dissimilarity)") +
    scale_color_manual(values = c("Urban_Stray" = "#D95F02", "Rural_Domestic" = "#1B9E77"))
  
  return(list(permanova = permanova_res, plot = p_ord))
}

# 4. Differential Abundance Analysis using ANCOM-BC2
run_ancombc2 <- function(physeq_obj) {
  out <- ancombc2(
    data = physeq_obj,
    fix_formula = "cohort + district",
    p_adj_method = "holm",
    pseudo_sens = TRUE,
    prv_cut = 0.10,
    lib_cut = 1000,
    group = "cohort",
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 4,
    verbose = TRUE
  )
  return(out)
}

# 5. Machine Learning Biomarker Identification (Random Forest)
train_biomarker_classifier <- function(feature_matrix, cohort_labels) {
  set.seed(42)
  rf_model <- randomForest(x = feature_matrix, y = as.factor(cohort_labels), ntree = 1000, importance = TRUE)
  
  imp <- importance(rf_model, type = 1) # Mean Decrease in Accuracy
  imp_df <- data.frame(Feature = rownames(imp), MeanDecreaseAccuracy = imp[, 1]) %>%
    arrange(desc(MeanDecreaseAccuracy)) %>%
    head(15)
  
  p_imp <- ggplot(imp_df, aes(x = reorder(Feature, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
    geom_col(fill = "#7570B3") +
    coord_flip() +
    theme_classic(base_size = 13) +
    labs(title = "Top 15 Predictive Biomarkers of Humanisation",
         x = "Microbial Taxa / Functional Orthologs", y = "Mean Decrease in Accuracy")
  
  return(list(model = rf_model, importance_plot = p_imp))
}

message("R Ecological and Statistical Analysis Pipeline Loaded.")

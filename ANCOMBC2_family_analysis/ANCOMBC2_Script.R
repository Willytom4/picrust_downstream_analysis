setwd("~/Downloads/Classes_2026/AdvancedMicrobiomeBioinf/Class14")


# ANCOM-BC2 Analysis. Rumen Family-Level Microbial Abundance by RFI Classification.

# 1. Install/load packages

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_bioc <- c("ANCOMBC", "phyloseq", "microbiome")

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

required_cran <- c("tidyverse")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ANCOMBC)
library(phyloseq)
library(microbiome)
library(tidyverse)

#install.packages("ANCOMBC")

# 2. Import data

counts <- read.csv(
  "family_counts_Rumen.csv",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "metadata_Rumen.csv",
  row.names = 1,
  check.names = FALSE
)



# 3. Check sample matching

cat("Number of families:", nrow(counts), "\n")
cat("Number of samples in count table:", ncol(counts), "\n")
cat("Number of samples in metadata:", nrow(metadata), "\n")

if (!all(colnames(counts) %in% rownames(metadata))) {
  missing_samples <- setdiff(colnames(counts), rownames(metadata))
  print(missing_samples)
  stop("Some samples in the count table are missing from the metadata.")
}

metadata <- metadata[colnames(counts), ]

if (!all(colnames(counts) == rownames(metadata))) {
  stop("Sample order still does not match after reordering metadata.")
}

cat("Sample names match correctly.\n")




# 4. Prepare RFI classification

metadata$RFIClass <- factor(
  metadata$RFIClass,
  levels = c("1LowRFI", "2MidRFI", "3HighRFI")
)

cat("RFI class distribution:\n")
print(table(metadata$RFIClass))


# 5. Create phyloseq object

otu <- otu_table(
  as.matrix(counts),
  taxa_are_rows = TRUE
)

sam <- sample_data(metadata)

ps <- phyloseq(otu, sam)



# 6. Run ANCOM-BC2

out <- ancombc2(
  data = ps,
  assay_name = "counts",
  tax_level = NULL,
  fix_formula = "RFIClass",
  rand_formula = NULL,
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  group = "RFIClass",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  global = TRUE,
  pairwise = TRUE,
  dunnet = FALSE,
  trend = FALSE,
  n_cl = 2,
  verbose = TRUE
)



# 7. Inspect output

cat("ANCOM-BC2 output components:\n")
print(names(out))

cat("Main result dimensions:\n")
print(dim(out$res))

cat("Global result dimensions:\n")
print(dim(out$res_global))

cat("Pairwise result dimensions:\n")
print(dim(out$res_pair))



# 8. Save result tables

write.csv(
  out$res,
  "ANCOMBC2_family_full_results.csv",
  row.names = FALSE
)

write.csv(
  out$res_global,
  "ANCOMBC2_family_global_test_FULL.csv",
  row.names = FALSE
)

write.csv(
  out$res_pair,
  "ANCOMBC2_family_pairwise_tests_FULL.csv",
  row.names = FALSE
)



# 9. Extract main result columns

res <- out$res

lfc_cols <- grep("^lfc", names(res), value = TRUE)
se_cols <- grep("^se", names(res), value = TRUE)
w_cols <- grep("^W", names(res), value = TRUE)
p_cols <- grep("^p", names(res), value = TRUE)
q_cols <- grep("^q", names(res), value = TRUE)
diff_cols <- grep("^diff", names(res), value = TRUE)

cat("Main q-value columns found:\n")
print(q_cols)

lfc_table <- res[, c("taxon", lfc_cols), drop = FALSE]
q_table <- res[, c("taxon", q_cols), drop = FALSE]
diff_table <- res[, c("taxon", diff_cols), drop = FALSE]

write.csv(
  lfc_table,
  "ANCOMBC2_family_log_fold_changes.csv",
  row.names = FALSE
)

write.csv(
  q_table,
  "ANCOMBC2_family_q_values.csv",
  row.names = FALSE
)

write.csv(
  diff_table,
  "ANCOMBC2_family_significant_results.csv",
  row.names = FALSE
)

alpha <- 0.05

if (length(q_cols) > 0) {
  
  significant_main <- res[
    rowSums(res[, q_cols, drop = FALSE] < alpha, na.rm = TRUE) > 0,
  ]
  
  write.csv(
    significant_main,
    "ANCOMBC2_family_main_significant_only.csv",
    row.names = FALSE
  )
  
  cat("Number of significant families in main results:",
      nrow(significant_main), "\n")
  
} else {
  
  warning("No q-value columns found in out$res.")
  
}



# 10. Save results

res_global <- out$res_global

global_q_cols <- grep("^q", names(res_global), value = TRUE)

cat("Global q-value columns found:\n")
print(global_q_cols)

if (length(global_q_cols) > 0) {
  
  significant_global <- res_global[
    rowSums(res_global[, global_q_cols, drop = FALSE] < alpha, na.rm = TRUE) > 0,
  ]
  
  write.csv(
    significant_global,
    "ANCOMBC2_family_global_significant_only.csv",
    row.names = FALSE
  )
  
  cat("Number of significant families in global test:",
      nrow(significant_global), "\n")
  
} else {
  
  warning("No q-value columns found in out$res_global.")
  
}

res_pair <- out$res_pair

pair_q_cols <- grep("^q", names(res_pair), value = TRUE)

cat("Pairwise q-value columns found:\n")
print(pair_q_cols)

if (length(pair_q_cols) > 0) {
  
  significant_pairwise <- res_pair[
    rowSums(res_pair[, pair_q_cols, drop = FALSE] < alpha, na.rm = TRUE) > 0,
  ]
  
  write.csv(
    significant_pairwise,
    "ANCOMBC2_family_pairwise_significant_only.csv",
    row.names = FALSE
  )
  
  cat("Number of significant families in pairwise tests:",
      nrow(significant_pairwise), "\n")
  
} else {
  
  warning("No q-value columns found in out$res_pair.")
  
}



# 11. Interpretation guide

cat("\nInterpretation guide:\n")
cat("RFIClass2MidRFI = 2MidRFI vs 1LowRFI\n")
cat("RFIClass3HighRFI = 3HighRFI vs 1LowRFI\n")
cat("Positive lfc = greater abundance compared with 1LowRFI\n")
cat("Negative lfc = lower abundance compared with 1LowRFI\n")
cat("q-value < 0.05 = significant after FDR correction\n")
cat("Global test = any difference among 1LowRFI, 2MidRFI, and 3HighRFI\n")
cat("Pairwise test = pairwise comparisons among RFI classes\n")
#!/usr/bin/env Rscript

# ============================================================
# Metabolomics Batch Effect Correction and sPCA Visualization
# ============================================================
# This script loads metabolomics data, performs log transformation,
# applies PCA and batch correction (ComBat), visualizes sPCA plots,
# and identifies outliers before and after correction.
# ============================================================

# Load Required Libraries
library(dplyr)
library(dbnorm)
library(readxl)
library(mixOmics)
library(sva)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# Load Data (Update file paths as needed)
# ------------------------------------------------------------
# Set working directory
# setwd("your/path/here")  # <- Customize this for public use or remove

# Read input data
TASOAC <- read_excel("TASOAC_dataset.xlsx")
Licofelone <- read_excel("Licofelone_dataset.xlsx")

# ------------------------------------------------------------
# Preprocess Data (Log2 transform, batch info)
# ------------------------------------------------------------
data2 <- TASOAC[6:109]
data3 <- log2(data2)
batch <- TASOAC[,1]
data4 <- cbind(batch, data3)

# ------------------------------------------------------------
# Tune PCA and Perform sPCA (Before Correction)
# ------------------------------------------------------------
explainedVariance <- tune.pca(data3, ncomp = 10, center = TRUE, scale = TRUE)
plot(explainedVariance)  # Choose ncomp visually

trans.spca <- spca(data3, ncomp = 2, center = TRUE, scale = TRUE)

scores_before <- as.data.frame(trans.spca$variates$X)
scores_before$BatchID <- as.factor(TASOAC$BatchID)

expl_var_before <- trans.spca$prop_expl_var$X * 100
xlab_before <- paste0("PC1: ", round(expl_var_before[1], 1), "% expl. var")
ylab_before <- paste0("PC2: ", round(expl_var_before[2], 1), "% expl. var")

plot_a <- ggplot(scores_before, aes(x = PC1, y = PC2, color = BatchID, label = rownames(scores_before))) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = BatchID), level = 0.99) +
  geom_text(nudge_x = 0.3, nudge_y = 0.3, size = 3, check_overlap = TRUE) +
  labs(x = xlab_before, y = ylab_before) +
  theme_minimal() +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")

# ------------------------------------------------------------
# Normality Check (Shapiro-Wilk Test)
# ------------------------------------------------------------
shapiro_tests <- apply(data3, 2, function(x) shapiro.test(x)$p.value)
shapiro_results <- data.frame(Feature = colnames(data3), P_Value = shapiro_tests)
non_normal_features <- shapiro_results[shapiro_results$P_Value < 0.05, ]

print(shapiro_results)
print(non_normal_features)

# ------------------------------------------------------------
# Batch Effect Correction (ComBat)
# ------------------------------------------------------------
data5 <- ComBat(dat = t(data3), batch = data4$BatchID, mod = NULL, par.prior = FALSE)
data5 <- as.data.frame(t(data5))

# ------------------------------------------------------------
# sPCA (After Batch Correction)
# ------------------------------------------------------------
trans.spca_batch_corrected <- spca(data5, ncomp = 2, center = TRUE, scale = TRUE)

scores_after <- as.data.frame(trans.spca_batch_corrected$variates$X)
scores_after$BatchID <- as.factor(data4$BatchID)

expl_var_after <- trans.spca_batch_corrected$prop_expl_var$X * 100
xlab_after <- paste0("PC1: ", round(expl_var_after[1], 1), "% expl. var")
ylab_after <- paste0("PC2: ", round(expl_var_after[2], 1), "% expl. var")

plot_b <- ggplot(scores_after, aes(x = PC1, y = PC2, color = BatchID, label = rownames(scores_after))) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = BatchID), level = 0.999) +
  geom_text(nudge_x = 0.3, nudge_y = 0.3, size = 3, check_overlap = TRUE) +
  labs(x = xlab_after, y = ylab_after) +
  theme_minimal() +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")

# Combine plots
combined_plot <- plot_a + plot_b + plot_layout(ncol = 2)
print(combined_plot)

# Save figure
ggsave("sPCA_before_and_after_batch_correction.png", plot = combined_plot,
       width = 16, height = 6, units = "in", dpi = 300)

# ------------------------------------------------------------
# Outlier Detection (Before and After Batch Correction)
# ------------------------------------------------------------

# Before correction
z_scores_before <- scale(trans.spca$variates$X)
outliers_before <- which(abs(z_scores_before) > 3, arr.ind = TRUE)
if (nrow(outliers_before) > 0) {
  cat("Outliers BEFORE correction:\n")
  print(rownames(z_scores_before)[unique(outliers_before[,1])])
  print(data.frame(
    Sample = rownames(z_scores_before)[outliers_before[,1]],
    Component = colnames(z_scores_before)[outliers_before[,2]],
    Z_score = z_scores_before[outliers_before]
  ))
} else {
  cat("No outliers detected BEFORE correction.\n")
}

# After correction
z_scores_after <- scale(trans.spca_batch_corrected$variates$X)
outliers_after <- which(abs(z_scores_after) > 3, arr.ind = TRUE)
if (nrow(outliers_after) > 0) {
  cat("Outliers AFTER correction:\n")
  print(rownames(z_scores_after)[unique(outliers_after[,1])])
  print(data.frame(
    Sample = rownames(z_scores_after)[outliers_after[,1]],
    Component = colnames(z_scores_after)[outliers_after[,2]],
    Z_score = z_scores_after[outliers_after]
  ))
} else {
  cat("No outliers detected AFTER correction.\n")
}

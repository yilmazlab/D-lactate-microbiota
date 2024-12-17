#!/bin/bash
# Title: Alpha and Beta Diversity Analysis
# Description: This script performs alpha and beta diversity analyses on microbial data imported from QIIME 2 into phyloseq.
# Usage: Execute this script in a Linux environment with R installed or in a managed R environment.

# Step 1: Activate R environment (optional)
# Uncomment the following line if using conda
# conda activate r_env

# Step 2: Execute R script
Rscript <<EOF
# Load required libraries
library(phyloseq)
library(ggplot2)
library(vegan)

# Step 1: Import QIIME 2 data into phyloseq
phyloseq_object <- qza_to_phyloseq(
  features = "table.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

# Step 2: Filter samples with at least 1000 reads
phyloseq_object <- prune_samples(sample_sums(phyloseq_object) >= 1000, phyloseq_object)

# Step 3: Alpha diversity analysis
# Plot Shannon and Simpson indices
alpha_plot <- plot_richness(phyloseq_object, x = "Group", measures = c("Shannon", "Simpson")) +
  theme_bw() +
  ggtitle("Alpha Diversity")
ggsave("alpha_diversity_plot.pdf", alpha_plot, width = 8, height = 6)

# Perform Wilcoxon test for alpha diversity
alpha_data <- estimate_richness(phyloseq_object, measures = c("Shannon", "Simpson"))
alpha_test <- pairwise.wilcox.test(alpha_data$Shannon, sample_data(phyloseq_object)$Group)
write.table(alpha_test$p.value, file = "alpha_diversity_stats.txt", sep = "\t", col.names = NA)
print(alpha_test)

# Step 4: Beta diversity analysis
# Calculate Bray-Curtis distances
bray_dist <- phyloseq::distance(phyloseq_object, method = "bray")

# Perform PCoA
ordination <- ordinate(phyloseq_object, method = "PCoA", distance = bray_dist)

# Plot PCoA
pcoa_plot <- plot_ordination(phyloseq_object, ordination, color = "Group") +
  geom_point(size = 4) +
  theme_bw() +
  ggtitle("Beta Diversity (PCoA)")
ggsave("beta_diversity_pcoa_plot.pdf", pcoa_plot, width = 8, height = 6)

# Perform PERMANOVA (Adonis) on beta diversity
adonis_results <- adonis(bray_dist ~ Group, data = as(sample_data(phyloseq_object), "data.frame"))
write.table(adonis_results$aov.tab, file = "beta_diversity_adonis_results.txt", sep = "\t", col.names = NA)
print(adonis_results)
EOF
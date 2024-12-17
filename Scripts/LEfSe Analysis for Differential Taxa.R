#!/bin/bash
# Title: LEfSe Analysis for Differential Taxa
# Description: This script performs LEfSe analysis on genus-level microbial data to identify taxa associated with responders and non-responders.
# Usage: Run this script with R installed or within a managed R environment.

# Step 1: Activate R environment (optional, uncomment if needed)
# conda activate r_env

# Step 2: Execute R script
Rscript <<EOF
# Load necessary libraries
library(phyloseq)
library(microbiomeMarker)
library(ggplot2)

# Step 1: Prepare data for LEfSe analysis
# Aggregate taxa to genus level
data_genus <- tax_glom(physeq_data, taxrank = "Genus")

# Add a response group column for LEfSe analysis
sample_data(data_genus)$Response_Group <- ifelse(
  sample_data(data_genus)$Group == "Treatment", "Responders", "NonResponders"
)

# Step 2: Run LEfSe analysis
lefse_results <- run_lefse(
  data_genus,
  group = "Response_Group", # Grouping variable
  taxa_rank = "Genus",      # Taxonomic level
  lda_cutoff = 2.0,         # Minimum LDA score
  wilcoxon_cutoff = 0.05    # Maximum p-value
)

# Step 3: Extract and visualize results
# Convert LEfSe results to a data frame
marker_table <- as.data.frame(lefse_results@marker_table)

# Save the marker table
write.csv(marker_table, file = "lefse_marker_table.csv", row.names = FALSE)

# Create and save a bar plot of LDA scores
lda_plot <- ggplot(marker_table, aes(x = reorder(feature, ef_lda), y = ef_lda, fill = enrich_group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "LEfSe Results", x = "Genera", y = "LDA Score") +
  theme_minimal() +
  scale_fill_manual(values = c("Responders" = "#1f78b4", "NonResponders" = "#e31a1c"))
ggsave("lefse_lda_plot.pdf", lda_plot, width = 10, height = 7)
EOF
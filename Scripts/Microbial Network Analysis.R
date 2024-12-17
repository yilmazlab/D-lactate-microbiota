#!/bin/bash
# Title: Microbial Network Analysis
# Description: This script analyzes microbial interactions in gut microbiota data at the genus level and calculates network metrics.
# Usage: Execute this script in an R-enabled environment.

# Step 1: Activate R environment (optional)
# Uncomment the following line if using conda or similar environments
# conda activate r_env

# Step 2: Execute R script
Rscript <<EOF
# Load necessary libraries
library(phyloseq)
library(Hmisc)
library(igraph)

# Step 1: Define the microbial network analysis function
analyze_microbial_network <- function(physeq, group, cor_threshold = 0.7, p_threshold = 0.05) {
  # Subset samples for the specified group
  subset_physeq <- prune_samples(sample_data(physeq)$Group == group, physeq)
  abundance_mat <- as.matrix(t(otu_table(subset_physeq)))

  # Filter taxa by genus assignment and prevalence
  genus_names <- tax_table(physeq)[, "Genus"]
  abundance_mat <- abundance_mat[, !is.na(genus_names) & apply(abundance_mat > 0, 2, mean) >= 0.2]
  genus_names <- genus_names[!is.na(genus_names) & apply(abundance_mat > 0, 2, mean) >= 0.2]

  # Calculate correlations and create an adjacency matrix
  cor_result <- rcorr(abundance_mat, type = "spearman")
  adj_matrix <- cor_result$r
  adj_matrix[abs(adj_matrix) < cor_threshold | cor_result$P > p_threshold] <- 0

  # Generate a graph with genus names and compute network metrics
  graph <- graph.adjacency(abs(adj_matrix), mode = "undirected", weighted = TRUE)
  V(graph)$name <- genus_names
  metrics_df <- data.frame(
    Genus = genus_names,
    Degree = degree(graph),
    Betweenness = betweenness(graph, weights = abs(E(graph)$weight)),
    Closeness = closeness(graph, weights = abs(E(graph)$weight)),
    EigenCentrality = eigen_centrality(graph, weights = abs(E(graph)$weight))$vector
  )
  return(list(graph = graph, node_metrics = metrics_df, abundance_mat = abundance_mat))
}

# Step 2: Example usage of the analyze_microbial_network function
# Load phyloseq object (replace 'your_phyloseq_object.rds' with your file path)
physeq <- readRDS("your_phyloseq_object.rds")

# Specify the group of interest (e.g., "Group A")
network_results <- analyze_microbial_network(physeq, group = "Group A")

# Save network metrics to a CSV file
write.csv(network_results$node_metrics, file = "network_metrics.csv", row.names = FALSE)

# Save the graph as a GraphML file for visualization
write_graph(network_results$graph, file = "microbial_network.graphml", format = "graphml")
EOF
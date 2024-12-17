#!/bin/bash
# Title: Random Forest Analysis for Microbial Data
# Description: This script runs a Random Forest analysis using R to identify microbial signatures associated with responsiveness at specific time points.
# Usage: Run this script in a Linux environment with R installed.

# Step 1: Activate R environment (if using a managed environment like conda)
# Uncomment the following line if R is part of a conda environment
# conda activate r_env

# Step 2: Execute the R script
Rscript <<EOF
# Load necessary libraries
library(randomForest)
library(caret)
library(pROC)
library(phyloseq)
library(ggplot2)

# Define Random Forest analysis function
perform_rf_analysis <- function(physeq_obj, time_point) {
  # Subset data for the specified time point
  physeq_filtered <- subset_samples(physeq_obj, Time == time_point)
  physeq_genus <- tax_glom(physeq_filtered, taxrank = "Genus")
  
  # Prepare feature (X) and response (y) variables
  otu_table_df <- as.data.frame(otu_table(physeq_genus))
  sample_data_df <- data.frame(sample_data(physeq_filtered))
  X <- t(otu_table_df)  # Features (OTU/genus table)
  y <- factor(sample_data_df$Responsive, levels = c("non_responsive", "responsive"))
  
  # Split data into training and test sets
  set.seed(123)
  train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
  X_train <- X[train_idx, ]
  X_test <- X[-train_idx, ]
  y_train <- y[train_idx]
  y_test <- y[-train_idx]
  
  # Train Random Forest model
  rf_model <- randomForest(X_train, y_train, importance = TRUE, ntree = 500)
  
  # Evaluate model performance
  predictions <- predict(rf_model, X_test, type = "prob")
  roc_curve <- roc(response = y_test, predictor = predictions[, 2])
  
  # Save model and performance metrics to files
  saveRDS(rf_model, file = paste0("rf_model_time", time_point, ".rds"))
  pdf(file = paste0("roc_curve_time", time_point, ".pdf"))
  plot(roc_curve, main = paste("ROC Curve for Time", time_point))
  dev.off()
  
  # Return results
  list(model = rf_model, roc_curve = roc_curve)
}

# Example usage for Time 0
# Replace 'Exp_mice' with your actual phyloseq object
results_time0 <- perform_rf_analysis(physeq_obj = Exp_mice, time_point = 0)

# Example usage for Time 10
results_time10 <- perform_rf_analysis(physeq_obj = Exp_mice, time_point = 10)
EOF
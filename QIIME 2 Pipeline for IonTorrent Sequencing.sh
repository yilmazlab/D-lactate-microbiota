#!/bin/bash
# Title: QIIME 2 Pipeline for IonTorrent Sequencing
# Description: This script processes raw FASTQ files using QIIME 2-2021.4, including importing, demultiplexing, quality control, denoising, and taxonomic classification.
# Usage: Run on UBELIX Linux cluster with QIIME 2-2021.4 environment activated.

# Activate QIIME 2 environment
conda activate qiime2-2021.4

# Step 1: Import raw FASTQ data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path example_data/reads.fastq.gz \
  --output-path demux.qza

# Step 2: Summarize demultiplexed data
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

# Step 3: Denoise and filter sequences using DADA2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trunc-len 120 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

# Step 4: Summarize feature table
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv

# Step 5: Taxonomic classification using pre-trained Silva classifier
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

# Step 6: Generate taxonomy bar plot
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --o-visualization taxa-bar-plots.qzv
# load packages needed for Aim 2 setup
library(tidyverse)
library(phyloseq)
library(ape)

# load the phyloseq object saved from Aim 1
ps <- readRDS("../data/ps_raw.rds")

# print the phyloseq object and basic dimensions
ps
cat("Samples:", nsamples(ps), "Taxa:", ntaxa(ps), "\n")

# keep only gastric cancer samples with intestinal or diffused Lauren subtype
ps2 <- subset_samples(
  ps,
  Group == "Gastric cancer (GC)" &
    `Histopathology.......lauren.classification.` %in% c("Intestinal type", "Diffused type")
)

# remove taxa with zero counts after subsetting samples
ps2 <- prune_taxa(taxa_sums(ps2) > 0, ps2)

# create a simpler subtype variable for downstream analysis
sample_data(ps2)$Subtype <- ifelse(
  sample_data(ps2)$`Histopathology.......lauren.classification.` == "Intestinal type",
  "Intestinal",
  "Diffuse"
)

# set subtype as a factor with intestinal as the reference level
sample_data(ps2)$Subtype <- factor(
  sample_data(ps2)$Subtype,
  levels = c("Intestinal", "Diffuse")
)

# check subtype counts and subset dimensions
table(sample_data(ps2)$Subtype)
cat("Subset samples:", nsamples(ps2), "Taxa:", ntaxa(ps2), "\n")


# load additional packages needed for Aim 2 analyses and plotting
library(ANCOMBC)
library(indicspecies)
library(pheatmap)
library(ggrepel)
library(scales)
library(vegan)

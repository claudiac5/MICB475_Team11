# Loads in relevant packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(ape)

# Designates file paths
feature_fp <- "feature_table_gc_intestinal_diffuse.tsv"
tax_fp     <- "taxonomy.tsv"
meta_fp    <- "metadata.tsv"
tree_fp    <- "tree.nwk"

# Creates feature table 
otu2 <- read.delim(
  feature_fp,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  quote = ""
)

# Creates Metadata variables (using sample-id as rownames)
meta2 <- read.delim(meta_fp, header = TRUE, sep = "\t", check.names = FALSE)
rownames(meta2) <- meta2[["sample-id"]]
meta2 <- meta2[, setdiff(names(meta2), "sample-id"), drop = FALSE]

# Taxonomy (split Taxon into ranks)
tax <- read_delim(tax_fp, delim = "\t", show_col_types = FALSE)

tax_clean <- tax %>%
  rename(FeatureID = 1) %>%
  select(-any_of("Confidence"))

tax_mat <- tax_clean %>%
  separate(
    col = Taxon,
    into = c("Domain","Phylum","Class","Order","Family","Genus","Species"),
    sep = ";\\s*",
    fill = "right"
  ) %>%
  mutate(across(-FeatureID, ~ gsub("^D_\\d+__|^[a-z]__", "", .x))) %>%
  tibble::column_to_rownames("FeatureID") %>%
  as.matrix()

# Reads in tree
tree <- read.tree(tree_fp)

# Build phyloseq pieces
OTU  <- otu_table(as.matrix(otu2), taxa_are_rows = TRUE)
SAMP <- sample_data(meta2)
TAX  <- tax_table(tax_mat)

# Build phyloseq object
ps <- phyloseq(OTU, SAMP, TAX, phy_tree(tree))

# Saves phyloseq file
ps

cat("Samples:", nsamples(ps), "Taxa:", ntaxa(ps), "\n")

saveRDS(ps, "ps_raw.rds")

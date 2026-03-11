# load packages needed for Aim 2 setup and analysis
library(tidyverse)
library(phyloseq)
library(ape)
library(ANCOMBC)
library(indicspecies)
library(pheatmap)
library(ggrepel)
library(scales)
library(vegan)

# create results and figures folders if they do not already exist
if (!dir.exists("../results")) dir.create("../results", recursive = TRUE)
if (!dir.exists("../figures")) dir.create("../figures", recursive = TRUE)

# designates file paths
feature_fp <- "feature_table_gc_intestinal_diffuse.tsv"
tax_fp     <- "taxonomy.tsv"
meta_fp    <- "metadata.tsv"
tree_fp    <- "tree.nwk"

# creates feature table 
otu2 <- read.delim(
  feature_fp,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  quote = ""
)

# creates Metadata variables (using sample-id as rownames)
meta2 <- read.delim(meta_fp, header = TRUE, sep = "\t", check.names = FALSE)
rownames(meta2) <- meta2[["sample-id"]]
meta2 <- meta2[, setdiff(names(meta2), "sample-id"), drop = FALSE]

# taxonomy (split Taxon into ranks)
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

# reads in tree
tree <- read.tree(tree_fp)

# build phyloseq pieces
OTU  <- otu_table(as.matrix(otu2), taxa_are_rows = TRUE)
SAMP <- sample_data(meta2)
TAX  <- tax_table(tax_mat)

# build phyloseq object
ps <- phyloseq(OTU, SAMP, TAX, phy_tree(tree))

# saves phyloseq file
ps

cat("Samples:", nsamples(ps), "Taxa:", ntaxa(ps), "\n")

saveRDS(ps, "ps_raw.rds")

# load the phyloseq object saved from Aim 1
ps <- readRDS("ps_raw.rds")

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

# inspect taxonomy for possible mitochondria and chloroplast labels
tax_df_check <- as.data.frame(tax_table(ps2))
unique(tax_df_check$Class)
unique(tax_df_check$Order)
unique(tax_df_check$Family)

# remove taxa annotated as mitochondria or chloroplast
tax_df_ps2 <- as.data.frame(tax_table(ps2))

keep_taxa_mc <- !(
  grepl("mitochondria", tax_df_ps2$Class, ignore.case = TRUE, perl = TRUE) |
    grepl("mitochondria", tax_df_ps2$Order, ignore.case = TRUE, perl = TRUE) |
    grepl("mitochondria", tax_df_ps2$Family, ignore.case = TRUE, perl = TRUE) |
    grepl("chloroplast", tax_df_ps2$Class, ignore.case = TRUE, perl = TRUE) |
    grepl("chloroplast", tax_df_ps2$Order, ignore.case = TRUE, perl = TRUE) |
    grepl("chloroplast", tax_df_ps2$Family, ignore.case = TRUE, perl = TRUE)
)

keep_taxa_mc[is.na(keep_taxa_mc)] <- TRUE
keep_taxa_mc[is.na(keep_taxa_mc)] <- TRUE

keep_taxa_mc[is.na(keep_taxa_mc)] <- TRUE

ps2 <- prune_taxa(keep_taxa_mc, ps2)

# remove taxa with zero counts after mitochondria/chloroplast filtering
ps2 <- prune_taxa(taxa_sums(ps2) > 0, ps2)

# print dimensions after contaminant filtering
cat("After mitochondria/chloroplast filtering - Samples:", nsamples(ps2), "Taxa:", ntaxa(ps2), "\n")

# save the filtered Aim 2 subset
saveRDS(ps2, "../results/ps_aim2_subset.rds")

# agglomerate taxa to the genus level
ps_genus <- tax_glom(ps2, taxrank = "Genus", NArm = FALSE)

# convert taxonomy table to a data frame and clean genus labels
tax_df <- as.data.frame(tax_table(ps_genus))
tax_df$Genus <- as.character(tax_df$Genus)
tax_df$Family <- as.character(tax_df$Family)

# replace missing genus labels with family-based placeholder labels
missing_genus <- is.na(tax_df$Genus) | tax_df$Genus == ""
tax_df$Genus[missing_genus] <- paste0("Unclassified_", tax_df$Family[missing_genus])

# put the cleaned taxonomy table back into the phyloseq object
tax_table(ps_genus) <- tax_table(as.matrix(tax_df))

# print and save the genus-level phyloseq object
cat("Genus-level taxa:", ntaxa(ps_genus), "\n")
saveRDS(ps_genus, "../results/ps_genus.rds")

# extract the OTU table as a matrix for prevalence filtering
otu_mat <- as(otu_table(ps_genus), "matrix")

# count the number of samples in which each genus is present
if (taxa_are_rows(ps_genus)) {
  prev_counts <- apply(otu_mat, 1, function(x) sum(x > 0))
} else {
  prev_counts <- apply(otu_mat, 2, function(x) sum(x > 0))
}

# define the minimum prevalence cutoff as 10 percent of samples
min_prev <- ceiling(0.10 * nsamples(ps_genus))
cat("Minimum prevalence cutoff:", min_prev, "samples\n")

# keep only genera present in at least the minimum number of samples
keep_taxa_prev <- names(prev_counts[prev_counts >= min_prev])

# prune low-prevalence genera from the phyloseq object
ps_genus_filt <- prune_taxa(keep_taxa_prev, ps_genus)

# print and save the filtered genus-level object
cat("Filtered genus-level taxa:", ntaxa(ps_genus_filt), "\n")
saveRDS(ps_genus_filt, "../results/ps_genus_filt.rds")

# create a relative abundance version for plotting only
ps_genus_rel <- transform_sample_counts(ps_genus_filt, function(x) x / sum(x))

# save the relative-abundance object
saveRDS(ps_genus_rel, "../results/ps_genus_rel.rds")

# make sure subtype is set correctly for downstream analyses
sample_data(ps_genus_filt)$Subtype <- factor(
  sample_data(ps_genus_filt)$Subtype,
  levels = c("Intestinal", "Diffuse")
)

# ANCOM-BC2 differential abundance analysis

# run ANCOM-BC2 on the filtered genus-level count data
ancom_out <- ancombc2(
  data = ps_genus_filt,
  tax_level = NULL,
  fix_formula = "Subtype",
  rand_formula = NULL,
  group = "Subtype",
  p_adj_method = "BH",
  prv_cut = 0,
  lib_cut = 0,
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 1,
  verbose = TRUE,
  global = FALSE,
  pairwise = TRUE,
  dunnet = FALSE,
  trend = FALSE
)

# inspect the ANCOM-BC2 output structure
names(ancom_out)
colnames(ancom_out$res)

# manually extract the taxonomy table from the phyloseq object as a dataframe
tax_map_df <- as.data.frame(as(tax_table(ps_genus_filt), "matrix")) %>%
  rownames_to_column("taxon") %>%
  select(taxon, Genus)

ancom_res <- as.data.frame(ancom_out$res) %>%
  left_join(tax_map_df, by = "taxon") %>%
  mutate(
    Genus = ifelse(is.na(Genus) | Genus == "", taxon, Genus),
    Genus = make.unique(as.character(Genus))
  )

# identify the relevant output columns automatically
lfc_col  <- grep("^lfc_.*Diffuse", colnames(ancom_res), value = TRUE)[1]
q_col    <- grep("^q_.*Diffuse", colnames(ancom_res), value = TRUE)[1]
p_col    <- grep("^p_.*Diffuse", colnames(ancom_res), value = TRUE)[1]
diff_col <- grep("^diff_.*Diffuse", colnames(ancom_res), value = TRUE)[1]
se_col   <- grep("^se_.*Diffuse", colnames(ancom_res), value = TRUE)[1]

# create a clean ANCOM-BC2 results table
ancom_tbl <- ancom_res %>%
  transmute(
    Genus = Genus,
    logFC = .data[[lfc_col]],
    SE = if (!is.na(se_col)) .data[[se_col]] else NA_real_,
    p = .data[[p_col]],
    q = .data[[q_col]],
    significant = .data[[diff_col]]
  ) %>%
  mutate(
    enriched_in = case_when(
      significant & logFC > 0 ~ "Diffuse",
      significant & logFC < 0 ~ "Intestinal",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(q)

# create a filtered table of significant ANCOM-BC2 genera
ancom_sig <- ancom_tbl %>%
  filter(significant == TRUE, q < 0.05)

sig_genera <- ancom_sig$Genus

# save ANCOM-BC2 results
write.csv(ancom_tbl, "../results/aim2_ancombc2_all_genera.csv", row.names = FALSE)
write.csv(ancom_sig, "../results/aim2_ancombc2_significant_genera.csv", row.names = FALSE)

# print the number of significant ANCOM-BC2 genera
cat("Number of significant ANCOM-BC2 genera:", nrow(ancom_sig), "\n")

# Indicator Species Analysis

# create a temporary object 
ps_temp <- ps_genus_rel

# rename taxa to Genus names, but use make.unique to prevent errors if 
# multiple ASVs have the same Genus name
taxa_names(ps_temp) <- make.unique(as.character(as.data.frame(tax_table(ps_temp))$Genus))

# convert to matrix (the column names will now be the unique Genus names)
ind_mat <- as(otu_table(ps_temp), "matrix")

if (taxa_are_rows(ps_temp)) {
  ind_mat <- t(ind_mat)
}

# extract subtype labels for indicator species analysis
group_vec <- sample_data(ps_genus_rel)$Subtype %>% as.character()

# run indicator species analysis using subtype groups
ind_out <- multipatt(
  x = ind_mat,
  cluster = group_vec,
  func = "IndVal.g",
  duleg = TRUE,
  control = how(nperm = 9999)
)

# extract indicator species results and adjust p-values with BH correction
ind_sign <- as.data.frame(ind_out$sign) %>%
  rownames_to_column("Genus") %>%
  mutate(q = p.adjust(p.value, method = "BH"))

# map indicator combination indices back to subtype labels
comb_df <- as.data.frame(ind_out$comb)
comb_df$index <- seq_len(nrow(comb_df))
grp_cols <- setdiff(colnames(comb_df), "index")
comb_df$associated_group <- apply(
  comb_df[, grp_cols, drop = FALSE],
  1,
  function(x) paste(names(x)[x == 1], collapse = "+")
)

# create clean indicator species result tables
ind_tbl <- ind_sign %>%
  left_join(comb_df[, c("index", "associated_group")], by = "index") %>%
  arrange(q)

ind_sig <- ind_tbl %>%
  filter(q < 0.05)

# save indicator species results
write.csv(ind_tbl, "../results/aim2_indicator_species_all_genera.csv", row.names = FALSE)
write.csv(ind_sig, "../results/aim2_indicator_species_significant_genera.csv", row.names = FALSE)

# print the number of significant indicator genera
cat("Number of significant indicator genera:", nrow(ind_sig), "\n")

# Combined summary table

# combine ANCOM-BC2 and indicator species results into one summary table
summary_tbl <- full_join(
  ancom_tbl %>% select(Genus, logFC, SE, p, q, significant, enriched_in),
  ind_tbl %>% select(Genus, stat, p.value, q, associated_group),
  by = "Genus",
  suffix = c("_ANCOM", "_IndVal")
)

# save the combined summary table
write.csv(summary_tbl, "../results/aim2_combined_summary.csv", row.names = FALSE)

# Barplot of top genera by subtype

# melt the relative-abundance genus-level object for plotting
plot_bar_df <- psmelt(ps_genus_rel) %>%
  group_by(Subtype, Genus) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop")

# identify the top 15 genera by overall mean relative abundance
top_genera <- plot_bar_df %>%
  group_by(Genus) %>%
  summarise(overall_mean = mean(mean_abundance), .groups = "drop") %>%
  arrange(desc(overall_mean)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

# keep only the top genera for the barplot
plot_bar_top <- plot_bar_df %>%
  filter(Genus %in% top_genera)

# create the stacked barplot of top genera by subtype
p_bar <- ggplot(plot_bar_top, aes(x = Subtype, y = mean_abundance, fill = Genus)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Top 15 genera by subtype",
    x = "Subtype",
    y = "Relative abundance"
  ) +
  theme_bw()

# save the barplot
ggsave("../figures/aim2_top15_genera_barplot.png", p_bar, width = 8, height = 6, dpi = 300)

# ===========================
# Heatmap of Significant Genera with q-value stars
# ===========================

if(length(sig_genera) > 0) {
  
  # 1. Ensure OTU names in phyloseq match Genus
  taxa_names(ps_genus_rel) <- make.unique(as.character(as.data.frame(tax_table(ps_genus_rel))$Genus))
  
  # 2. Keep only significant genera present in phyloseq
  sig_genera_checked <- sig_genera[sig_genera %in% taxa_names(ps_genus_rel)]
  
  if(length(sig_genera_checked) == 0) stop("No significant genera found in the phyloseq object.")
  
  # 3. Melt phyloseq and filter
  heat_df_base <- psmelt(ps_genus_rel) %>%
    filter(OTU %in% sig_genera_checked)
  
  # 4. Combine q-values from ANCOM and Indicator Species, take minimum
  q_df <- data.frame(OTU = sig_genera_checked) %>%
    left_join(ancom_tbl %>% select(OTU = Genus, q_ancom = q), by = "OTU") %>%
    left_join(ind_tbl %>% select(OTU = Genus, q_ind = q), by = "OTU") %>%
    rowwise() %>%
    mutate(
      best_q = min(c(q_ancom, q_ind), na.rm = TRUE),
      sig_stars = case_when(
        best_q < 0.001 ~ "***",
        best_q < 0.01  ~ "**",
        best_q < 0.05  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    ungroup() %>%
    select(OTU, sig_stars)
  
  # 5. Add stars to Genus names
  heat_df <- heat_df_base %>%
    left_join(q_df, by = "OTU") %>%
    mutate(
      sig_stars = ifelse(is.na(sig_stars), "", sig_stars),
      Genus_labeled = paste0(Genus, sig_stars)
    ) %>%
    group_by(Sample, Subtype, Genus_labeled) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  # 6. Pivot to wide matrix
  heat_mat <- heat_df %>%
    pivot_wider(
      names_from = Genus_labeled,
      values_from = Abundance,
      values_fill = 0
    ) %>%
    column_to_rownames("Sample") %>%
    select(-Subtype) %>%
    as.matrix()
  
  if(nrow(heat_mat) < 1 || ncol(heat_mat) < 1) stop("No data left to plot heatmap.")
  
  # 7. Log-transform safely
  heat_mat_log <- log10(heat_mat + 1e-6)
  heat_mat_log[!is.finite(heat_mat_log)] <- 0
  
  # 8. Scale columns and add tiny noise to fully constant columns
  heat_mat_scaled <- apply(heat_mat_log, 2, function(x) {
    if(sd(x, na.rm = TRUE) == 0) {
      x <- rep(0, length(x))
    } else {
      x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    }
    if(all(x == 0)) x <- x + rnorm(length(x), 0, 1e-6)
    x
  })
  heat_mat_scaled <- as.matrix(heat_mat_scaled)
  heat_mat_scaled[!is.finite(heat_mat_scaled)] <- 0
  
  # 9. Build annotation dataframe
  annotation_df_safe <- heat_df %>%
    select(Sample, Subtype) %>%
    distinct() %>%
    column_to_rownames("Sample")
  
  # 10. Safe clustering flags
  cluster_rows_flag <- nrow(heat_mat_scaled) >= 2
  cluster_cols_flag <- ncol(heat_mat_scaled) >= 2
  
  # 11. Assign colors for factors
  annotation_colors <- list()
  if(ncol(annotation_df_safe) > 0) {
    for(col_name in colnames(annotation_df_safe)) {
      if(is.factor(annotation_df_safe[[col_name]])) {
        levels_col <- levels(annotation_df_safe[[col_name]])
        annotation_colors[[col_name]] <- scales::hue_pal()(length(levels_col))
        names(annotation_colors[[col_name]]) <- levels_col
      }
    }
  } else {
    annotation_df_safe <- NULL
    annotation_colors <- NULL
  }
  
  # 12. Draw heatmap
  p_heat <- pheatmap(
    mat = heat_mat_scaled,
    annotation_row = annotation_df_safe,
    annotation_colors = annotation_colors,
    scale = "none",
    clustering_method = "complete",
    cluster_rows = cluster_rows_flag,
    cluster_cols = cluster_cols_flag,
    fontsize_row = 8,
    fontsize_col = 8,
    main = "Significant Genera Across Intestinal and Diffuse Samples (* q<0.05, ** q<0.01)"
  )
  
  # 13. Save
  ggsave("../figures/aim2_significant_genera_heatmap.png", plot = p_heat$gtable, width = 8, height = 8, dpi = 300)
  
}
# Volcano plot

# create a volcano plot data frame from ANCOM-BC2 results
volcano_df <- ancom_tbl %>%
  mutate(
    neg_log10_q = -log10(q),
    label = ifelse(significant == TRUE & q < 0.05, Genus, NA)
  )

# create the volcano plot
p_volcano <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_q, color = enriched_in)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "solid") +
  geom_text_repel(aes(label = label), na.rm = TRUE, size = 3, max.overlaps = 20) +
  labs(
    title = "ANCOM-BC2 differential abundance by subtype",
    x = "Log fold change (Diffuse vs Intestinal)",
    y = "-log10(q-value)",
    color = "Direction"
  ) +
  theme_bw()

# save the volcano plot
ggsave("../figures/aim2_volcano_plot.png", p_volcano, width = 8, height = 6, dpi = 300)

# Final checks

# print final subtype counts and key output summaries
cat("Final subtype counts:\n")
print(table(sample_data(ps_genus_filt)$Subtype))
cat("Final filtered genus-level taxa:", ntaxa(ps_genus_filt), "\n")
cat("Significant ANCOM-BC2 genera:", nrow(ancom_sig), "\n")
cat("Significant indicator genera:", nrow(ind_sig), "\n")


# update Barplot with Significance Asterisks
p_bar <- ggplot(plot_bar_top, aes(x = Subtype, y = mean_abundance, fill = Genus)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Top 15 genera by subtype", x = "Subtype", y = "Relative abundance") +
  theme_bw()


print(p_bar)
ggsave("../figures/aim2_top15_genera_barplot.png", p_bar, width = 8, height = 6, dpi = 300)
# update Volcano Plot with Explicit q-values
volcano_df <- ancom_tbl %>%
  mutate(
    neg_log10_q = -log10(q),
    label = ifelse(significant == TRUE & q < 0.05, paste0(Genus, "\n(q = ", signif(q, 2), ")"), NA)
  )

# no significant bacteria in top 15 most abundant 
intersect(
  top_genera,
  ancom_tbl %>% filter(q < 0.05) %>% pull(Genus)
)
p_volcano_sig <- ggplot(volcano_df, aes(x = logFC, y = neg_log10_q, color = enriched_in)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "solid") +
  geom_text_repel(aes(label = label), na.rm = TRUE, size = 3, max.overlaps = 20, lineheight = 0.8) +
  labs(
    title = "ANCOM-BC2 differential abundance by subtype",
    x = "Log fold change (Diffuse vs Intestinal)",
    y = "-log10(q-value)",
    color = "Direction"
  ) +
  theme_bw()

# display and save volcano plot interactively
print(p_volcano_sig)
ggsave("../figures/aim2_volcano_plot.png", p_volcano_sig, width = 8, height = 6, dpi = 300)

# Update Heatmap with Significance Asterisks
if (length(sig_genera) > 0) {
  
  # 1. Melt phyloseq object and filter significant genera
  heat_df_base <- psmelt(ps_genus_rel) %>% filter(OTU %in% sig_genera)
  
  # 2. Extract q-values from both ANCOM-BC2 and Indicator Species tests
  q_ancom <- ancom_tbl %>% select(OTU = Genus, q_ancom = q)
  q_ind   <- ind_tbl    %>% select(OTU = Genus, q_ind   = q)
  
  # 3. Combine q-values and determine lowest q-value for each taxa
  q_dict <- data.frame(OTU = sig_genera) %>%
    left_join(q_ancom, by = "OTU") %>%
    left_join(q_ind, by = "OTU") %>%
    rowwise() %>%
    mutate(
      best_q = min(c(q_ancom, q_ind), na.rm = TRUE),
      sig_stars = case_when(
        best_q < 0.001 ~ "***",
        best_q < 0.01  ~ "**",
        best_q < 0.05  ~ "*",
        TRUE           ~ ""
      )
    ) %>%
    ungroup() %>%
    select(OTU, sig_stars)
  
  # 4. Inject stars into genus names and aggregate duplicates
  heat_df <- heat_df_base %>%
    left_join(q_dict, by = "OTU") %>%
    mutate(
      sig_stars = ifelse(is.na(sig_stars), "", sig_stars),
      Genus_labeled = paste0(as.character(Genus), sig_stars)
    ) %>%
    group_by(Sample, Subtype, Genus_labeled) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = Genus_labeled,
      values_from = Abundance,
      values_fill = 0
    )
  
  # 5. Build heat matrix
  heat_mat <- heat_df %>%
    column_to_rownames("Sample") %>%
    select(-Subtype) %>%
    as.matrix()
  
  # 6. Log-transform safely
  heat_mat_log <- log10(heat_mat + 1e-6)
  
  # 7. Build row annotation
  annotation_df <- heat_df %>%
    select(Sample, Subtype) %>%
    distinct() %>%
    column_to_rownames("Sample")
  
  # 8. Draw heatmap
  p_heat <- pheatmap(
    mat = heat_mat_log,
    annotation_row = annotation_df,
    scale = "column",
    clustering_method = "complete",
    fontsize_row = 8,
    fontsize_col = 8,
    main = "Significant Genera Across Intestinal and Diffuse Samples (* q<0.05, ** q<0.01)"
  )
  
  # 9. Save heatmap
  ggsave("../figures/aim2_significant_genera_heatmap.png", plot = p_heat$gtable,
         width = 8, height = 8, dpi = 300)
}


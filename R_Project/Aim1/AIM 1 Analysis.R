# Loads in relevant packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)

# Creates variable for phyloseq object from Phyloseq.R script
ps <- readRDS("ps_raw.rds")

# Extracts metadata from phyloseq into dataframes
md <- as.data.frame(sample_data(ps))
md$SampleID <- rownames(md)

# Cleans up columns for Aim 1
grp <- trimws(md$Group)
lauren <- trimws(md[["Histopathology.......lauren.classification."]])

# Creates subtype variables for plotting and stats
md$Subtype <- NA
md$Subtype[lauren == "Intestinal type"] <- "Intestinal"
md$Subtype[lauren == "Diffused type"]   <- "Diffuse"

# Keeps only Gastric cancer subtypes that are intestinal or diffuse
# i.e. drops mixed and those without cancer
keep_a1 <- md$SampleID[grp == "Gastric cancer (GC)" & md$Subtype %in% c("Intestinal","Diffuse")]
cat("keep_a1:", length(keep_a1), "\n")

# Prunes and subsets phlyoseq object for Aim 1
# Drops any taxa that may become zero after filtering
ps_a1 <- prune_samples(keep_a1, ps)
ps_a1 <- prune_taxa(taxa_sums(ps_a1) > 0, ps_a1)

# Replaces sample_data with only the columns we'll need
# Group, for gastric cancers only
# Gender, for gender analyses
# Subtype, for comparing intestinal vs diffuse
md_a1 <- md[md$SampleID %in% sample_names(ps_a1), c("Group","Gender","Subtype")]
sample_data(ps_a1) <- sample_data(md_a1)

# Performs alpha rarefaction before diversity analyses
depth <- 8386  # sets to rarefaction depth lower than proposed initially (for consistnecy with aim 1b)

# If any samples are below the depth, drops them first
low_n <- sum(sample_sums(ps_a1) < depth)
if (low_n > 0) {
  message("Dropping ", low_n, " samples with < ", depth, " reads before rarefaction.")
  ps_a1 <- prune_samples(sample_sums(ps_a1) >= depth, ps_a1)
  ps_a1 <- prune_taxa(taxa_sums(ps_a1) > 0, ps_a1)
}

set.seed(1)
ps_a1 <- rarefy_even_depth(ps_a1, sample.size = depth, rngseed = 1,
                           replace = FALSE, trimOTUs = TRUE, verbose = FALSE)

# Verifies depth of alpha rarefaction, should be 8386
# Length of unique samples should be 1
s <- sample_sums(ps_a1)
summary(s)
length(unique(s))

# -----
# Aim 1A: alpha diversity
# -----
# We quantify diversity within each sample, then test whether
# Intestinal vs Diffuse GC differs in Shannon diversity after rarefaction

# Computes shannon for each sample
# Adds subtype labels for plotting and sta
alpha_a1 <- estimate_richness(ps_a1, measures = c("Shannon"))
alpha_a1$Subtype <- sample_data(ps_a1)$Subtype

# Violin plot with jittered overlay of Shannon by subtype
p_alpha <- ggplot(alpha_a1, aes(x = Subtype, y = Shannon)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.12, alpha = 0.7) +
  theme_bw() +
  labs(title = "Shannon Alpha Diversity by Lauren Subtype")

# Plot and saves PNG file
p_alpha
ggsave("aim1_alpha_shannon.png", p_alpha, width = 6, height = 4, dpi = 300)

# Wilcox rank-sum test (non-parametric)
wilcox.test(Shannon ~ Subtype, data = alpha_a1)

# -----
# Aim 1A: beta diversity
# -----
# We quantify how different microbial communities are between samples,
# visualize with PCoA plots, and then test if the subtypes explain
# community differences with PERMANOVA. Dispersion is also tested
# in case of confounding if groups have of within-group variance

# Extracts sample metadata with subtype as a categorical factor
meta_a1 <- data.frame(sample_data(ps_a1))
meta_a1$Subtype <- factor(meta_a1$Subtype)

run_beta <- function(ps_obj, dist_method, label) {
  
  # Pairwise comparison of sample distance
  d <- phyloseq::distance(ps_obj, method = dist_method)
  
  # PCoA ordintion for visualization
  ord <- ordinate(ps_obj, method = "PCoA", distance = d)
  p <- plot_ordination(ps_obj, ord, color = "Subtype") +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = paste(label, "PCoA"))
  
  # Saves plot for report
  ggsave(paste0("aim1_pcoa_", label, ".png"), p, width = 6, height = 4, dpi = 300)
  
  # Checks dispersion (within-group variance)
  # If there is significant dispersion, PERMANOVA can be misleading 
  bd <- betadisper(as.dist(d), group = meta_a1$Subtype)
  bd_p <- anova(bd)$`Pr(>F)`[1]
  
  # PERMANOVA with permutations to report p-value
  ad <- adonis2(as.dist(d) ~ Subtype, data = meta_a1, permutations = 9999)
  
  # Returns list
  list(dist = dist_method, betadisper_p = bd_p, permanova = ad)
}

# Run beta diversity for the three distance metrics used in Aim 1
res_bray <- run_beta(ps_a1, "bray",    "bray")
res_wuf  <- run_beta(ps_a1, "wunifrac","weighted_unifrac")
res_uuf  <- run_beta(ps_a1, "unifrac", "unweighted_unifrac")

# Print key results for reporting
# PERMANOVA, to check whether subtype explains community differences
# betadisper, to check whether within-group dispersion differs (potential PERMANOVA confound)
res_bray$permanova; res_bray$betadisper_p
res_wuf$permanova;  res_wuf$betadisper_p
res_uuf$permanova;  res_uuf$betadisper_p


# -----
# Aim 1B setup
# -----
# Aim 1B compares male vs female within each subtype separately.

ps_a1b <- ps_a1

# -----
# Aim 1B: alpha diversity within subtype
# -----
# Compare Shannon between male and female samples separately within:
#   1) Diffuse
#   2) Intestinal

run_a1b_alpha_gender_within_subtype <- function(ps_obj, subtype_label, out_prefix) {
  
  # Keep only one subtype at a time
  ps_sub <- prune_samples(sample_data(ps_obj)$Subtype == subtype_label, ps_obj)
  
  # Compute Shannon
  alpha_a1b <- estimate_richness(ps_sub, measures = c("Shannon"))
  alpha_a1b$Gender <- factor(trimws(sample_data(ps_sub)$Gender))
  
  # Plot Shannon by Gender within subtype
  p_alpha_a1b <- ggplot(alpha_a1b, aes(x = Gender, y = Shannon)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.12, alpha = 0.7) +
    theme_bw() +
    labs(title = paste("Shannon Alpha Diversity by Gender within", subtype_label),
         x = "Gender", y = "Shannon")
  
  ggsave(paste0(out_prefix, "_alpha_shannon_gender.png"),
         p_alpha_a1b, width = 6, height = 4, dpi = 300)
  
  # Wilcox rank-sum test for male vs female within subtype
  print(wilcox.test(Shannon ~ Gender, data = alpha_a1b, exact = FALSE))
}

# Run alpha comparisons within each subtype
run_a1b_alpha_gender_within_subtype(ps_a1b, "Diffuse",    "aim1b_diffuse")
run_a1b_alpha_gender_within_subtype(ps_a1b, "Intestinal", "aim1b_intestinal")

# -----
# Aim 1B: beta diversity within subtype
# -----
# Compare male vs female community composition separately within:
#   1) Diffuse
#   2) Intestinal

run_a1b_beta_gender_within_subtype <- function(ps_obj, subtype_label, dist_method, out_prefix) {
  
  # Keep only one subtype at a time
  ps_sub <- prune_samples(sample_data(ps_obj)$Subtype == subtype_label, ps_obj)
  
  # Metadata for stats
  meta_a1b <- data.frame(sample_data(ps_sub))
  meta_a1b$Gender <- factor(trimws(meta_a1b$Gender))
  
  # Pairwise distance matrix
  d_a1b <- phyloseq::distance(ps_sub, method = dist_method)
  
  # PCoA plot colored by Gender
  ord_a1b <- ordinate(ps_sub, method = "PCoA", distance = d_a1b)
  p_beta_a1b <- plot_ordination(ps_sub, ord_a1b, color = "Gender") +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = paste(subtype_label, "-", dist_method, "PCoA by Gender"))
  
  ggsave(paste0(out_prefix, "_pcoa_", dist_method, "_gender.png"),
         p_beta_a1b, width = 6, height = 4, dpi = 300)
  
  # Dispersion test for Gender within subtype
  bd_a1b <- betadisper(as.dist(d_a1b), group = meta_a1b$Gender)
  bd_p_a1b <- anova(bd_a1b)$`Pr(>F)`[1]
  
  # PERMANOVA for Gender within subtype
  ad_a1b <- adonis2(as.dist(d_a1b) ~ Gender, data = meta_a1b, permutations = 9999)
  
  list(dist = dist_method, betadisper_p = bd_p_a1b, permanova = ad_a1b)
}

# Diffuse: male vs female
a1b_diff_bray <- run_a1b_beta_gender_within_subtype(ps_a1b, "Diffuse", "bray",    "aim1b_diffuse")
a1b_diff_wuf  <- run_a1b_beta_gender_within_subtype(ps_a1b, "Diffuse", "wunifrac","aim1b_diffuse")
a1b_diff_uuf  <- run_a1b_beta_gender_within_subtype(ps_a1b, "Diffuse", "unifrac", "aim1b_diffuse")

# Intestinal: male vs female
a1b_int_bray <- run_a1b_beta_gender_within_subtype(ps_a1b, "Intestinal", "bray",    "aim1b_intestinal")
a1b_int_wuf  <- run_a1b_beta_gender_within_subtype(ps_a1b, "Intestinal", "wunifrac","aim1b_intestinal")
a1b_int_uuf  <- run_a1b_beta_gender_within_subtype(ps_a1b, "Intestinal", "unifrac", "aim1b_intestinal")

# Print results for Diffuse
a1b_diff_bray$permanova; a1b_diff_bray$betadisper_p
a1b_diff_wuf$permanova;  a1b_diff_wuf$betadisper_p
a1b_diff_uuf$permanova;  a1b_diff_uuf$betadisper_p

# Print results for Intestinal
a1b_int_bray$permanova; a1b_int_bray$betadisper_p
a1b_int_wuf$permanova;  a1b_int_wuf$betadisper_p
a1b_int_uuf$permanova;  a1b_int_uuf$betadisper_p

# -----
# Aim 1C setup
# -----
# Aim 1C looks to whether patterns will differ by Gender, 
# and if subtype effect depends on Gender (Subtype X Gender)

# Reuses our Aim 1 filtered and rareified phyloseq object
ps_a1c <- ps_a1

# Extracts sample metadata and ensures its categorical 

meta_a1c <- data.frame(sample_data(ps_a1c))
meta_a1c$Subtype <- factor(meta_a1c$Subtype)
meta_a1c$Gender  <- factor(trimws(meta_a1c$Gender))  #trims spaces if present

table(meta_a1c$Subtype, meta_a1c$Gender)

# -----
# Aim 1C: alpha diversity (interaction)
# -----
# Computes Shannon per sample and tests for:
# effects of subtype 
# effects of gender
# effect of subtype and gender interaction

alpha_a1c <- estimate_richness(ps_a1c, measures = c("Shannon"))
alpha_a1c$Subtype <- meta_a1c$Subtype
alpha_a1c$Gender  <- meta_a1c$Gender

# Plots Shannon diversity with trend lines and jittered points
p_alpha1c <- ggplot(alpha_a1c, aes(x = Subtype, y = Shannon, color = Gender, group = Gender)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.2)) +
  stat_summary(fun = mean, geom = "line") +
  geom_jitter(width = 0.12, alpha = 0.6) +
  theme_bw() +
  labs(title = "Shannon Alpha Diversity (Subtype x Gender)")

# Plot and save image
p_alpha1c
ggsave("aim1c_alpha_shannon_interaction.png", p_alpha1c, width = 6, height = 4, dpi = 300)

# Two-way ANOVA
fit_a1c_alpha <- aov(Shannon ~ Subtype * Gender, data = alpha_a1c)
summary(fit_a1c_alpha)


run_beta_a1c <- function(ps_obj, dist_method, label) {
  
  # Computes pairwise distance between samples
  d <- phyloseq::distance(ps_obj, method = dist_method)
  
  # PCoA ordination for visualization
  ord <- ordinate(ps_obj, method = "PCoA", distance = d)
  
  # Plots ordination, color by Subtype and shape by Gender
  p <- plot_ordination(ps_obj, ord, color = "Subtype", shape = "Gender") +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = paste(label, "PCoA (Subtype x Gender)"))
  
  # Saves plot
  ggsave(paste0("aim1c_pcoa_", label, ".png"), p, width = 6, height = 4, dpi = 300)
  
  # Extracts metadata for PERMANOVA and dispersion tests
  meta <- data.frame(sample_data(ps_obj))
  meta$Subtype <- factor(meta$Subtype)
  meta$Gender  <- factor(trimws(meta$Gender))
  
  # Dispersion check across subtype x gender groups
  combo <- interaction(meta$Subtype, meta$Gender, drop = TRUE)
  bd <- betadisper(as.dist(d), group = combo)
  bd_p <- anova(bd)$`Pr(>F)`[1]
  
  # PERMANOVA with term-by-term testing
  ad <- adonis2(
    as.dist(d) ~ Subtype * Gender,
    data = meta,
    permutations = 9999,
    by = "terms"
  )
  
  # Returns results
  list(dist = dist_method, betadisper_p = bd_p, permanova = ad)
}

# Runs Aim 1C beta diversity using same metrics as Aim 1A
a1c_bray <- run_beta_a1c(ps_a1c, "bray", "bray")
a1c_wuf  <- run_beta_a1c(ps_a1c, "wunifrac", "weighted_unifrac")
a1c_uuf  <- run_beta_a1c(ps_a1c, "unifrac", "unweighted_unifrac")

# Prints PERMANOVA tables
a1c_bray$permanova
a1c_wuf$permanova
a1c_uuf$permanova

# Prints dispersion p-values
a1c_bray$betadisper_p
a1c_wuf$betadisper_p
a1c_uuf$betadisper_p

# PERMANOVA shows significance, but so is dispersion
# May be driven by within group variance, check with TA
# weighted UniFrac gave a negative squared distances warning

# -----
# Aim 1C: EXTRA VISUALS
# -----
# Additional PCoAs faceted by Gender
# Allows easy visual comparison within each gender

# Bray PCoA faceted by Gender 
ord_bray_a1c <- ordinate(ps_a1c, method = "PCoA", distance = "bray")
p_pcoa_bray_facet <- plot_ordination(ps_a1c, ord_bray_a1c, color = "Subtype") +
  geom_point(size = 3) +
  facet_wrap(~ Gender) +
  theme_bw() +
  labs(title = "Bray PCoA by subtype, faceted by gender")
ggsave("aim1c_pcoa_bray_facet_gender.png", p_pcoa_bray_facet, width = 8, height = 4, dpi = 300)

# Unweighted UniFrac PCoA faceted by Gender 
ord_uuf_a1c <- ordinate(ps_a1c, method = "PCoA", distance = "unifrac")
p_pcoa_uuf_facet <- plot_ordination(ps_a1c, ord_uuf_a1c, color = "Subtype") +
  geom_point(size = 3) +
  facet_wrap(~ Gender) +
  theme_bw() +
  labs(title = "Unweighted UniFrac PCoA by subtype, faceted by gender")
ggsave("aim1c_pcoa_unweighted_unifrac_facet_gender.png", p_pcoa_uuf_facet, width = 8, height = 4, dpi = 300)

# -----
# Taxa barplot stratified by Gender (Genus level, top 10)
# -----
# Collapses features to the Genus level and converts counts to relative abundance per sample.
# Plots the top 10 genera overall and groups everything else into "Other".
# Faceted by Gender to compare patterns within Gender.

# Taxonomy aggregation and normalization
ps_genus <- tax_glom(ps_a1c, taxrank = "Genus", NArm = TRUE)
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

# Finds top 10 genus by abundance
topN <- 10
top_taxa <- names(sort(taxa_sums(ps_genus_rel), decreasing = TRUE))[1:topN]

# Grabs top genus, then groups everyting else as "others" 
# Makes it a 1 row OTU table to merge back with plot data
ps_top <- prune_taxa(top_taxa, ps_genus_rel)
ps_other <- prune_taxa(!taxa_names(ps_genus_rel) %in% top_taxa, ps_genus_rel)
other_vec <- sample_sums(ps_other)
other_mat <- matrix(other_vec, nrow = 1)
rownames(other_mat) <- "Other"
colnames(other_mat) <- sample_names(ps_genus_rel)

# Creates minimal taxonomy row for "Other" so phyloseq can treat it as a taxon
OTU_other <- otu_table(other_mat, taxa_are_rows = TRUE)
TAX_other <- tax_table(matrix(
  c("Other", NA, NA, NA, NA, "Other", NA),
  nrow = 1,
  dimnames = list("Other", colnames(tax_table(ps_genus_rel)))
))

# Merges samples and creates bar plot for abundance stratified by gender
ps_bar <- merge_phyloseq(ps_top, OTU_other, TAX_other, sample_data(ps_genus_rel))
p_bar <- plot_bar(ps_bar, x = "Subtype", fill = "Genus") +
  facet_wrap(~ Gender) +
  theme_bw() +
  labs(title = "Top genera relative abundance by subtype, stratified by gender")

# Saves image for report 
ggsave("aim1c_taxa_barplot_top10_genus_facet_gender.png", p_bar, width = 10, height = 5, dpi = 300)







# -----
# FINAL FIGURE CLEANUP
# -----


library(scales)

# -----------------------------
# Shared styling and helpers
# -----------------------------

fig_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text   = element_text(face = "bold")
  )

subtype_cols <- c(
  "Diffuse"    = "#D55E00",
  "Intestinal" = "#0072B2"
)

gender_cols <- c(
  "Female" = "#CC79A7",
  "Male"   = "#009E73"
)

gender_shapes <- c(
  "Female" = 16,
  "Male"   = 17
)

fmt_p <- function(p) {
  format.pval(p, digits = 3, eps = 1e-04)
}

get_axis_pct <- function(ord, axis = 1) {
  vals <- ord$values$Relative_eig
  if (is.null(vals)) {
    vals <- ord$values$Eigenvalues / sum(ord$values$Eigenvalues)
  }
  round(as.numeric(vals[axis]) * 100, 1)
}

clean_gender <- function(x) {
  x <- trimws(as.character(x))
  x <- ifelse(tolower(x) == "female", "Female",
              ifelse(tolower(x) == "male", "Male", x))
  factor(x, levels = c("Female", "Male"))
}

# -----------------------------
# Helper: alpha plot with violin
# -----------------------------
make_alpha_violin <- function(df, xvar, title_text, subtitle_text,
                              xlab_text, fill_vals, file_name,
                              caption_text = NULL) {
  
  p <- ggplot(df, aes(x = .data[[xvar]], y = Shannon, fill = .data[[xvar]])) +
    geom_violin(trim = FALSE, alpha = 0.25, color = "black") +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6, color = "black") +
    geom_jitter(width = 0.08, size = 2.6, alpha = 0.8, color = "black") +
    scale_fill_manual(values = fill_vals) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      caption = caption_text,
      x = xlab_text,
      y = "Shannon Diversity"
    ) +
    fig_theme +
    theme(legend.position = "none")
  
  ggsave(file_name, p, width = 7, height = 5, dpi = 300)
  p
}

# -----------------------------
# Helper: alpha plot for small n
# -----------------------------
make_alpha_box_small <- function(df, xvar, title_text, subtitle_text,
                                 xlab_text, fill_vals, file_name,
                                 caption_text = NULL) {
  
  p <- ggplot(df, aes(x = .data[[xvar]], y = Shannon, fill = .data[[xvar]])) +
    geom_boxplot(width = 0.35, outlier.shape = NA, alpha = 0.5, color = "black") +
    geom_jitter(width = 0.08, size = 3, alpha = 0.8, color = "black") +
    scale_fill_manual(values = fill_vals) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      caption = caption_text,
      x = xlab_text,
      y = "Shannon Diversity"
    ) +
    fig_theme +
    theme(legend.position = "none")
  
  ggsave(file_name, p, width = 7, height = 5, dpi = 300)
  p
}

# -----------------------------
# Helper: subtype-only PCoA
# -----------------------------
make_clean_subtype_pcoa <- function(ps_obj, dist_method, title_text, file_name, stat_res) {
  
  d <- phyloseq::distance(ps_obj, method = dist_method)
  ord <- ordinate(ps_obj, method = "PCoA", distance = d)
  
  plot_df <- plot_ordination(ps_obj, ord, justDF = TRUE)
  plot_df$Subtype <- factor(plot_df$Subtype, levels = c("Diffuse", "Intestinal"))
  
  x_pct <- get_axis_pct(ord, 1)
  y_pct <- get_axis_pct(ord, 2)
  
  perm_p <- fmt_p(stat_res$permanova$`Pr(>F)`[1])
  r2     <- round(stat_res$permanova$R2[1], 3)
  disp_p <- fmt_p(stat_res$betadisper_p)
  
  p <- ggplot(plot_df, aes(x = Axis.1, y = Axis.2)) +
    stat_ellipse(aes(color = Subtype, group = Subtype),
                 linewidth = 0.8, linetype = 2, type = "norm",
                 show.legend = FALSE) +
    geom_point(aes(fill = Subtype),
               size = 4, shape = 21, color = "black", stroke = 0.4, alpha = 0.9) +
    scale_fill_manual(values = subtype_cols) +
    scale_color_manual(values = subtype_cols) +
    labs(
      title = title_text,
      subtitle = paste0("PERMANOVA: p = ", perm_p, ", R2 = ", r2),
      caption = paste0("Betadisper: p = ", disp_p),
      x = paste0("PCoA1 (", x_pct, "%)"),
      y = paste0("PCoA2 (", y_pct, "%)"),
      fill = "Lauren Subtype"
    ) +
    fig_theme
  
  ggsave(file_name, p, width = 7.2, height = 5.2, dpi = 300)
  p
}

# -----------------------------
# Helper: gender-only PCoA
# -----------------------------
make_clean_gender_pcoa <- function(ps_obj, dist_method, title_text, file_name,
                                   stat_res, caption_extra = NULL) {
  
  d <- phyloseq::distance(ps_obj, method = dist_method)
  ord <- ordinate(ps_obj, method = "PCoA", distance = d)
  
  plot_df <- plot_ordination(ps_obj, ord, justDF = TRUE)
  plot_df$Gender <- clean_gender(plot_df$Gender)
  
  x_pct <- get_axis_pct(ord, 1)
  y_pct <- get_axis_pct(ord, 2)
  
  perm_p <- fmt_p(stat_res$permanova$`Pr(>F)`[1])
  r2     <- round(stat_res$permanova$R2[1], 3)
  disp_p <- fmt_p(stat_res$betadisper_p)
  
  cap_text <- paste0("Betadisper: p = ", disp_p)
  if (!is.null(caption_extra)) {
    cap_text <- paste0(cap_text, "; ", caption_extra)
  }
  
  p <- ggplot(plot_df, aes(x = Axis.1, y = Axis.2, fill = Gender)) +
    geom_point(size = 4, shape = 21, color = "black", stroke = 0.4, alpha = 0.9) +
    scale_fill_manual(values = gender_cols) +
    labs(
      title = title_text,
      subtitle = paste0("PERMANOVA: p = ", perm_p, ", R2 = ", r2),
      caption = cap_text,
      x = paste0("PCoA1 (", x_pct, "%)"),
      y = paste0("PCoA2 (", y_pct, "%)"),
      fill = "Gender"
    ) +
    fig_theme
  
  ggsave(file_name, p, width = 7.2, height = 5.2, dpi = 300)
  p
}

# -----------------------------
# Helper: interaction PCoA
# -----------------------------
make_clean_interaction_pcoa <- function(ps_obj, dist_method, title_text, file_name, stat_res) {
  
  d <- phyloseq::distance(ps_obj, method = dist_method)
  ord <- ordinate(ps_obj, method = "PCoA", distance = d)
  
  plot_df <- plot_ordination(ps_obj, ord, justDF = TRUE)
  plot_df$Subtype <- factor(plot_df$Subtype, levels = c("Diffuse", "Intestinal"))
  plot_df$Gender  <- clean_gender(plot_df$Gender)
  
  x_pct <- get_axis_pct(ord, 1)
  y_pct <- get_axis_pct(ord, 2)
  
  perm_p <- fmt_p(stat_res$permanova$`Pr(>F)`[1])
  r2     <- round(stat_res$permanova$R2[1], 3)
  disp_p <- fmt_p(stat_res$betadisper_p)
  
  p <- ggplot(plot_df, aes(x = Axis.1, y = Axis.2, color = Subtype, shape = Gender)) +
    geom_point(size = 4, alpha = 0.9) +
    scale_color_manual(values = subtype_cols) +
    scale_shape_manual(values = gender_shapes) +
    labs(
      title = title_text,
      subtitle = paste0("Overall PERMANOVA Model: p = ", perm_p, ", R2 = ", r2),
      caption = paste0("Betadisper: p = ", disp_p),
      x = paste0("PCoA1 (", x_pct, "%)"),
      y = paste0("PCoA2 (", y_pct, "%)"),
      color = "Lauren Subtype",
      shape = "Gender"
    ) +
    fig_theme
  
  ggsave(file_name, p, width = 7.6, height = 5.6, dpi = 300)
  p
}

# -----------------------------
# Helper: faceted PCoA
# -----------------------------
make_clean_facet_pcoa <- function(ps_obj, dist_method, title_text, file_name,
                                  stat_res, caption_extra = NULL) {
  
  d <- phyloseq::distance(ps_obj, method = dist_method)
  ord <- ordinate(ps_obj, method = "PCoA", distance = d)
  
  plot_df <- plot_ordination(ps_obj, ord, justDF = TRUE)
  plot_df$Subtype <- factor(plot_df$Subtype, levels = c("Diffuse", "Intestinal"))
  plot_df$Gender  <- clean_gender(plot_df$Gender)
  
  x_pct <- get_axis_pct(ord, 1)
  y_pct <- get_axis_pct(ord, 2)
  
  perm_p <- fmt_p(stat_res$permanova$`Pr(>F)`[1])
  r2     <- round(stat_res$permanova$R2[1], 3)
  disp_p <- fmt_p(stat_res$betadisper_p)
  
  cap_text <- paste0("Betadisper: p = ", disp_p)
  if (!is.null(caption_extra)) {
    cap_text <- paste0(cap_text, "; ", caption_extra)
  }
  
  p <- ggplot(plot_df, aes(x = Axis.1, y = Axis.2, fill = Subtype)) +
    geom_point(size = 4, shape = 21, color = "black", stroke = 0.4, alpha = 0.9) +
    facet_wrap(~ Gender) +
    scale_fill_manual(values = subtype_cols) +
    labs(
      title = title_text,
      subtitle = paste0("Overall PERMANOVA Model: p = ", perm_p, ", R2 = ", r2),
      caption = cap_text,
      x = paste0("PCoA1 (", x_pct, "%)"),
      y = paste0("PCoA2 (", y_pct, "%)"),
      fill = "Lauren Subtype"
    ) +
    fig_theme
  
  ggsave(file_name, p, width = 9.2, height = 5.6, dpi = 300)
  p
}

# =========================================================
# AIM 1A
# =========================================================

# -----------------------------
# Aim 1A: cleaned alpha by subtype
# -----------------------------
alpha_a1_clean <- estimate_richness(ps_a1, measures = "Shannon")
alpha_a1_clean$Subtype <- factor(data.frame(sample_data(ps_a1))$Subtype,
                                 levels = c("Diffuse", "Intestinal"))

p_a1_alpha_clean <- make_alpha_violin(
  df = alpha_a1_clean,
  xvar = "Subtype",
  title_text = "Shannon Diversity by Lauren Subtype",
  subtitle_text = paste0(
    "Wilcoxon Rank-Sum: p = ",
    fmt_p(wilcox.test(Shannon ~ Subtype, data = alpha_a1_clean, exact = FALSE)$p.value)
  ),
  xlab_text = "Lauren Subtype",
  fill_vals = subtype_cols,
  file_name = "clean_aim1_alpha_shannon_subtype.png"
)

# -----------------------------
# Aim 1A: cleaned beta PCoAs by subtype
# -----------------------------
p_a1_bray_clean <- make_clean_subtype_pcoa(
  ps_a1, "bray",
  "Bray-Curtis PCoA",
  "clean_aim1_pcoa_bray_subtype.png",
  res_bray
)

p_a1_wuf_clean <- make_clean_subtype_pcoa(
  ps_a1, "wunifrac",
  "Weighted UniFrac PCoA",
  "clean_aim1_pcoa_weighted_unifrac_subtype.png",
  res_wuf
)

p_a1_uuf_clean <- make_clean_subtype_pcoa(
  ps_a1, "unifrac",
  "Unweighted UniFrac PCoA",
  "clean_aim1_pcoa_unweighted_unifrac_subtype.png",
  res_uuf
)

# =========================================================
# AIM 1B - DIFFUSE
# =========================================================

# -----------------------------
# Aim 1B Diffuse: cleaned alpha by gender
# -----------------------------
ps_diffuse_clean <- prune_samples(sample_data(ps_a1b)$Subtype == "Diffuse", ps_a1b)

alpha_diff_clean <- estimate_richness(ps_diffuse_clean, measures = "Shannon")
alpha_diff_clean$Gender <- clean_gender(data.frame(sample_data(ps_diffuse_clean))$Gender)

p_diff_alpha_clean <- make_alpha_violin(
  df = alpha_diff_clean,
  xvar = "Gender",
  title_text = "Shannon Diversity by Gender Within Diffuse Subtype",
  subtitle_text = paste0(
    "Wilcoxon Rank-Sum: p = ",
    fmt_p(wilcox.test(Shannon ~ Gender, data = alpha_diff_clean, exact = FALSE)$p.value)
  ),
  xlab_text = "Gender",
  fill_vals = gender_cols,
  file_name = "clean_aim1b_diffuse_alpha_shannon_gender.png"
)

# -----------------------------
# Aim 1B Diffuse: cleaned beta PCoAs by gender
# -----------------------------
p_diff_bray_clean <- make_clean_gender_pcoa(
  ps_diffuse_clean, "bray",
  "Diffuse Subtype: Bray-Curtis PCoA by Gender",
  "clean_aim1b_diffuse_pcoa_bray_gender.png",
  a1b_diff_bray
)

p_diff_wuf_clean <- make_clean_gender_pcoa(
  ps_diffuse_clean, "wunifrac",
  "Diffuse Subtype: Weighted UniFrac PCoA by Gender",
  "clean_aim1b_diffuse_pcoa_weighted_unifrac_gender.png",
  a1b_diff_wuf
)

p_diff_uuf_clean <- make_clean_gender_pcoa(
  ps_diffuse_clean, "unifrac",
  "Diffuse Subtype: Unweighted UniFrac PCoA by Gender",
  "clean_aim1b_diffuse_pcoa_unweighted_unifrac_gender.png",
  a1b_diff_uuf
)

# =========================================================
# AIM 1B - INTESTINAL
# =========================================================

# -----------------------------
# Aim 1B Intestinal: cleaned alpha by gender
# -----------------------------
ps_intestinal_clean <- prune_samples(sample_data(ps_a1b)$Subtype == "Intestinal", ps_a1b)

alpha_int_clean <- estimate_richness(ps_intestinal_clean, measures = "Shannon")
alpha_int_clean$Gender <- clean_gender(data.frame(sample_data(ps_intestinal_clean))$Gender)

n_female_int <- sum(alpha_int_clean$Gender == "Female")
n_male_int   <- sum(alpha_int_clean$Gender == "Male")

p_int_alpha_clean <- make_alpha_box_small(
  df = alpha_int_clean,
  xvar = "Gender",
  title_text = "Shannon Diversity by Gender Within Intestinal Subtype",
  subtitle_text = paste0(
    "Wilcoxon Rank-Sum: p = ",
    fmt_p(wilcox.test(Shannon ~ Gender, data = alpha_int_clean, exact = FALSE)$p.value)
  ),
  xlab_text = "Gender",
  fill_vals = gender_cols,
  file_name = "clean_aim1b_intestinal_alpha_shannon_gender.png",
  caption_text = paste0("Female n = ", n_female_int, ", Male n = ", n_male_int)
)

# -----------------------------
# Aim 1B Intestinal: cleaned beta PCoAs by gender
# -----------------------------
int_caption <- paste0("Female n = ", n_female_int, ", Male n = ", n_male_int)

p_int_bray_clean <- make_clean_gender_pcoa(
  ps_intestinal_clean, "bray",
  "Intestinal Subtype: Bray-Curtis PCoA by Gender",
  "clean_aim1b_intestinal_pcoa_bray_gender.png",
  a1b_int_bray,
  caption_extra = int_caption
)

p_int_wuf_clean <- make_clean_gender_pcoa(
  ps_intestinal_clean, "wunifrac",
  "Intestinal Subtype: Weighted UniFrac PCoA by Gender",
  "clean_aim1b_intestinal_pcoa_weighted_unifrac_gender.png",
  a1b_int_wuf,
  caption_extra = int_caption
)

p_int_uuf_clean <- make_clean_gender_pcoa(
  ps_intestinal_clean, "unifrac",
  "Intestinal Subtype: Unweighted UniFrac PCoA by Gender",
  "clean_aim1b_intestinal_pcoa_unweighted_unifrac_gender.png",
  a1b_int_uuf,
  caption_extra = int_caption
)

# =========================================================
# AIM 1C
# =========================================================

# -----------------------------
# Aim 1C: cleaned alpha interaction plot
# -----------------------------
alpha_a1c_clean <- estimate_richness(ps_a1c, measures = "Shannon")
meta_a1c_clean <- data.frame(sample_data(ps_a1c))
meta_a1c_clean$Subtype <- factor(meta_a1c_clean$Subtype, levels = c("Diffuse", "Intestinal"))
meta_a1c_clean$Gender  <- clean_gender(meta_a1c_clean$Gender)

alpha_a1c_clean$Subtype <- meta_a1c_clean$Subtype
alpha_a1c_clean$Gender  <- meta_a1c_clean$Gender

fit_a1c_alpha_clean <- aov(Shannon ~ Subtype * Gender, data = alpha_a1c_clean)
aov_tab <- summary(fit_a1c_alpha_clean)[[1]]

p_sub  <- fmt_p(aov_tab["Subtype", "Pr(>F)"])
p_gen  <- fmt_p(aov_tab["Gender", "Pr(>F)"])
p_int  <- fmt_p(aov_tab["Subtype:Gender", "Pr(>F)"])

alpha_sum_a1c <- alpha_a1c_clean %>%
  group_by(Subtype, Gender) %>%
  summarise(
    mean_shannon = mean(Shannon),
    se_shannon   = sd(Shannon) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p_a1c_alpha_clean <- ggplot() +
  geom_point(
    data = alpha_a1c_clean,
    aes(x = Subtype, y = Shannon, color = Gender),
    position = position_jitter(width = 0.08, height = 0),
    alpha = 0.45,
    size = 2.7
  ) +
  geom_line(
    data = alpha_sum_a1c,
    aes(x = Subtype, y = mean_shannon, color = Gender, group = Gender),
    linewidth = 1
  ) +
  geom_point(
    data = alpha_sum_a1c,
    aes(x = Subtype, y = mean_shannon, color = Gender),
    size = 3.7
  ) +
  geom_errorbar(
    data = alpha_sum_a1c,
    aes(x = Subtype,
        ymin = mean_shannon - se_shannon,
        ymax = mean_shannon + se_shannon,
        color = Gender),
    width = 0.08,
    linewidth = 0.8
  ) +
  scale_color_manual(values = gender_cols) +
  labs(
    title = "Shannon Diversity by Subtype and Gender",
    subtitle = paste0(
      "Two-Way ANOVA: subtype p = ", p_sub,
      ", gender p = ", p_gen,
      ", interaction p = ", p_int
    ),
    caption = "Interpret cautiously: Intestinal Female n = 3",
    x = "Lauren Subtype",
    y = "Shannon Diversity",
    color = "Gender"
  ) +
  fig_theme

ggsave("clean_aim1c_alpha_shannon_interaction.png",
       p_a1c_alpha_clean, width = 7.8, height = 5.6, dpi = 300)

# -----------------------------
# Aim 1C: cleaned interaction PCoAs (Subtype x Gender)
# -----------------------------
p_a1c_bray_clean <- make_clean_interaction_pcoa(
  ps_a1c, "bray",
  "Bray-Curtis PCoA (Subtype x Gender)",
  "clean_aim1c_pcoa_bray_interaction.png",
  a1c_bray
)

p_a1c_wuf_clean <- make_clean_interaction_pcoa(
  ps_a1c, "wunifrac",
  "Weighted UniFrac PCoA (Subtype x Gender)",
  "clean_aim1c_pcoa_weighted_unifrac_interaction.png",
  a1c_wuf
)

p_a1c_uuf_clean <- make_clean_interaction_pcoa(
  ps_a1c, "unifrac",
  "Unweighted UniFrac PCoA (Subtype x Gender)",
  "clean_aim1c_pcoa_unweighted_unifrac_interaction.png",
  a1c_uuf
)

# -----------------------------
# Aim 1C: cleaned faceted PCoAs by gender
# -----------------------------
p_a1c_bray_facet_clean <- make_clean_facet_pcoa(
  ps_a1c, "bray",
  "Bray-Curtis PCoA by Subtype, Faceted by Gender",
  "clean_aim1c_pcoa_bray_facet_gender.png",
  a1c_bray,
  caption_extra = "Intestinal Female n = 3"
)

p_a1c_uuf_facet_clean <- make_clean_facet_pcoa(
  ps_a1c, "unifrac",
  "Unweighted UniFrac PCoA by Subtype, Faceted by Gender",
  "clean_aim1c_pcoa_unweighted_unifrac_facet_gender.png",
  a1c_uuf,
  caption_extra = "Intestinal Female n = 3"
)

# -----------------------------
# Aim 1C: cleaned top genera barplot
# -----------------------------
ps_genus_clean <- tax_glom(ps_a1c, taxrank = "Genus", NArm = TRUE)
ps_genus_rel_clean <- transform_sample_counts(ps_genus_clean, function(x) x / sum(x))

bar_df <- psmelt(ps_genus_rel_clean)
bar_df$Genus <- as.character(bar_df$Genus)
bar_df$Genus[is.na(bar_df$Genus)] <- "Unclassified"
bar_df$Gender <- clean_gender(bar_df$Gender)
bar_df$Subtype <- factor(bar_df$Subtype, levels = c("Diffuse", "Intestinal"))

top_genera <- bar_df %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 8) %>%
  pull(Genus)

bar_df_clean <- bar_df %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other")) %>%
  group_by(Gender, Subtype, Genus) %>%
  summarise(MeanAbundance = mean(Abundance), .groups = "drop")

genus_order <- bar_df_clean %>%
  group_by(Genus) %>%
  summarise(total_abund = sum(MeanAbundance), .groups = "drop") %>%
  arrange(desc(total_abund)) %>%
  pull(Genus)

genus_order <- c(setdiff(genus_order, "Other"), "Other")
bar_df_clean$Genus <- factor(bar_df_clean$Genus, levels = genus_order)

p_a1c_bar_clean <- ggplot(bar_df_clean, aes(x = Subtype, y = MeanAbundance, fill = Genus)) +
  geom_col(color = "black", width = 0.7) +
  facet_wrap(~ Gender) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Mean Genus-Level Relative Abundance by Subtype and Gender",
    subtitle = "Top 8 Genera Shown; Remaining Genera Grouped as Other",
    x = "Lauren Subtype",
    y = "Mean Relative Abundance",
    fill = "Genus"
  ) +
  fig_theme

ggsave("clean_aim1c_taxa_barplot_top_genera.png",
       p_a1c_bar_clean, width = 10.2, height = 5.8, dpi = 300)
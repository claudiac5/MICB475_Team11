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

# Adds a nonrareified version of phyloseq object for Aim 3 later
# See Aim 3 comments
ps_a1_unrarefied <- ps_a1

# Performs alpha rarefaction before diversity analyses
depth <- 10222  # set to rarefaction depth as outlined in proposal (should lose 2 samples)

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

# Verifies depth of alpha rarefaction, should be 10222
# Length of unique samples should be 1
s <- sample_sums(ps_a1)
summary(s)
length(unique(s))

# -----
# Aim 1: alpha diversity
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

# T test for significance (parametric assumption)
t.test(Shannon ~ Subtype, data = alpha_a1) 

# signfiicant, P = 0.01315

# Wilcox rank-sum test (non-parametric)
# wilcox.test(Shannon ~ Subtype, data = alpha_a1)

# Actually shows huge difference in significance?
# P = 0.2219 

# -----
# Aim 1: beta diversity
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

# PERMANOVA shows significance for Bray-Curtis, significant (p = 0.0025) and explains modest variance
# Betadisper (p = 0.795), dispersion isn't very different so PERMANOVA isn't driven by unequal spread

# Both weighted and unweighted unifrac are not significant
# Weighted unifrac shows no subtype signal by abundance
# Unweighted is close to signifiance, (p = 0.052), weak trend in presence/absence but doesn't cross 0.05


# -----
# Aim 3 setup
# -----
# Aim 3 looks to whether patterns will differ by Gender, 
# and if subtype effect depends on Gender (Subtype X Gender)

# Reuses our Aim 1 filtered and rareified phyloseq object
ps_a3 <- ps_a1

# Extracts sample metadata and ensures its categorical 
meta_a3 <- data.frame(sample_data(ps_a3))
meta_a3$Subtype <- factor(meta_a3$Subtype)
meta_a3$Gender  <- factor(trimws(meta_a3$Gender))  #trims spaces if present

table(meta_a3$Subtype, meta_a3$Gender)
# ONLY 1 FEMALE INTESTINAL LEFT? RAREFACTION EXCLUSIVELY DROPS FEMALE INTESTINAL
# SEEK HELP
# STATS POWER IS ENOUGH?

# Performs alpha rarefaction but with lower depth for aim 3,
# to maintain all 3 female intestinal samples
depth <- 8386  # set to lower depth now

ps_a3 <- ps_a1_unrarefied

# If any samples are below the depth, drops them first
low_n <- sum(sample_sums(ps_a3) < depth)
if (low_n > 0) {
  message("Dropping ", low_n, " samples with < ", depth, " reads before rarefaction.")
  ps_a3 <- prune_samples(sample_sums(ps_a3) >= depth, ps_a3)
  ps_a3 <- prune_taxa(taxa_sums(ps_a3) > 0, ps_a3)
}

set.seed(1)
ps_a3 <- rarefy_even_depth(ps_a3, sample.size = depth, rngseed = 1,
                           replace = FALSE, trimOTUs = TRUE, verbose = FALSE)

ps_a3

# Repeats code just with lower rarefaction depth
meta_a3 <- data.frame(sample_data(ps_a3))
meta_a3$Subtype <- factor(meta_a3$Subtype)
meta_a3$Gender  <- factor(trimws(meta_a3$Gender))  #trims spaces if present

table(meta_a3$Subtype, meta_a3$Gender)


# -----
# Aim 3: alpha diversity (interaction)
# -----
# Computes Shannon per sample and tests for:
# effects of subtype 
# effects of gender
# effect of subtype and gender interaction

alpha_a3 <- estimate_richness(ps_a3, measures = c("Shannon"))
alpha_a3$Subtype <- meta_a3$Subtype
alpha_a3$Gender  <- meta_a3$Gender

# Plots Shannon diversity with trend lines and splits by gender 
p_alpha3 <- ggplot(alpha_a3, aes(x = Subtype, y = Shannon, color = Gender, group = Gender)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.2)) +
  stat_summary(fun = mean, geom = "line") +
  geom_jitter(width = 0.12, alpha = 0.6) +
  theme_bw() +
  labs(title = "Shannon Alpha Diversity (Subtype x Gender)")

# Plot and saves image
p_alpha3
ggsave("aim3_alpha_shannon_interaction.png", p_alpha3, width = 6, height = 4, dpi = 300)

fit_a3_alpha <- aov(Shannon ~ Subtype * Gender, data = alpha_a3)
summary(fit_a3_alpha)
# Two-way ANOVA shows significant effect of subtype and gender on Shannon
# Careful with sample size of 3 for female intestinal though

# -----
# Aim 3: beta diversity (Subtype * Gender)
# -----
# Beta diversity comparison now for Gender and Subtype interactions
# We also check dispersion because PERMANOVA again

run_beta_a3 <- function(ps_obj, dist_method, label) {
  
  # Computes pairwise distance between samples
  d <- phyloseq::distance(ps_obj, method = dist_method)
  
  # PCoA ordination for visualization
  ord <- ordinate(ps_obj, method = "PCoA", distance = d)
  
  # Plots ordination, color with Subtype, and shape with Gender
  p <- plot_ordination(ps_obj, ord, color = "Subtype", shape = "Gender") +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = paste(label, "PCoA (Subtype x Gender)"))
  
  # Saves plot for report
  ggsave(paste0("aim3_pcoa_", label, ".png"), p, width = 6, height = 4, dpi = 300)
  
  # Extracts metadata for PERMANOVA and dispersion tests
  meta <- data.frame(sample_data(ps_obj))
  meta$Subtype <- factor(meta$Subtype)
  meta$Gender  <- factor(trimws(meta$Gender))
  
  # Dispersion check across subtype/gender
  combo <- interaction(meta$Subtype, meta$Gender, drop = TRUE)
  bd <- betadisper(as.dist(d), group = combo)
  bd_p <- anova(bd)$`Pr(>F)`[1]
  
  # PERMANOVA with permutations to report p-value 
  ad <- adonis2(as.dist(d) ~ Subtype * Gender, data = meta, permutations = 9999)
  
  # Returns results for reporting
  list(dist = dist_method, betadisper_p = bd_p, permanova = ad)
  
  
}

# Runs Aim 3 beta diversity using same metrics as Aim 1
a3_bray <- run_beta_a3(ps_a3, "bray", "bray")
a3_wuf  <- run_beta_a3(ps_a3, "wunifrac", "weighted_unifrac")
a3_uuf  <- run_beta_a3(ps_a3, "unifrac", "unweighted_unifrac")

# Prints PERMANOVA outputs
a3_bray$permanova
a3_wuf$permanova
a3_uuf$permanova

# Prints dispersion P values
a3_bray$betadisper_p
a3_wuf$betadisper_p
a3_uuf$betadisper_p

# PERMANOVA shows significance, but so is dispersion
# May be driven by within group variance, check with TA
# Bray: p ≈ 0.001
# weighted UniFrac: p ≈ 0.017 (also gave a negative squared distances warning)
# unweighted UniFrac: p ≈ 0.033

# -----
# Aim 3: EXTRA VISUALS
# -----
# Additional PCoAs faceted by Gender
# Allows easy visual comparison within each gender

# Bray PCoA faceted by Gender 
ord_bray_a3 <- ordinate(ps_a3, method = "PCoA", distance = "bray")
p_pcoa_bray_facet <- plot_ordination(ps_a3, ord_bray_a3, color = "Subtype") +
  geom_point(size = 3) +
  facet_wrap(~ Gender) +
  theme_bw() +
  labs(title = "Bray PCoA by subtype, faceted by gender")
ggsave("aim3_pcoa_bray_facet_gender.png", p_pcoa_bray_facet, width = 8, height = 4, dpi = 300)

# Unweighted UniFrac PCoA faceted by Gender 
ord_uuf_a3 <- ordinate(ps_a3, method = "PCoA", distance = "unifrac")
p_pcoa_uuf_facet <- plot_ordination(ps_a3, ord_uuf_a3, color = "Subtype") +
  geom_point(size = 3) +
  facet_wrap(~ Gender) +
  theme_bw() +
  labs(title = "Unweighted UniFrac PCoA by subtype, faceted by gender")
ggsave("aim3_pcoa_unweighted_unifrac_facet_gender.png", p_pcoa_uuf_facet, width = 8, height = 4, dpi = 300)

# -----
# Taxa barplot stratified by Gender (Genus level, top 10)
# -----
# Collapses features to the Genus level and converts counts to relative abundance per sample.
# Plots the top 10 genera overall and groups everything else into "Other".
# Faceted by Gender to compare patterns within Gender.

# Taxonomy aggregation and normalization
ps_genus <- tax_glom(ps_a3, taxrank = "Genus", NArm = TRUE)
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
ggsave("aim3_taxa_barplot_top10_genus_facet_gender.png", p_bar, width = 10, height = 5, dpi = 300)

# -----
# Taxa heatmap stratified by Gender
# -----

# Heatmap for genus level abundance alos stratified by gender
p_heat <- plot_heatmap(ps_genus_rel, taxa.label = "Genus", sample.label = "Subtype") +
  facet_wrap(~ Gender) +
  theme_bw() +
  labs(title = "Genus-level heatmap (relative abundance), stratified by gender")

# Saves image for report
ggsave("aim3_taxa_heatmap_genus_facet_gender.png", p_heat, width = 10, height = 6, dpi = 300)

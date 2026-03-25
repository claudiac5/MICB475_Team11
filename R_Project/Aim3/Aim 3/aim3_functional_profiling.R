# Loading Libraries
library(phyloseq)
library(tidyverse)
library(ANCOMBC)
library(pheatmap)

# Loading Files
pathway_data <- read.table("path_abun_unstrat.tsv.gz", 
                           header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

metadata <- read.delim("metadata.tsv", row.names = 1) %>%
  rename(Lauren_classification = `Histopathology.......lauren.classification.`)

# Filtering for our 45 GC samples
metadata_filtered <- metadata %>%
  filter(Group == "Gastric cancer (GC)") %>%
  filter(Lauren_classification %in% c("Intestinal type", "Diffused type"))

# Matching the columns from metadata and the pathway data
common_samples <- intersect(colnames(pathway_data), rownames(metadata_filtered))

pathway_final <- pathway_data[, common_samples]

metadata_final <- metadata_filtered[common_samples, ]

# Creating the phyloseq object
phyloseq_functional <- phyloseq(otu_table(as.matrix(pathway_final), taxa_are_rows = TRUE),
                                sample_data(metadata_final))



#Ancom-BC2 Analysis:
ancom_results <- ancombc2(data = phyloseq_functional, 
                          fix_formula = "Lauren_classification",
                          p_adj_method = "BH",
                          group = "Lauren_classification")

# Making the results table
ancom_results_table <- ancom_results$res

# Filtering for the significant pathways
sig_pathways <- ancom_results_table %>%
  filter(`q_Lauren_classificationIntestinal type` < 0.05) %>%
  arrange(`q_Lauren_classificationIntestinal type`)

# Saving significant pathways as a CSV
write.csv(sig_pathways, "Aim3_Functional_Results.csv")

# Choosing the 20 most significant pathways so our heatmap isn't clustered
top_pathways <- sig_pathways$taxon[1:20]

# Preparing data for our heatmap, using log-transformation
plot_data <- pathway_final[top_pathways, ]

plot_data_log <- log10(as.matrix(plot_data) + 1)

# Generating our heatmap

# There was an error with NAs, so these lines of code should clear it up
plot_matrix <- as.matrix(plot_data_log)

plot_matrix <- plot_matrix[apply(plot_matrix, 1, var) > 0, ]

# Changing the pathway names
pathway_names <- c(
  "PWY-6654" = "Phosphopantothenate Biosynthesis III",
  "PWY-6167" = "Flavin Biosynthesis II",
  "PWY-6143" = "CMP-Pseudaminate Biosynthesis",
  "PWY-7373" = "Demethylmenaquinol-6 Biosynthesis II",
  "PWY-6305" = "Putrescine Biosynthesis",
  "LIPASYN-PWY" = "Phospholipid Biosynthesis I",
  "PWY-6803" = "Phosphatidylcholine Acyl Editing",
  "RIBOSYN2-PWY" = "Flavin Biosynthesis I",
  "PWY-7357" = "Thiamine Diphosphate Formation",
  "PWY-6897" = "Thiamine Diphosphate Salvage II",
  "PWY-6147" = "Pterin Diphosphate Biosynthesis",
  "PWY0-1241" = "ADP-heptose Biosynthesis",
  "PWY-6703" = "PreQ0 Biosynthesis",
  "NAGLIPASYN-PWY" = "Lipid IVA Biosynthesis (E. coli)",
  "PWY-8073" = "Lipid IVA Biosynthesis (P. putida)",
  "PWY-1269" = "CMP-KDO Biosynthesis",
  "PWY-6467" = "KDO Transfer to Lipid IVA",
  "PWY-7356" = "Thiamine Diphosphate Salvage IV (Yeast)",
  "PWY-I9" = "L-Cysteine Biosynthesis VI",
  "PWY0-162" = "Pyrimidine Ribonucleotide Synthesis",
  
  
  # These are for the sex-based ones!
  "COA-PWY-1" = "Coenzyme A Biosynthesis III",
  "DAPLYSINESYN-PWY" = "L-Lysine Biosynthesis I",
  "GLUCOSE1PMETAB-PWY" = "Glucose & Glucose-1-P Degradation",
  "PENTOSE-P-PWY" = "Pentose Phosphate Pathway",
  "PEPTIDOGLYCANSYN-PWY" = "Peptidoglycan Biosynthesis I",
  "PHOSLIPSYN-PWY" = "Phospholipid Biosynthesis III",
  "POLYISOPRENSYN-PWY" = "Polyisoprenoid Biosynthesis",
  "PWY-2942" = "L-Lysine Biosynthesis III",
  "PWY-5097" = "L-Lysine Biosynthesis VI",
  "PWY-5686" = "UMP Biosynthesis I",
  "PWY-6151" = "S-adenosyl-L-methionine salvage I",
  "PWY-6385" = "Peptidoglycan Biosynthesis III",
  "PWY-6386" = "UDP-N-acetylmuramoyl-pentapeptide Biosynthesis II ",
  "PWY-6387" = "UDP-N-acetylmuramoyl-pentapeptide biosynthesis I",
  "PWY-7221" = "Guanosine Nucleotide Biosynthesis",
  "PWY-724" = "L-Lysine, L-Threonine & L-Methionine Biosynthesis II",
  "ARGDEG-IV-PWY" = "L-Arginine Degradation VIII"
)

rownames(plot_matrix) <- pathway_names[rownames(plot_matrix)]

# This is the annotation
anno_df <- as.data.frame(metadata_final[, "Lauren_classification", drop=FALSE])

heatmap_plot <- pheatmap(plot_matrix,
                         annotation_col = anno_df, 
                         show_colnames = FALSE,
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         scale = "row",
                         main = "Aim 3: Functional Microbiome Signatures: Intestinal vs. Diffused Gastric Cancer (Top 20 Pathways, q < 0.05)",
                         color = colorRampPalette(c("blue", "white", "red"))(100))

# Saving the heatmap plot as a PNG
png("Aim3_Functional_Heatmap.png", width = 16, height = 10, units = "in", res = 800)

print(heatmap_plot)

dev.off()



# Volcano plot (making this cause i feel like just the heatmap might not be enough)
library(ggplot2)
library(ggrepel)

# Data for the volcano plot

volcano_data <- ancom_results_table %>%
  mutate(
    pathway_label = pathway_names[as.character(taxon)],
    neg_log10_q = -log10(`q_Lauren_classificationIntestinal type`),
    lfc=`lfc_Lauren_classificationIntestinal type`,
    Significance = case_when(
      lfc > 0.5 & `q_Lauren_classificationIntestinal type` < 0.05 ~ "Higher in Intestinal",
      lfc < -0.5 & `q_Lauren_classificationIntestinal type` < 0.05 ~ "Higher in Diffused",
      TRUE ~ "Not Significant"
      )
    )

# Creating the volcano plot
volcano_plot <- ggplot(volcano_data, 
                       aes(x = lfc, y = neg_log10_q, color = Significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(data = head(volcano_data %>% arrange(`q_Lauren_classificationIntestinal type`), 4),
                  aes(label = pathway_label), size = 3, color = "black", box.padding = 1, point.padding = 0.5, max.overlaps = Inf) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  scale_color_manual(values = c("Higher in Intestinal" = "salmon", 
                                "Higher in Diffused" = "turquoise3", 
                                "Not Significant" = "grey80")) +
  theme_minimal() +
  labs(title = "Aim 3: Metabolic Pathway Distribution",
       subtitle = "Intestinal vs. Diffused Gastric Cancer (q < 0.05)",
       x = "Log Fold Change (Intestinal / Diffused)",
       y = "-log10(q-value)") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")   


print(volcano_plot)


# Saving the plot
ggsave("Aim3_Volcano_Plot.png", plot = volcano_plot, width = 10, height = 7, dpi = 800)





# LOOKING AT MALE VS FEMALE IN DIFFUSED GROUP

# Filtering phyloseq objext by diffused
phyloseq_diffused <- subset_samples(phyloseq_functional, Lauren_classification == "Diffused type")

# ANCOMBC2 for the biological gender
ancom_sex_results <- ancombc2(data = phyloseq_diffused,
                              fix_formula = "Gender",
                              p_adj_method = "BH",
                              group = "Gender")

sex_results_table <-ancom_sex_results$res

# Filtering for the significant pathways based on Sex (male compared to female)
sig_pathways_sex <- sex_results_table %>%
  filter(q_Gendermale < 0.05) %>%
  arrange(q_Gendermale)

# Saving significant pathways as a CSV
write.csv(sig_pathways_sex, "Aim3_Functional_Sex_Results.csv")


# These are the top 20 for the heatmap
top_pathways_sex <- sig_pathways_sex$taxon[1:20]

# Preparing the heatmap
plot_data_sex <- pathway_final[top_pathways_sex, ]

plot_data_sex_log <- log10(as.matrix(plot_data_sex) + 1)


# Generating the heatmap
plot_matrix_sex <- as.matrix(plot_data_sex_log)

plot_matrix_sex <- plot_matrix_sex[apply(plot_matrix_sex, 1, var) > 0, ]

# Renaming the pathways

rownames(plot_matrix_sex) <- pathway_names[rownames(plot_matrix_sex)]

# This the annotation for gender
anno_sex <- as.data.frame(metadata_final[, "Gender", drop=FALSE])

heatmap_plot_sex <- pheatmap(plot_matrix_sex,
                         annotation_col = anno_sex, 
                         show_colnames = FALSE,
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         scale = "row",
                         main = "Aim 3: Gender-Based Functional Microbiome Signatures in Diffused Subtype (Top 20 Pathways, q < 0.05)",
                         color = colorRampPalette(c("blue", "white", "red"))(100))

# Saving the heatmap plot as a PNG
png("Aim3_Functional_Gender_Heatmap.png", width = 16, height = 10, units = "in", res = 800)

print(heatmap_plot_sex)

dev.off()




# Volcano plot stuff

# Preparing the data
volcano_data_sex <- sex_results_table %>%
  mutate(
    pathway_label = pathway_names[as.character(taxon)], 
    neg_log10_q = -log10(q_Gendermale),
    lfc = lfc_Gendermale,
    Significance = case_when(
      lfc > 0.5 & q_Gendermale < 0.05 ~ "Higher in Male",
      lfc < -0.5 & q_Gendermale < 0.05 ~ "Higher in Female",
      TRUE ~ "Not Significant"
      )
    )

# Creating the volcano plot
volcano_plot_sex <- ggplot(volcano_data_sex,
                           aes(x = lfc, y = neg_log10_q, color = Significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(
    data = head(volcano_data_sex %>% arrange(desc(lfc)), 2),
    aes(label = pathway_label), size = 3.5, color = "black", box.padding = 1.3, 
    point.padding = 0.5, max.overlaps = Inf) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  scale_color_manual(values = c("Higher in Male" = "salmon", 
                                "Higher in Female" = "turquoise", 
                                "Not Significant" = "grey80")) +
  theme_minimal() +
  labs(title = "Aim 3: Functional Pathway Distribution by Biological Gender",
       subtitle = "Comparison within Diffused Type Gastric Cancer (q < 0.05)",
       x = "Log Fold Change (Male / Female)",
       y = "-log10(q-value)") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")

print(volcano_plot_sex)


# saving the volcano plot
ggsave("Aim3_Gender_Volcano_Plot.png", plot = volcano_plot_sex, width = 10, height = 7, dpi = 800)



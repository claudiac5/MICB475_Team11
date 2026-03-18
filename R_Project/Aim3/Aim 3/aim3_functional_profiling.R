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

# This is the annotation
anno_df <- as.data.frame(metadata_final[, "Lauren_classification", drop=FALSE])

pheatmap(plot_matrix,
         annotation_col = anno_df, 
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         main = "Aim 3: Functional Microbiome Signatures: Intestinal vs. Diffused Gastric Cancer",
         color = colorRampPalette(c("blue", "white", "red"))(100))

# Saving the heatmap plot as a PNG
png("Aim3_Functional_Heatmap.png", width = 16, height = 10, units = "in", res = 300)

print(pheatmap(plot_matrix,
               annotation_col = anno_df,
               show_colnames = FALSE,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               scale = "row",
               main = "Aim 3: Functional Microbiome Signatures: Intestinal vs. Diffused Gastric Cancer (Top 20 Pathways, q < 0.05)",
               color = colorRampPalette(c("blue", "white", "red"))(100)))

dev.off()





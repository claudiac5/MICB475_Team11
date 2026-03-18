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

heatmap_plot <- pheatmap(plot_matrix,
                         annotation_col = anno_df, 
                         show_colnames = FALSE,
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         scale = "row",
                         main = "Aim 3: Functional Microbiome Signatures: Intestinal vs. Diffused Gastric Cancer (Top 20 Pathways, q < 0.05)",
                         color = colorRampPalette(c("blue", "white", "red"))(100))

# Saving the heatmap plot as a PNG
png("Aim3_Functional_Heatmap.png", width = 16, height = 10, units = "in", res = 300)

print(heatmap_plot)

dev.off()



# Volcano plot (making this cause i feel like just the heatmap might not be enough)
library(ggplot2)
library(ggrepel)

# Data for the volcano plot

volcano_data <- ancom_results_table %>%
  mutate(neg_log10_q = -log10(`q_Lauren_classificationIntestinal type`),
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
  geom_text_repel(data = head(volcano_data %>% arrange(`q_Lauren_classificationIntestinal type`), 1),
                  aes(label = taxon), size = 3, color = "black", box.padding = 0.5, max.overlaps = Inf) +
  scale_color_manual(values = c("Higher in Intestinal" = "salmon", 
                                "Higher in Diffused" = "turquoise3", 
                                "Not Significant" = "grey80")) +
  theme_minimal() +
  labs(title = "Aim 3: Metabolic Pathway Distribution",
       subtitle = "Intestinal vs. Diffused Gastric Cancer (q < 0.05)",
       x = "Log Fold Change (LFC)",
       y = "-log10(q-value)") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")   


print(volcano_plot)


# Saving the plot
ggsave("Aim3_Volcano_Plot.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)



# Research paper that aligns with our aim 3 findings: 
# id.elsevier.com/as/gh9ccxxZhR/resume/as/authorization.ping
# "Patients with metastatic gastric cancer of the diffuse type have the worst survival."

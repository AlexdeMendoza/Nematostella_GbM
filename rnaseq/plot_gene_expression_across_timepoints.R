# Plot expression of DEGs across timepoints

# Load required packages
library(dplyr)
library(DESeq2)
library(plotR)
library(pheatmap)

# Load a list matching genes to transcripts
genes_to_transcript_list = readRDS("../nematostella_genome/genes_to_transcript_list.rds")

# Load list of methylated and unmethylated genes
unmethylated_and_methylated_genes_list = readRDS("../nematostella_methylation/unmethylated_and_methylated_genes_list.rds")

# Load tables with transcript counts
external_transcript_counts = round(data.frame(data.table::fread("external_rnaseq_kallisto_quantification/kallisto_transcript_counts_merged.tsv.gz", data.table = F), row.names = 1))
june_2024_transcript_counts = round(data.frame(data.table::fread("rnaseq_june_2024/kallisto_quantification/kallisto_transcript_counts_merged.tsv.gz", data.table = F), row.names = 1))
march_2025_transcript_counts = round(data.frame(data.table::fread("rnaseq_march_2025/kallisto_quantification/kallisto_transcript_counts_merged.tsv.gz", data.table = F), row.names = 1))[c(3, 4, 1, 2)]

# Order names of external_transcript_counts 
external_transcript_counts_names_df = data.frame(
  name = names(external_transcript_counts),
  study = gsub("_.*", "", names(external_transcript_counts)),
  timepoint = as.numeric(gsub("Hpf_.*", "", gsub("[CR]_", "", names(external_transcript_counts)))),
  replicate = gsub(".*_", "", names(external_transcript_counts))
  )
external_transcript_counts_names = arrange(external_transcript_counts_names_df, timepoint, study, replicate)$name
external_transcript_counts = external_transcript_counts[external_transcript_counts_names]

# Combine tables
combined_transcript_counts = cbind(external_transcript_counts, june_2024_transcript_counts, march_2025_transcript_counts)

# Combine transcript counts for genes
system.time({nvec_gene_counts = methodical::sumTranscriptValuesForGenes(combined_transcript_counts, genes_to_transcript_list)})

# Create a metadata table for nvec_gene_counts
coldata = data.frame(sample = names(nvec_gene_counts), study = gsub("_.*", "", names(nvec_gene_counts)))
coldata$study[coldata$study == "Nvec"] = "M"
coldata$study[coldata$study == "Nv"] = "Cross"

# Normalize nvec_gene_counts using DESeq2
nvec_gene_counts_dds = DESeqDataSetFromMatrix(countData = nvec_gene_counts, colData = coldata, design = ~ study)
nvec_gene_counts_dds  = estimateSizeFactors(nvec_gene_counts_dds) 
nvec_gene_counts_normalized_deseq = data.frame(counts(nvec_gene_counts_dds, normalized = T))

# Normalize each gene by dividing by its maximum value across samples
nvec_gene_counts_normalized = data.frame(t(apply(nvec_gene_counts_normalized_deseq, 1, function(x) x/max(x))))

# Convert nvec_gene_counts_normalized to long format
nvec_gene_counts_normalized_long = tidyr::pivot_longer(data = tibble::rownames_to_column(nvec_gene_counts_normalized, "gene_id"), 
  -gene_id, names_to = "timepoint", values_to = "normalized_expression")

# Order timepoint and add column indicating if it is in the gastrula stage
nvec_gene_counts_normalized_long$timepoint = factor(nvec_gene_counts_normalized_long$timepoint, 
  levels = c(names(external_transcript_counts)[1:10], names(june_2024_transcript_counts), 
    names(march_2025_transcript_counts), names(external_transcript_counts)[11:28]))
nvec_gene_counts_normalized_long$status = case_when(
  grepl("24Hpf", nvec_gene_counts_normalized_long$timepoint) ~ "Gastrula", 
  grepl("GSK", nvec_gene_counts_normalized_long$timepoint) ~ "GSK",
  grepl("DMSO", nvec_gene_counts_normalized_long$timepoint) ~ "DMSO",
  grepl("Cr._", nvec_gene_counts_normalized_long$timepoint) ~ "GSK F1",
  grepl("Cr.._", nvec_gene_counts_normalized_long$timepoint) ~ "Morpholino F1",
  TRUE ~ "Other")

# Get names of different classes of genes based on whether they were 
# upregulated or downregulated after GSK treatment and whether they are methylated or unmethylated
downregulated_methylated_genes = data.table::fread("deg_beds/DEG_downregulated_methylated.bed")[[4]]
downregulated_unmethylated_genes = data.table::fread("deg_beds/DEG_downregulated_unmethylated.bed")[[4]]
upregulated_methylated_genes = data.table::fread("deg_beds/DEG_upregulated_methylated.bed")[[4]]
upregulated_unmethylated_genes = data.table::fread("deg_beds/DEG_upregulated_unmethylated.bed")[[4]]
nodiff_methylated_genes = data.table::fread("deg_beds/DEG_nodiff_methylated.bed")[[4]]
nodiff_unmethylated_genes = data.table::fread("deg_beds/DEG_nodiff_unmethylated.bed")[[4]]
methylated_genes = c(downregulated_methylated_genes, upregulated_methylated_genes, nodiff_methylated_genes)
unmethylated_genes = c(downregulated_unmethylated_genes, upregulated_unmethylated_genes, nodiff_unmethylated_genes)

# Create a vector for matching genes names to their class
gene_class_list = list(
  upregulated_methylated = upregulated_methylated_genes,
  downregulated_methylated = downregulated_methylated_genes,
  nodiff_methylated = nodiff_methylated_genes,
  upregulated_unmethylated = upregulated_unmethylated_genes,
  downregulated_unmethylated = downregulated_unmethylated_genes,
  nodiff_unmethylated = nodiff_unmethylated_genes
)
gene_class_vec = setNames(rep(names(gene_class_list), sapply(gene_class_list, length)), unlist(gene_class_list))

# Annotate gene class and remove genes which don't belong to any class
nvec_gene_counts_normalized_long$gene_class = gene_class_vec[nvec_gene_counts_normalized_long$gene_id]
nvec_gene_counts_normalized_long = dplyr::filter(nvec_gene_counts_normalized_long, !is.na(gene_class))
nvec_gene_counts_normalized_long$gene_class = factor(nvec_gene_counts_normalized_long$gene_class, levels = c(
  "upregulated_methylated", "upregulated_unmethylated",
  "downregulated_methylated", "downregulated_unmethylated",
  "nodiff_methylated", "nodiff_unmethylated"
))

# Add a column indicating if genes are methylated or unmethylated
nvec_gene_counts_normalized_long$methylation = 
  ifelse(nvec_gene_counts_normalized_long$gene_id %in% unmethylated_and_methylated_genes_list$methylated, "Methylated", "Unmethylated")

# Make boxplots of in external datasets for normalized expression values at different timepoints of different gene groups
nvec_gene_counts_normalized_long_external = filter(nvec_gene_counts_normalized_long, status %in% c("Other", "Gastrula"))
external_boxplots = ggplot(nvec_gene_counts_normalized_long_external, mapping = aes(x = timepoint, y = normalized_expression, fill = status)) + 
  geom_boxplot()
external_boxplots = customize_ggplot_theme(external_boxplots, 
  title = "Transcription Across Developmental Timepoints", xlab = "Timepoint", ylab = "Max Normalized\nDESeq2 Counts", 
  plot_title_size = 30, axis_title_size = 24, legend_key_size = 2, legend_text_size = 24, x_labels_angle = 45) +
  theme(strip.background = element_blank(), plot.margin = margin(t = 5, r = 5, b = 5, l = 20))
external_boxplots
ggsave(plot = external_boxplots, "deseq2_counts_timepoint_external_boxplots.pdf", width = 16, height = 9)

# Facet plots by whether genes are methylated or unmethylated
external_boxplots_faceted_methylation = customize_ggplot_theme(external_boxplots, 
  title = "Transcription Across Developmental Timepoints", xlab = "Timepoint", ylab = "Max Normalized\nDESeq2 Counts", 
  plot_title_size = 30, axis_title_size = 24, legend_key_size = 2, legend_text_size = 24,
  facet = "methylation", x_labels_angle = 45, facet_nrow = 3) +
  theme(strip.background = element_blank(), plot.margin = margin(t = 5, r = 5, b = 5, l = 20))
external_boxplots_faceted_methylation
ggsave(plot = external_boxplots_faceted_methylation, "deseq2_counts_timepoint_external_boxplots_faceted_methylation.pdf", width = 24, height = 13.5)


# Facet plots by gene class
external_boxplots_faceted = customize_ggplot_theme(external_boxplots, 
  title = "Transcription Across Developmental Timepoints", xlab = "Timepoint", ylab = "Max Normalized\nDESeq2 Counts", 
  plot_title_size = 30, axis_title_size = 24, legend_key_size = 2, legend_text_size = 24,
  facet = "gene_class", x_labels_angle = 45, facet_nrow = 3,
  facet_labels = c("Upregulated Methylated", "Upregulated Unmethylated",
  "Downregulated Methylated", "Downregulated Unmethylated",
  "Unchanged Methylated", "Unchanged Unmethylated")) +
  theme(strip.background = element_blank(), plot.margin = margin(t = 5, r = 5, b = 5, l = 20))
external_boxplots_faceted
ggsave(plot = external_boxplots_faceted, "deseq2_counts_timepoint_external_boxplots_faceted.pdf", width = 32, height = 18)

# Make boxplots of in external datasets for normalized expression values in the gastrula stage of different gene groups
nvec_gene_counts_normalized_long_gastrula = dplyr::filter(nvec_gene_counts_normalized_long, !status %in% c("Gastrula", "Other"))
gastrula_boxplots = ggplot(nvec_gene_counts_normalized_long_gastrula, mapping = aes(x = timepoint, y = normalized_expression, fill = status)) + 
  geom_boxplot()
gastrula_boxplots = customize_ggplot_theme(gastrula_boxplots, 
  title = "Impact of DNMT1 Inhibition on Transcription in Gastrula", xlab = "Timepoint", ylab = "Max Normalized\nDESeq2 Counts", 
  fill_labels = c("DMSO", "GSK", "GSK F1", "Morpholino F1"), 
  plot_title_size = 30, axis_title_size = 24, legend_key_size = 2, legend_text_size = 24,
  facet = "gene_class", x_labels_angle = 45, facet_nrow = 3, fill_colors = plotR::colour_list$rainbox_six[c(6, 4, 1, 2)],
  facet_labels = c("Upregulated Methylated", "Upregulated Unmethylated",
  "Downregulated Methylated", "Downregulated Unmethylated",
  "Unchanged Methylated", "Unchanged Unmethylated")) +
  theme(strip.background = element_blank(), plot.margin = margin(t = 5, r = 5, b = 5, l = 20))
gastrula_boxplots
ggsave(plot = gastrula_boxplots, "deseq2_counts_timepoint_gastrula_boxplots.pdf", width = 32, height = 18)

# Make a clustered heatmap for the gastrula samples
heatmap_df = nvec_gene_counts_normalized[grep("24[^0]|25[^0]", names(nvec_gene_counts_normalized))] 
heatmap_df = heatmap_df[complete.cases(heatmap_df), ]
annotation = data.frame(tibble::column_to_rownames(unique(select(nvec_gene_counts_normalized_long_gastrula, timepoint, status)), "timepoint"))
names(annotation) = "Status"
annotation$Status[annotation$Status == "Gastrula"] = "External"
annotation_colors = list(Status = setNames(plotR::colour_list$rainbox_six[c(6, 5, 4, 1, 2)], sort(unique(annotation$Status))))
set.seed(123)
pheatmap(heatmap_df, kmeans_k = 100, show_rownames = F, clustering_method = "ward.D", 
  treeheight_row = 0, annotation_col = annotation, annotation_names_col = F, 
  annotation_colors = annotation_colors, fontsize = 14, fontsize_col = 10, 
  main = "Clustering of Samples by RNA-Seq", filename = "gastrula_samples_clustered_heatmap.pdf", width = 16, height = 9)

# Make boxplots of normalized gene expression across timepoints for different categories of differentially expressed and differentially methylated genes

nvec_gene_counts_normalized_boxplots = ggplot(nvec_gene_counts_normalized_long, mapping = aes(x = timepoint, y = normalized_expression, fill = status)) + 
  geom_boxplot()
nvec_gene_counts_normalized_boxplots = customize_ggplot_theme(nvec_gene_counts_normalized_boxplots, xlab = "Timepoint", ylab = "Max Normalized\nDESeq2 Counts", 
  facet = "gene_class", x_labels_angle = 45, facet_nrow = 3, fill_colors = plotR::colour_list$rainbox_six[c(6, 5, 4, 1, 2)],
  facet_labels = c("Upregulated Methylated", "Upregulated Unmethylated",
  "Downregulated Methylated", "Downregulated Unmethylated",
  "Unchanged Methylated", "Unchanged Unmethylated")) +
  theme(strip.background = element_blank(), plot.margin = margin(t = 5, r = 5, b = 5, l = 20))
nvec_gene_counts_normalized_boxplots
ggsave(plot = nvec_gene_counts_normalized_boxplots, "plots/deseq2_counts_timepoint_nvec_gene_counts_normalized_boxplots_with_crossings.pdf", width = 32, height = 18)

# Make a correlation heatmap of our samples
plotR::heatmap_without_clustering(cor(select(nvec_gene_counts_normalized, starts_with("N")), use = "c"), 
  filename = "plots/gsk_dmso_samples_transcript_cor_heatmap.pdf", file_dimensions = c(16, 16))

# Make a dendogram of our samples
pdf("plots/gsk_dmso_samples_transcript_cor_dendogram.pdf", width = 16, height = 9)
plot(hclust(as.dist(1 - cor(select(nvec_gene_counts_normalized, starts_with("N")), use = "c")), method = "ward.D"))
dev.off()

### Test differential expression for DMSO Vs GSK, GSK F1 and Morpholino F1

# Select our samples
de_nvec_gene_counts = select(nvec_gene_counts, starts_with("N"))

# Remove Nvec_24hpf_GSK_rep3 so that there are two samples each for GSK, GSK F1 and Morpholino F1
de_nvec_gene_counts = select(de_nvec_gene_counts, -Nvec_24hpf_GSK_rep3)

# Create a metadata table for nvec_gene_counts
de_coldata = data.frame(sample = names(de_nvec_gene_counts), 
  treatment = c(rep("DMSO", 3), rep("GSK", 2), rep("Morpholino_F1", 2), rep("GSK_F1", 2)))

# Normalize nvec_gene_counts using DESeq2
nvec_gene_counts_dds = DESeqDataSetFromMatrix(countData = de_nvec_gene_counts, colData = de_coldata, design = ~ treatment)
nvec_gene_counts_dds  = estimateSizeFactors(nvec_gene_counts_dds) 
system.time({nvec_gene_counts_dds = DESeq(nvec_gene_counts_dds)})

# Get results for GSK, GSK F1 and Morpholino F1
deseq_results_gsk_vs_dmso = results(nvec_gene_counts_dds, contrast = c("treatment", "GSK", "DMSO"), alpha = 0.05)
deseq_results_morpholino_f1_vs_dmso = results(nvec_gene_counts_dds, contrast = c("treatment", "Morpholino_F1", "DMSO"), alpha = 0.05)
deseq_results_gsk_f1_vs_dmso = results(nvec_gene_counts_dds, contrast = c("treatment", "GSK_F1", "DMSO"), alpha = 0.05)

# Create a list with upregulated and downregulated genes for each condition and create and Upset plot
de_list = list(
  gsk_vs_dmso = row.names(deseq_results_gsk_vs_dmso[which(deseq_results_gsk_vs_dmso$padj < 0.05), ]),
  morpholino_f1_vs_dmso = row.names(deseq_results_morpholino_f1_vs_dmso[which(deseq_results_morpholino_f1_vs_dmso$padj < 0.05), ]),
  gsk_f1_vs_dmso = row.names(deseq_results_gsk_f1_vs_dmso[which(deseq_results_gsk_f1_vs_dmso$padj < 0.05), ])
  )
make_upset_plot(input_list = de_list, filename = "plots/de_upset_plot.jpeg", title = "Upregulated and Downregulated Genes")

# Make a data.frame with the log2 fold change for each gene in each condition along with the associated q-value
de_results_df = data.frame(
  gene = row.names(deseq_results_gsk_vs_dmso),
  gsk_vs_dmso_log2fc = deseq_results_gsk_vs_dmso$log2FoldChange,
  gsk_vs_dmso_qval = deseq_results_gsk_vs_dmso$padj,
  gsk_f1_vs_dmso_log2fc = deseq_results_gsk_f1_vs_dmso$log2FoldChange,
  gsk_f1_vs_dmso_qval = deseq_results_gsk_f1_vs_dmso$padj,
  morpholino_f1_vs_dmso_log2fc = deseq_results_morpholino_f1_vs_dmso$log2FoldChange,
  morpholino_f1_vs_dmso_qval = deseq_results_morpholino_f1_vs_dmso$padj)

# Add a column indicating if genes are methylated in control samples
de_results_df$methylation = case_when(
  de_results_df$gene %in% unmethylated_and_methylated_genes_list$methylated ~ "Methylated",
  de_results_df$gene %in% unmethylated_and_methylated_genes_list$unmethylated ~ "Unmethylated",
  TRUE ~ NA)
de_results_df = filter(de_results_df, !is.na(methylation))

# Create tables for DEGs in GSK or GSKF1 and DEGs in GSK or morpholino
degs_gsk_gskf1_df = filter(de_results_df, gsk_vs_dmso_qval < 0.05 & abs(gsk_vs_dmso_log2fc) >= 1 | gsk_f1_vs_dmso_qval < 0.05 & abs(gsk_f1_vs_dmso_log2fc) >= 1)
degs_gsk_morpholino_df = filter(de_results_df, gsk_vs_dmso_qval < 0.05  & abs(gsk_vs_dmso_log2fc) >= 1 | morpholino_f1_vs_dmso_qval < 0.05  & abs(morpholino_f1_vs_dmso_log2fc) >= 1)

# Create a scatter plot showing the log2FC between GSK/DMSO vs GSKF1/DMSO
gsk_vs_gskf1_scatter = ggplot(de_results_df, aes(x = gsk_vs_dmso_log2fc, y = gsk_f1_vs_dmso_log2fc, fill = methylation)) +
  geom_point(alpha = 0.5, shape = 21, size = 5) +
  geom_text(label = paste("Pearson's r =", round(cor(de_results_df$gsk_vs_dmso_log2fc, de_results_df$gsk_f1_vs_dmso_log2fc, use = "c"),  2)), x = 0, y = -20, size = 5)
gsk_vs_gskf1_scatter = customize_ggplot_theme(gsk_vs_gskf1_scatter, xlab = "GSK vs DMSO log2FC", ylab = "GSK F1 vs DMSO log2FC")
gsk_vs_gskf1_scatter = ggExtra::ggMarginal(gsk_vs_gskf1_scatter, type = "histogram")
gsk_vs_gskf1_scatter
ggsave(plot = gsk_vs_gskf1_scatter, filename = "plots/gsk_and_gsk_f1_log2fc_scatter_plot.pdf", width = 16, height = 9)

# Create a scatter plot showing the log2FC between GSK/DMSO vs morpholino/DMSO
gsk_morpholino_scatter = ggplot(de_results_df, aes(x = gsk_vs_dmso_log2fc, y = morpholino_f1_vs_dmso_log2fc, fill = methylation)) +
  geom_point(alpha = 0.5, shape = 21, size = 5) +
  geom_text(label = paste("Pearson's r =", round(cor(de_results_df$gsk_vs_dmso_log2fc, de_results_df$morpholino_f1_vs_dmso_log2fc, use = "c"), 2)), x = 0, y = -20, size = 5)
gsk_morpholino_scatter = customize_ggplot_theme(gsk_morpholino_scatter, xlab = "GSK vs DMSO log2FC", ylab = "Morpholino F1 vs DMSO log2FC")
gsk_morpholino_scatter = ggExtra::ggMarginal(gsk_morpholino_scatter, type = "histogram")
gsk_morpholino_scatter
ggsave(plot = gsk_morpholino_scatter, filename = "plots/gsk_and_morpholino_f1_log2fc_scatter_plot.pdf", width = 16, height = 9)
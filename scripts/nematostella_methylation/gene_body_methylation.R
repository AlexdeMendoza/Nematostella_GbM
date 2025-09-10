# Calculate methylation of protein-coding genes in Nematostella and define unmethylated and methylated genes

# Load required packages
library(methodical)
library(dplyr)
source("../auxillary_scripts/plotting_functions.R")

# Load methrix for Nematostella and convert to a meth RSE
nvec_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("nematostella_complete_meth_rse/")

# Get GRanges with genes from Nematostella
nvev_pc_genes_gr = readRDS("~/nematostella_project/nematostella_genome/nvec_pc_genes_gr.rds")

# Get mean methylation for each gene. 
bpparam = BiocParallel::MulticoreParam(2)
system.time({gene_meth = summarizeRegionMethylation(meth_rse = nvec_meth_rse, genomic_regions = nvev_pc_genes_gr, 
  genomic_region_names = names(nvev_pc_genes_gr), BPPARAM = bpparam)})
data.table::fwrite(tibble::rownames_to_column(gene_meth, "gene_id"), "nvec_gene_meth_df.tsv.gz", sep = "\t")

# Filter for DMSO treated samples and get mean methylation of each gene. 63 genes are missing methylation in all DMSO samples
mean_gene_meth_dmso_df = data.frame(mean_meth = rowMeans(select(gene_meth, starts_with("DMSO")), na.rm = T))

# Create a density plot of gene methylation in DMSO samples
gene_meth_density_plot_dmso = ggplot(data.frame(mean_meth = mean_gene_meth_dmso_df), aes(x = mean_meth)) +
  geom_density(fill = "#7B5C90", color = "black") 
gene_meth_density_plot_dmso = customize_ggplot_theme(gene_meth_density_plot_dmso, 
  xlab = "Mean Methylation", ylab = "Density", title = "Gene Body Methylation in Nematostella") + 
  scale_x_log10(labels = scales::label_number(drop0trailing=TRUE), expand = expansion(mult = c(0, 0), add = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01), add = c(0, 0))) + 
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  theme(plot.margin = margin(r = 20))
gene_meth_density_plot_dmso 
ggsave(plot = gene_meth_density_plot_dmso, "nematostella_gene_meth_density_plot_dmso.pdf", width = 16, height = 9)

# Use 0.1 as a threshold to define methylated and unmethylated genes
unmethylated_genes = row.names(filter(mean_gene_meth_dmso_df, mean_meth < 0.1))
methylated_genes = row.names(filter(mean_gene_meth_dmso_df, mean_meth >= 0.1))

# Create a list with unmethylated and methylated genes and save
unmethylated_and_methylated_genes_list = list(
  unmethylated = unmethylated_genes,
  methylated = methylated_genes
)
saveRDS(unmethylated_and_methylated_genes_list, "unmethylated_and_methylated_genes_list.rds")

# Load a list of genes and their associated transcripts
genes_to_transcript_list = readRDS("~/genomes/nematostella/genes_to_transcript_list.rds")

# Create a list with transcripts associated with unmethylated and methylated genes and save
unmethylated_and_methylated_transcripts_list = lapply(unmethylated_and_methylated_genes_list, function(x)
  unlist(genes_to_transcript_list[x]))
saveRDS(unmethylated_and_methylated_transcripts_list, "unmethylated_and_methylated_transcripts_list.rds")

# Create a data.frame with the number of unmethylated and methylated genes
methylation_status_df = data.frame(
  status = names(unmethylated_and_methylated_genes_list),
  number = lengths(unmethylated_and_methylated_genes_list)
  )

# Make a barplot of the number of unmethylated and methylated genes
gene_barplot = ggplot(methylation_status_df, aes(x = status, y = number, fill = status)) + geom_col(color = "black")
gene_barplot = customize_ggplot_theme(gene_barplot, title = "Number of Methylated and Unmethylated Genes", 
  x_labels = c("Methylated", "Unmethylated"), fill_colors = c("#CD2626", "#53868B"), show_legend = F)
gene_barplot
ggsave(plot = gene_barplot, "methylation_status_barplot.pdf", width = 9, height = 9)

# Filter for GSK treated samples and get mean methylation of each gene
mean_gene_meth_gsk_df = data.frame(mean_meth = rowMeans(select(gene_meth, starts_with("gsk")), na.rm = T))

# Combine mean_gene_meth_dmso_df and mean_gene_meth_gsk_df
combined_mean_gene_meth_df = bind_rows(list(
  dmso = mean_gene_meth_dmso_df,
  gsk = mean_gene_meth_gsk_df),
  .id = "treatment"
)

# Create a density plot of gene methylation in GSK samples compared to DMSO samples
gene_meth_density_plot_gsk = ggplot(combined_mean_gene_meth_df, aes(x = mean_meth, fill = treatment)) +
  geom_density(color = "black", alpha = 0.4) 
gene_meth_density_plot_gsk = customize_ggplot_theme(gene_meth_density_plot_gsk, 
  xlab = "Mean Methylation", ylab = "Density", title = "Gene Body Methylation in Nematostella After GSK Treatment",
  fill_colors = c("#7B5C90", "#bfab25"), fill_labels = c("DMSO", "GSK")) + 
  scale_x_log10(labels = scales::label_number(drop0trailing=TRUE), expand = expansion(mult = c(0, 0), add = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01), add = c(0, 0))) + 
  theme(plot.margin = margin(r = 20))
gene_meth_density_plot_gsk 
ggsave(plot = gene_meth_density_plot_gsk, "nematostella_gene_meth_density_plot_gsk.pdf", width = 16, height = 9)

# Turn unmethylated_and_methylated_genes_list inside out
unmethylated_and_methylated_genes_list = unlist(with(reshape2::melt(readRDS("unmethylated_and_methylated_genes_list.rds")), split(L1, value)))

# Load table with raw gene counts and convert to TPM
complete_gene_counts_raw = data.frame(data.table::fread("rnaseq/nvec_complete_gene_counts_raw.tsv.gz"), row.names = 1)
complete_gene_counts_tpm = data.frame(apply(complete_gene_counts_raw, 2, function(x) x/sum(x)*1e6))

# Create data.frames with the mean TPM in DMSO and GSK samples
dmso_mean_gene_counts_tpm_df = data.frame(
  gene = row.names(complete_gene_counts_tpm),
  mean_tpm = rowMeans(select(complete_gene_counts_tpm, starts_with("DMSO"))),
  methylated_status = unmethylated_and_methylated_genes_list[row.names(complete_gene_counts_tpm)],
  row.names = NULL)
dmso_mean_gene_counts_tpm_df = dmso_mean_gene_counts_tpm_df[complete.cases(dmso_mean_gene_counts_tpm_df), ]
gsk_mean_gene_counts_tpm_df = data.frame(
  gene = row.names(complete_gene_counts_tpm),
  mean_tpm = rowMeans(select(complete_gene_counts_tpm, starts_with("gsk"))),
  methylated_status = unmethylated_and_methylated_genes_list[row.names(complete_gene_counts_tpm)],
  row.names = NULL)
gsk_mean_gene_counts_tpm_df = gsk_mean_gene_counts_tpm_df[complete.cases(gsk_mean_gene_counts_tpm_df), ]

# Combine tables
combined_mean_gene_counts_tpm_df = bind_rows(list(
  dmso = dmso_mean_gene_counts_tpm_df,
  gsk = gsk_mean_gene_counts_tpm_df),
  .id = "treatment")

# Create boxplots for mean log10 TPM values in DMSO samples
dmso_gene_counts_boxplot = ggplot(dmso_mean_gene_counts_tpm_df, aes(x = methylated_status, y = log10(mean_tpm), fill = methylated_status)) +
  geom_boxplot()
dmso_gene_counts_boxplot = customize_ggplot_theme(dmso_gene_counts_boxplot, 
  title = "Expression of Unmethylated and Methylated Genes", ylab = "Log10 TPM", x_labels = c("Methylated", "Unmethylated"), 
  fill_colors = c("#CD2626", "#53868B"), show_legend = F)
dmso_gene_counts_boxplot
ggsave(plot = dmso_gene_counts_boxplot, "dmso_gene_counts_boxplot.pdf", width = 9, height = 9)

# Create boxplots for mean log10 TPM values in DMSO and GSK samples
combined_gene_counts_boxplot = ggplot(combined_mean_gene_counts_tpm_df, aes(x = methylated_status, y = log10(mean_tpm), fill = treatment)) +
  geom_boxplot()
combined_gene_counts_boxplot = customize_ggplot_theme(combined_gene_counts_boxplot, 
  title = "Expression of Unmethylated and Methylated Genes After GSK Treatment", ylab = "Log10 TPM", x_labels = c("Methylated", "Unmethylated"), 
  fill_colors = c("#7B5C90", "#bfab25"), fill_labels = c("DMSO", "GSK"))
combined_gene_counts_boxplot
ggsave(plot = combined_gene_counts_boxplot, "combined_gene_counts_boxplot.pdf", width = 16, height = 9)

# Load table with differential expression results and indicate whether genes are methylated or not and upregulated or downregulated
deseq_gsk_vs_dmso_significant_results = data.table::fread("rnaseq/deseq_gsk_vs_dmso_significant_results.tsv.gz")
deseq_gsk_vs_dmso_significant_results$direction = ifelse(deseq_gsk_vs_dmso_significant_results$log2FoldChange > 0, "Upregulated", "Downregulated")
deseq_gsk_vs_dmso_significant_results$methylation_status = unmethylated_and_methylated_genes_list[deseq_gsk_vs_dmso_significant_results$gene_id]
deseq_gsk_vs_dmso_significant_results = deseq_gsk_vs_dmso_significant_results[complete.cases(deseq_gsk_vs_dmso_significant_results), ]

# Create a data.frame summarizing the DEG results
deseq_gsk_vs_dmso_significant_results_summary = summarize(group_by(deseq_gsk_vs_dmso_significant_results, direction, methylation_status),
  count = n())

# Create a barplot for the number of upregulated and downreulated genes by methylation status
deg_barplot = ggplot(deseq_gsk_vs_dmso_significant_results_summary, aes(x = direction, y = count, fill = methylation_status)) +
  geom_col(position = "dodge", color = "black")
deg_barplot = customize_ggplot_theme(deg_barplot, ylab = "Number of Genes", title = "Differential Expression Results After GSK Treatment",
  fill_colors = c("#CD2626", "#53868B"), fill_labels = c("Methylated", "Unmethylated"))
deg_barplot
ggsave(plot = deg_barplot, "deg_barplot.pdf", width = 16, height = 9)
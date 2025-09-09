# Fit linear models to methylation recovery in Nematostella

# Load required packages
library(dplyr)
library(methodical)
source("../auxillary_scripts/plotting_functions.R")
source("../auxillary_scripts/bigwig_summarize_over_regions.R")

# Load DMRs
dmr_gr = readRDS("../nematostella_methylation/dmso_dmrs.rds")

# Get GRanges with TSS
tss_gr = readRDS("~/nematostella_project/nematostella_genome/nvec_pc_transcripts_tss_gr.rds")

# Load protein-coding genes and split into 200 bp regions
gene_regions_gr = readRDS("../nematostella_genome/nvec_pc_genes_gr.rds")
gene_regions_gr = reduce(gene_regions_gr, ignore.strand = T)
gene_regions_gr_bins = unlist(tile(gene_regions_gr, width = 200))
names(gene_regions_gr_bins) = paste0("region_", seq_along(gene_regions_gr_bins))

# Indicate if gene_regions_gr_bins overlap DMRs and what is their status. 
# A small number of regions overlap slow and fast DMRs so their status is set to fast
gene_regions_gr_bins$dmr_status = "No DMR"
gene_regions_gr_bins$dmr_status[gene_regions_gr_bins %over% dmr_gr[dmr_gr$dmr_status == "Slow"]] = "Slow"
gene_regions_gr_bins$dmr_status[gene_regions_gr_bins %over% dmr_gr[dmr_gr$dmr_status == "Fast"]] = "Fast"

# Get paths to bigwig files for histone marks and extract name of marks
histone_bigwigs = list.files("navarrete_bigwigs", full.names = T)
histone_marks = gsub(".log2r.input.bw", "", gsub(".*Nvec_", "", histone_bigwigs))

# Load meth RSE object for Nematostella
nematostella_complete_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../nematostella_methylation/nematostella_complete_meth_rse")

# Define Treated F0 gastrula and F0 sperm samples for calculating mean methylation recovery
gsk_f0_samples = c("GSK_rep1", "GSK_rep2", "GSK_rep3")
gsk_sperm_F0_samples = c("NvIT17B14_sperm", "NvIT17B19_sperm", "NvIT17C6_sperm", "NvIT16K4_sperm_Oct31")

# Summarize histone mark values for each DMR in dmr_gsk_vs_dmso_gr. Took 20 seconds
system.time({gene_regions_histone_values = bigwig_summarize_over_regions(bw_filepaths = histone_bigwigs, 
  gr = gene_regions_gr_bins, column_names = histone_marks, statistic = "mean0")})

# Add distance to nearest TSS for each DMR
gene_regions_histone_values$`TSS Distance` = distance(gene_regions_gr_bins, tss_gr[nearest(gene_regions_gr_bins, tss_gr)])

# Get mean expression of transcripts in DMSO samples
nvec_complete_transcript_counts_normalized = rowMeans(select(
  data.frame(data.table::fread("~/nematostella_project/rnaseq/nvec_complete_transcript_counts_normalized.tsv.gz"), row.names = 1), starts_with("DMSO")), na.rm = T)

# Add activity of nearest TSS to each DMR
gene_regions_histone_values$`TSS Activity` = nvec_complete_transcript_counts_normalized[names(tss_gr[nearest(gene_regions_gr_bins, tss_gr)])]

# Get mean methylation for promoters, exons and introns in all samples. Took 10 seconds
system.time({gene_regions_methylation = summarizeRegionMethylation(meth_rse = nematostella_complete_meth_rse, 
  genomic_regions = gene_regions_gr_bins, BPPARAM = BiocParallel::SerialParam())})

# Calculate the proportion of methylation recovered in regions in sperm
gene_regions_histone_values$gene_regions_meth_recovery = rowMeans(dplyr::select(gene_regions_methylation, all_of(gsk_sperm_F0_samples)), na.rm = T) - 
  rowMeans(dplyr::select(gene_regions_methylation, all_of(gsk_f0_samples)), na.rm = T)

# Remove any NA values
gene_regions_histone_values = gene_regions_histone_values[complete.cases(gene_regions_histone_values), ]

# Make a heatmap of correlation values
heatmap_without_clustering(cor(gene_regions_histone_values))

# Fit a linear model using all the histone marks and save
gene_regions_lm = lm(gene_regions_meth_recovery ~ ., gene_regions_histone_values)
saveRDS(gene_regions_lm, "gene_regions_meth_histone_mark_lm.rds")
gene_regions_lm = readRDS("gene_regions_meth_histone_mark_lm.rds")

# Get the adjusted R-squared (0.61)
gene_regions_lm_adj_r_squared = round(summary(gene_regions_lm)$adj.r.squared, 2)

# Get the relative importance of the features as a data.frame
gene_regions_lm_feature_importance = sort(relaimpo::calc.relimp(gene_regions_lm, type = "lmg", rela = TRUE)$lmg, decreasing = T)
gene_regions_lm_feature_importance = data.frame(histone_mark = names(gene_regions_lm_feature_importance), 
  importance = gene_regions_lm_feature_importance, row.names = NULL)
gene_regions_lm_feature_importance$histone_mark = factor(gene_regions_lm_feature_importance$histone_mark, levels = gene_regions_lm_feature_importance$histone_mark)

# Create a barplot with the relative feature importances
gene_regions_lm_feature_importance_barplot = ggplot(gene_regions_lm_feature_importance, aes(x = histone_mark, y = importance)) +
  geom_col(color = "black", fill = "#A50F15")
gene_regions_lm_feature_importance_barplot = customize_ggplot_theme(gene_regions_lm_feature_importance_barplot, 
  xlab = "Histone Mark", ylab = "Relative Importance", x_labels_angle = 45, show_legend = F)
gene_regions_lm_feature_importance_barplot
ggsave(plot = gene_regions_lm_feature_importance_barplot, "gene_regions_lm_feature_importance_barplot.pdf", width = 9, height = 9)

# Create a data.frame with the the actual vs predicted methylation recovery values and add gene_regions status
gene_regions_lm_predictions_df = data.frame(predicted_recovery = predict(gene_regions_lm), actual_recovery = gene_regions_histone_values$gene_regions_meth_recovery)
gene_regions_lm_predictions_df$dmr_status = gene_regions_gr_bins[names(predict(gene_regions_lm))]$dmr_status

# Sample 10,000 random observations from gene_regions_lm_predictions_df
set.seed(123)
gene_regions_lm_predictions_df_sample = gene_regions_lm_predictions_df[sample(nrow(gene_regions_lm_predictions_df), 20000), ]

# create a scatter plot of the actual vs predicted methylation recovery values
gene_regions_lm_plot = ggplot(gene_regions_lm_predictions_df_sample, aes(x = predicted_recovery, y = actual_recovery)) +
  geom_point(mapping = aes(fill = dmr_status), shape = 21, size = 2, alpha = 0.5)
gene_regions_lm_plot = customize_ggplot_theme(gene_regions_lm_plot, xlab = "Predicted Recovery", ylab = "Actual Recovery",
  fill_colors = c("#F8766D", "white", "#00BFC4")) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") + 
  geom_text(label = paste("Adj. R-squared =", gene_regions_lm_adj_r_squared), x = 0.7, y = -0.5, size = 4)
gene_regions_lm_plot
ggsave(plot = gene_regions_lm_plot, "gene_regions_lm_actual_vs_predicted_plot.pdf", width = 9, height = 9)

### Fit a random forest instead of a linear model

# Fit random forest. Has an R-squared value of 0.55
set.seed(123)
system.time({dmr_rf = ranger::ranger(dependent.variable.name = "gene_regions_meth_recovery", data = gene_regions_histone_values, importance = "permutation")})

# Create a data.frame with variable importance for features
rf_variable_importance_df = data.frame(
  feature = names(dmr_rf$variable.importance),
  importance = dmr_rf$variable.importance,
  row.names = NULL
)

# Add column with importance scores normalized by the maximum score
rf_variable_importance_df$importance_normalized = rf_variable_importance_df$importance/max(rf_variable_importance_df$importance)

# Arrange by importance 
rf_variable_importance_df = arrange(rf_variable_importance_df, desc(importance_normalized))
rf_variable_importance_df$feature = factor(rf_variable_importance_df$feature, levels = rf_variable_importance_df$feature)

# Plot importance
rf_variable_importance_plot = ggplot(rf_variable_importance_df, aes(x = feature, y = importance_normalized)) +
  geom_col(color = "black", fill = "#A50F15")
rf_variable_importance_plot =  customize_ggplot_theme(rf_variable_importance_plot, xlab = "Variable Name", ylab = "Variable Importance", x_labels_angle = 45)
rf_variable_importance_plot

### Repeat for DMRs

# Summarize histone mark values for each DMR in dmr_gsk_vs_dmso_gr. Took 20 seconds
system.time({dmr_histone_values = genomicTools::bigwig_summarize_over_regions(bw_filepaths = histone_bigwigs, 
  gr = dmr_gr, column_names = histone_marks, statistic = "mean0")})

# Get mean methylation for DMRs
system.time({dmr_methylation = summarizeRegionMethylation(meth_rse = nematostella_complete_meth_rse, 
  genomic_regions = dmr_gr, BPPARAM = BiocParallel::SerialParam())})

# Calculate the proportion of methylation recovered in regions in sperm
dmr_histone_values$dmr_meth_recovery = rowMeans(select(dmr_methylation, all_of(gsk_sperm_F0_samples)), na.rm = T) - 
  rowMeans(select(dmr_methylation, all_of(gsk_f0_samples)), na.rm = T)

# Add distance to nearest TSS for each DMR
dmr_histone_values$`TSS Distance` = distance(dmr_gr, tss_gr[nearest(dmr_gr, tss_gr)])

# Add activity of nearest TSS to each DMR
dmr_histone_values$`TSS Activity` = nvec_complete_transcript_counts_normalized[names(tss_gr[nearest(dmr_gr, tss_gr)])]

# Remove any NA values
dmr_histone_values = dmr_histone_values[complete.cases(dmr_histone_values), ]

# Fit a random forest. Has an R-squared value of 0.55
set.seed(123)
system.time({dmr_rf = ranger::ranger(dependent.variable.name = "dmr_meth_recovery", data = dmr_histone_values, importance = "permutation")})

# Create a data.frame with variable importance for features
rf_variable_importance_df = data.frame(
  feature = names(dmr_rf$variable.importance),
  importance = dmr_rf$variable.importance,
  row.names = NULL
)

# Add column with importance scores normalized by the maximum score
rf_variable_importance_df$importance_normalized = rf_variable_importance_df$importance/max(rf_variable_importance_df$importance)

# Arrange by importance 
rf_variable_importance_df = arrange(rf_variable_importance_df, desc(importance_normalized))
rf_variable_importance_df$feature = factor(rf_variable_importance_df$feature, levels = rf_variable_importance_df$feature)

# Plot importance
rf_variable_importance_plot = ggplot(rf_variable_importance_df, aes(x = feature, y = importance_normalized)) +
  geom_col(color = "black", fill = "#A50F15")
rf_variable_importance_plot =  customize_ggplot_theme(rf_variable_importance_plot, xlab = "Variable Name", ylab = "Variable Importance", x_labels_angle = 45)
rf_variable_importance_plot
ggsave(plot = rf_variable_importance_plot, "dmrs_rf_feature_importance_barplot.pdf", width = 9, height = 9)

# Fit a linear model using all the histone marks
dmr_lm = lm(dmr_meth_recovery ~ ., dmr_histone_values)
saveRDS(dmr_lm, "dmr_meth_histone_mark_lm.rds")
dmr_lm = readRDS("dmr_meth_histone_mark_lm.rds")

# Get the adjusted R-squared. 0.53
dmr_lm_adj_r_squared = round(summary(dmr_lm)$adj.r.squared, 2)

# Get the relative importance of the features as a data.frame
dmrs_lm_feature_importance = sort(relaimpo::calc.relimp(dmr_lm, type = "lmg", rela = TRUE)$lmg, decreasing = T)
dmrs_lm_feature_importance = data.frame(histone_mark = names(dmrs_lm_feature_importance), 
  importance = dmrs_lm_feature_importance, row.names = NULL)

# Add correlation of histone mark with methylation gain
dmrs_lm_feature_importance$correlation = cor(dmr_histone_values, use = "p")["dmr_meth_recovery", dmrs_lm_feature_importance$histone_mark]

dmrs_lm_feature_importance$histone_mark = factor(dmrs_lm_feature_importance$histone_mark, levels = dmrs_lm_feature_importance$histone_mark)

# Create a barplot with the relative feature importances
dmrs_lm_feature_importance_barplot = ggplot(dmrs_lm_feature_importance, aes(x = histone_mark, y = importance, fill = correlation)) +
  geom_col(color = "black")
dmrs_lm_feature_importance_barplot = customize_ggplot_theme(dmrs_lm_feature_importance_barplot, 
  xlab = "Histone Mark", ylab = "Relative Importance", x_labels_angle = 45, show_legend = T) + scale_fill_gradient2(low = "#a61919ff", mid = "white", high = "#2b2e83ff")
dmrs_lm_feature_importance_barplot
ggsave(plot = dmrs_lm_feature_importance_barplot, "dmrs_lm_feature_importance_barplot.pdf", width = 9, height = 9)

# Create a data.frame with the the actual vs predicted methylation recovery values and add DMR status
dmr_lm_predictions_df = data.frame(predicted_recovery = predict(dmr_lm), actual_recovery = dmr_histone_values$dmr_meth_recovery)
dmr_lm_predictions_df$dmr_status = dmr_gr[names(predict(dmr_lm))]$dmr_status

# Sample 5,000 random observations from dmr_lm_predictions_df
set.seed(123)
dmr_lm_predictions_df_sample = dmr_lm_predictions_df[sample(nrow(dmr_lm_predictions_df), 5000), ]

# create a scatter plot of the actual vs predicted methylation recovery values
dmr_lm_plot = ggplot(dmr_lm_predictions_df_sample, aes(x = predicted_recovery, y = actual_recovery, fill = dmr_status)) +
  geom_point(shape = 21, size = 2, alpha = 0.5)
dmr_lm_plot = customize_ggplot_theme(dmr_lm_plot, xlab = "Predicted Recovery", ylab = "Actual Recovery") +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") + 
  geom_text(label = paste("Adj. R-squared =", dmr_lm_adj_r_squared), x = 0.7, y = -0.5, size = 6) +
  scale_y_continuous(limits = c(-0.5, 1))
dmr_lm_plot
ggsave(plot = dmr_lm_plot, "dmr_lm_actual_vs_predicted_plot_sample.pdf", width = 9, height = 9)
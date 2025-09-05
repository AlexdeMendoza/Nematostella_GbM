# Test differential methylation between different conditions

# Load required packages
library(DSS)
library(methodical)
library(ComplexHeatmap)
library(dplyr)
library(DESeq2)
library(plotR)
source("../auxillary_scripts/granges_functions.R")
source("../auxillary_scripts/bigwig_summarize_over_regions.R")

### Identify DMRs between DMSO and GSK treated samples

# Load BSSeq object for Nematostella
nematostella_complete_bsseq = HDF5Array::loadHDF5SummarizedExperiment("nematostella_complete_bsseq")

# Define groups to test
dmso_f0_samples = c("DMSO_rep1", "DMSO_rep2", "DMSO_rep3")
gsk_f0_samples = c("GSK_rep1", "GSK_rep2", "GSK_rep3")
gsk_sperm_F0_samples = c("NvIT17B14_sperm", "NvIT17B19_sperm", "NvIT17C6_sperm", "NvIT16K4_sperm_Oct31")
gsk_f1_samples = c("Nv_Cr1", "Nv_Cr2")

# Find DMLs between DMSO and GSK samples. Took 8 minutes with 4 cores.
system.time({dmso_dmls = DMLtest(nematostella_complete_bsseq, smoothing = T,
  group1 = dmso_f0_samples, group2 = gsk_f0_samples, ncores = 1)})
saveRDS(dmso_dmls, "dmso_dmls.rds")

# Merge differentially metyhlyated CpGs into regions with a change greater than 0.5.
# 0.5 chosen based on hist(abs(callDMR(dmso_dmls, delta = 0)$diff.Methy))
# All DMSO DMRs gain methylation
dmso_dmrs = callDMR(dmso_dmls, delta = 0.5)
dmso_dmrs = makeGRangesFromDataFrame(dmso_dmrs, keep.extra.columns = T)
dmso_dmrs = sort(dmso_dmrs)
names(dmso_dmrs) = paste0("dmr_", seq_along(dmso_dmrs))
saveRDS(dmso_dmrs, "dmso_dmrs.rds")
dmso_dmrs = readRDS("dmso_dmrs.rds")
rtracklayer::export.bed(dmso_dmrs, "dmso_dmrs.bed")

### Make boxplots and heatmap for the recovery of DMSO DMR methylation across developmental timepoints

# Load meth RSE object for Nematostella
nematostella_complete_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("nematostella_complete_meth_rse")

# Get mean methylation for DMSO DMRs in all samples. Took 10 seconds
system.time({dmso_dmr_methylation = summarizeRegionMethylation(meth_rse = nematostella_complete_meth_rse, 
  genomic_regions = dmso_dmrs, BPPARAM = BiocParallel::SerialParam())})

# Select samples representing timepoints of interest
timepoint_samples = c(dmso_f0_samples, gsk_f0_samples, "skimming_three_mpf_15X", "skimming_five_mpf_11X", gsk_sperm_F0_samples, gsk_f1_samples)

# Create a vector with the names of timepoints and name this vector with corresponding sample names
timepoints = c(rep("Untreated\nGastrula", 3), rep("Treated F0\nGastrula", 3), "Treated F0\n≤ 3 Months", "Treated F0\n≥ 5 Months", rep("Treated F0\nSperm", 4), 
  rep("Treated F1\nGastrula", 2))
timepoints = setNames(timepoints, timepoint_samples)

# Create a long-format data.frame with DMSO DMR methylation for the selected samples
dmso_dmr_methylation_timepoints_df_long = tidyr::pivot_longer(tibble::rownames_to_column(dmso_dmr_methylation[timepoint_samples], "dmr_id"), -dmr_id)
dmso_dmr_methylation_timepoints_df_long$timepoint = timepoints[dmso_dmr_methylation_timepoints_df_long$name]

# Put timepoint levels in correct order
dmso_dmr_methylation_timepoints_df_long$timepoint = factor(dmso_dmr_methylation_timepoints_df_long$timepoint, 
  levels = unique(dmso_dmr_methylation_timepoints_df_long$timepoint))

# Create boxplot for methylation recovery in DMSO samples
dmr_meth_recovery_plot = ggplot(dmso_dmr_methylation_timepoints_df_long, aes(x = timepoint, y = value)) +
  geom_boxplot(outlier.shape = NA, fill = "#D2C465")
dmr_meth_recovery_plot = customize_ggplot_theme(dmr_meth_recovery_plot, title = "Recovery of Methylation in DMRs", 
  xlab = "\nDevelopmental Timepoint", ylab = "Mean Methylation")
dmr_meth_recovery_plot
ggsave(plot = dmr_meth_recovery_plot, "plots/dmso_dmr_meth_recovery_boxplots.pdf", width = 16, height = 9)

### Create a heatmap of methylation of DMRs across timepoints

# Subset dmso_dmr_methylation_timepoint_samples for timepoint samples (except DMSO)
dmso_dmr_methylation_timepoint_samples = dmso_dmr_methylation[, timepoint_samples]

# Merge the samples from the same timepoint by taking their mean DMR methylation
merged_dmso_dmr_methylation_timepoint_samples = data.frame(
  "Treated F0 Gastrula" = rowMeans(dmso_dmr_methylation_timepoint_samples[gsk_f0_samples], na.rm = T),
  "Treated F0 3 Months" = rowMeans(dmso_dmr_methylation_timepoint_samples["skimming_three_mpf_15X"], na.rm = T),
  "Treated F0 5 Months" = rowMeans(dmso_dmr_methylation_timepoint_samples["skimming_five_mpf_11X"], na.rm = T),
  "Treated F0 Sperm" = rowMeans(dmso_dmr_methylation_timepoint_samples[gsk_sperm_F0_samples], na.rm = T),
  "Treated F1 Gastrula" = rowMeans(dmso_dmr_methylation_timepoint_samples[gsk_f1_samples], na.rm = T)
)

# Impute missing values
set.seed(123)
merged_dmso_dmr_methylation_timepoint_samples = impute::impute.knn(as.matrix(merged_dmso_dmr_methylation_timepoint_samples))$data

# Create column labels for the heatmap
heatmap_labels = gsub("Treated ", "", gsub("\n", " ", unique(timepoints)[-1]))

# Create a heatmap showing clustering of DMRs by methylation across timepoints into two groups using KNN 
# and then cluster these groups using hierarchical clustering with Ward's D2 method
set.seed(123)
pdf(file = "plots/dmso_dmr_meth_recovery_heatmap.pdf", width = 9, height = 9, bg = "white")
#png(filename = "plots/dmso_dmr_meth_recovery_heatmap.png", width = 16, height = 9, units = "in", res = 300)
system.time({dmr_heatmap = draw(Heatmap(matrix = as.matrix(merged_dmso_dmr_methylation_timepoint_samples), col = c("#4B878BFF", "white", "#D01C1FFF"), border = T,
  column_title = "Clustering of DMRs by Methylation Recovery", column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  clustering_method_rows = "ward.D2", row_km = 2, show_row_names = F, row_title = NULL, row_dend_width = unit(25, "mm"), gap = unit(5, "mm"), 
  cluster_columns = F, , column_labels = heatmap_labels, column_names_rot = 0, column_names_centered = T, column_names_gp = gpar(fontsize = 10), 
  name = "DMR\nMethylation", heatmap_legend_param = list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14), legend_height = unit(6, "cm"))))})
dev.off()

# Extract the DMRs from the two clusters and create separate GRanges with slow and fast recovery DMRs
dmr_clusters = row_order(dmr_heatmap)
slow_recovery_dmr_names = row.names(merged_dmso_dmr_methylation_timepoint_samples)[dmr_clusters[[1]]]
fast_recovery_dmr_names = row.names(merged_dmso_dmr_methylation_timepoint_samples)[dmr_clusters[[2]]]

# Check that DMR clusters are behaving as expected and save them
colMeans(merged_dmso_dmr_methylation_timepoint_samples[slow_recovery_dmr_names, ], na.rm = T)
colMeans(merged_dmso_dmr_methylation_timepoint_samples[fast_recovery_dmr_names, ], na.rm = T)

# Indicate if dmso_dmrs are fast or slow recovery and resave dmso_dmrs
dmso_dmrs$dmr_status = ifelse(names(dmso_dmrs) %in% slow_recovery_dmr_names, "Slow", "Fast")
saveRDS(dmso_dmrs, "dmso_dmrs.rds")

# Export BED files for fast and slow DMRs
export.bed(dmso_dmrs[dmso_dmrs$dmr_status == "Fast"], "dmso_dmrs_fast.bed")
export.bed(dmso_dmrs[dmso_dmrs$dmr_status == "Slow"], "dmso_dmrs_slow.bed")

# Create plot including DMSO samples but with the same clustering

# Merge the samples from the same timepoint by taking their mean DMR methylation with untreated samples
merged_dmso_dmr_methylation_timepoint_samples_with_untreated = data.frame(
  "Untreated Gastrula" = rowMeans(dmso_dmr_methylation_timepoint_samples[dmso_f0_samples], na.rm = T),
  "Treated F0 Gastrula" = rowMeans(dmso_dmr_methylation_timepoint_samples[gsk_f0_samples], na.rm = T),
  "Treated F0 3 Months" = rowMeans(dmso_dmr_methylation_timepoint_samples["skimming_three_mpf_15X"], na.rm = T),
  "Treated F0 5 Months" = rowMeans(dmso_dmr_methylation_timepoint_samples["skimming_five_mpf_11X"], na.rm = T),
  "Treated F0 Sperm" = rowMeans(dmso_dmr_methylation_timepoint_samples[gsk_sperm_F0_samples], na.rm = T),
  "Treated F1 Gastrula" = rowMeans(dmso_dmr_methylation_timepoint_samples[gsk_f1_samples], na.rm = T)
)

# Impute missing values
set.seed(123)
merged_dmso_dmr_methylation_timepoint_samples_with_untreated = impute::impute.knn(as.matrix(merged_dmso_dmr_methylation_timepoint_samples_with_untreated))$data

# Create column labels for the heatmap
heatmap_labels = gsub("Treated ", "", gsub("\n", " ", unique(timepoints)))

# Create a heatmap showing clustering of DMRs by methylation across timepoints into two groups using KNN 
# and then cluster these groups using hierarchical clustering with Ward's D2 method
set.seed(123)
pdf(file = "plots/dmso_dmr_meth_recovery_heatmap_with_untreated_samples.pdf", width = 9, height = 9, bg = "white")
#png(filename = "plots/dmso_dmr_meth_recovery_heatmap.png", width = 16, height = 9, units = "in", res = 300)
system.time({dmr_heatmap2 = draw(Heatmap(matrix = as.matrix(merged_dmso_dmr_methylation_timepoint_samples_with_untreated), col = c("#4B878BFF", "white", "#D01C1FFF"), border = T,
  column_title = "Clustering of DMRs by Methylation Recovery", column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  split = seq_along(row.names(merged_dmso_dmr_methylation_timepoint_samples_with_untreated)) %in% dmr_clusters[[1]], 
  show_row_names = F, row_title = NULL, row_dend_width = unit(25, "mm"), gap = unit(5, "mm"), 
  cluster_columns = F, , column_labels = heatmap_labels, column_names_rot = 0, column_names_centered = T, column_names_gp = gpar(fontsize = 10), 
  name = "DMR\nMethylation", heatmap_legend_param = list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14), legend_height = unit(6, "cm"))))})
dev.off()

### Plot the mean methylation of the different DMR groups for the different timepoints

# Get mean methylation for DMSO DMRs
dmso_dmr_timepoint_mean_methylation_df = data.frame(
  slow_dmrs = colMeans(dmso_dmr_methylation[slow_recovery_dmr_names, timepoint_samples], na.rm = T),
  fast_dmrs = colMeans(dmso_dmr_methylation[fast_recovery_dmr_names, timepoint_samples], na.rm = T),
  timepoint = unname(timepoints))

# Convert dmso_dmr_timepoint_mean_methylation_df to long format
dmso_dmr_timepoint_mean_methylation_df_long = tidyr::pivot_longer(dmso_dmr_timepoint_mean_methylation_df, -timepoint, 
  names_to = "dmr_class", values_to = "dmr_methylation")

# Convert timepoint to a factor in the right order
dmso_dmr_timepoint_mean_methylation_df_long$timepoint = 
  factor(dmso_dmr_timepoint_mean_methylation_df_long$timepoint, levels = unique(timepoints))

# Create a plot of the points and add a smoothed curve
set.seed(123)
timecourse_plot = ggplot(dmso_dmr_timepoint_mean_methylation_df_long , aes(x = timepoint, y = dmr_methylation, fill = dmr_class)) +
  geom_jitter(size = 8, alpha = 0.75, shape = 21, width = 0.05) +
  geom_smooth(aes(x = as.numeric(timepoint), y = dmr_methylation, color = dmr_class), se = F, show.legend = FALSE)
timecourse_plot = customize_ggplot_theme(timecourse_plot, xlab = "Timepoint", ylab = "Mean Methylation", 
  title = "Recovery of Methylation in DMR Groups", 
  fill_title = "DMR Group", fill_labels = c("Fast Recovery", "Slow Recovery"))
timecourse_plot
ggsave(plot = timecourse_plot, "plots/dmr_meth_recovery_timecourse_plot.pdf", width = 16, height = 9)

### Make boxplots and PCA plot of histone marks in DMRs

# Get paths to bigwig files for histone marks and extract name of marks
histone_bigwigs = list.files("../histone_marks/navarrete_bigwigs/", full.names = T)
histone_marks = gsub(".log2r.input.bw", "", gsub(".*Nvec_", "", histone_bigwigs))

# Summarize histone mark values for each DMR in dmr_gsk_vs_dmso_gr. Took 12 seconds
system.time({dmso_dmr_histone_values = bigwig_summarize_over_regions(bw_filepaths = histone_bigwigs, 
  gr = dmso_dmrs, column_names = histone_marks, statistic = "mean0")})

# Add DMR status to dmso_dmr_histone_values
dmso_dmr_histone_values$dmr_status = dmso_dmrs$dmr_status

# Test if there is a difference in medians for each histone mark between the slow and fast DMRs using Mood's median test
# mood_test_results = lapply(histone_marks, function(x)
#   data.frame(p_value = coin::pvalue(coin::median_test(dmso_dmr_histone_values[[x]] ~ factor(dmso_dmr_histone_values$dmr_status)))))
wilcoxon_test_results = lapply(histone_marks, function(x)
  data.frame(p_value = wilcox.test(dmso_dmr_histone_values[[x]] ~ factor(dmso_dmr_histone_values$dmr_status))$p.value))

# Convert results to a data.frame and correct p-values
wilcoxon_test_results = bind_rows(setNames(wilcoxon_test_results, histone_marks), .id = "histone_mark")
wilcoxon_test_results$q_value = p.adjust(wilcoxon_test_results$p_value, method = "BH")
wilcoxon_test_results$significance = plotR::sig_sym(wilcoxon_test_results$q_value)

# Convert table to long format
dmso_dmr_histone_values_long = tidyr::pivot_longer(dmso_dmr_histone_values, -dmr_status, names_to = "histone_mark")

# Make a vector matching histone marks to their effect and add effects to dmso_dmr_histone_values_long
histones_to_effects  = setNames(c("Activating Marks", "Repressive Marks", "Gene Body Marks", "Activating Marks", "Activating Marks", "Gene Body Marks", 
  "Gene Body Marks", "Gene Body Marks", "Activating Marks", "Repressive Marks", "Repressive Marks", "Activating Marks"), histone_marks)
dmso_dmr_histone_values_long$effect = histones_to_effects[dmso_dmr_histone_values_long$histone_mark]
wilcoxon_test_results$effect = histones_to_effects[wilcoxon_test_results$histone_mark]

# Make boxplots of histone marks for DMRs
histone_mark_boxplots = ggplot() +
  geom_boxplot(data = dmso_dmr_histone_values_long, mapping = aes(x = histone_mark, y = value, fill = dmr_status)) +
  geom_text(data = wilcoxon_test_results, mapping = aes(label = significance, x = histone_mark), y = 3, size = 10)
histone_mark_boxplots = customize_ggplot_theme(histone_mark_boxplots, title = "Histone Marks at DMRs", 
  xlab = "Histone Mark", ylab = "Value", x_labels_angle = 45, facet = "effect", facet_scales = "free_x")
histone_mark_boxplots
ggsave(plot = histone_mark_boxplots, "plots/dmr_histone_mark_boxplots.pdf", width = 16, height = 9)

# Create a PCA plot of DMR histone marks 
dmr_pca = prcomp(as.matrix(select(dmso_dmr_histone_values, -dmr_status)), , scale. = T, center = T)
dmr_pca_plot = pca_ggplot(dmr_pca, colour_groups = dmso_dmrs$dmr_status,
  alpha = 0.5, point_size = 2, title = "PCA of DMR Histone Marks", colour_legend_title = "DMR Group") + theme(legend.position = c(0.90, 0.85))
dmr_pca_plot
ggsave(plot = dmr_pca_plot, "plots/dmr_histone_pca_plot.pdf", width = 9, height = 9)

# Extract the absolute loading scores for the 2nd PC 
pc2_loadings = sort(abs(dmr_pca$rotation[, 2]), decreasing = T)

# Create a barplot of the absolute loading scores for PC2
pc2_loadings_df = data.frame(histone_mark = names(pc2_loadings), score = pc2_loadings, row.names = NULL)
pc2_loadings_df$histone_mark = factor(pc2_loadings_df$histone_mark, levels = pc2_loadings_df$histone_mark)
pc2_loadings_barplot = ggplot(pc2_loadings_df, aes(x = histone_mark, y = score)) + geom_col(fill = "#A50F15", color = "black")
pc2_loadings_barplot = customize_ggplot_theme(pc2_loadings_barplot, title = "Contribution of Histone Marks to PC2",
  xlab = "Histone Mark", ylab = "Absolute Loading Scores", x_labels_angle = 45)
pc2_loadings_barplot
ggsave(plot = pc2_loadings_barplot, "plots/dmr_histone_pc2_barplot.pdf", width = 9, height = 9)

### Plot the distribution of DMRs around TSS

# Load Nvec genes and get their TSS
nvec_gene_gr = readRDS("../nematostella_genome/nvec_pc_genes_gr.rds")
nvec_gene_tss_gr = resize(nvec_gene_gr, 1, fix = "start")

# Find bidirectional TSS defined as TSS with another TSS within 1KB on the opposite strand
plus_tss = nvec_gene_tss_gr[strand(nvec_gene_tss_gr) == "+"] 
minus_tss = nvec_gene_tss_gr[strand(nvec_gene_tss_gr) == "-"] 
plus_minus_overlaps = data.frame(findOverlaps(expand_granges(plus_tss, 500, 0), expand_granges(minus_tss, 500, 0), ignore.strand = T))
bidirectional_promoter_genes = c(plus_tss$gene[plus_minus_overlaps$queryHits], minus_tss$gene[plus_minus_overlaps$subjectHits])

# Add 20 KB upstream and downstream of each TSS
nvec_gene_tss_gr_expanded = methodical::expand_granges(nvec_gene_tss_gr, 20000, 20000)

# Find which genes overlap each DMR
dmr_gene_overlaps = data.frame(findOverlaps(dmso_dmrs, nvec_gene_tss_gr_expanded))

# Create a GRanges for the DMRs for each gene
dmrs_for_genes = dmso_dmrs[dmr_gene_overlaps$queryHits]
dmrs_for_genes$gene = names(nvec_gene_tss_gr_expanded)[dmr_gene_overlaps$subjectHits]

# Convert the DMRs to relative ranges
relative_ranges_dmrs = methodical:::rangesRelativeToTSS(dmrs_for_genes, tss_gr = nvec_gene_tss_gr[dmrs_for_genes$gene])

# Bin the relative ranges into 500 bp windows
binned_ranges_dmrs = bin_relative_ranges(relative_ranges = relative_ranges_dmrs, bin_start = -19750, 
  bin_end = 19750, bin_step = 500, category = relative_ranges_dmrs$dmr_status)
  
# Convert to long format with separate rows for negative and positive dmrs
binned_ranges_dmrs = reshape2::melt(binned_ranges_dmrs, id.vars = "bin_center")

# Plot the distribution of DMRs around TSS
dmr_distribution_plot = ggplot(binned_ranges_dmrs, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge", color = "black")
dmr_distribution_plot = customize_ggplot_theme(plot = dmr_distribution_plot, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of DMRs",
  title = NULL, fill_title = "DMR Status", axis_text_size = 14,
  fill_labels = c("Fast", "Slow"), scale_x = scale_x_continuous(breaks = seq(-20000, 20000, 2500), expand = c(0, 0), labels = scales::comma)) +
  geom_vline(xintercept = 0, linetype = "dotted")
dmr_distribution_plot
ggsave(plot = dmr_distribution_plot, "plots/dmr_distribution_around_tss_barplot.pdf", width = 16, height = 9)

### Create plot only for genes with unidirectional promoters

# Convert the DMRs to relative ranges
relative_ranges_dmrs_unidirectional = relative_ranges_dmrs[!relative_ranges_dmrs$gene %in% bidirectional_tss]

# Bin the relative ranges into 500 bp windows
binned_ranges_dmrs_unidirectional = bin_relative_ranges(relative_ranges = relative_ranges_dmrs_unidirectional, bin_start = -19750, 
  bin_end = 19750, bin_step = 500, category = relative_ranges_dmrs_unidirectional$dmr_status)
  
# Convert to long format with separate rows for negative and positive dmrs
binned_ranges_dmrs_unidirectional = reshape2::melt(binned_ranges_dmrs_unidirectional, id.vars = "bin_center")

# Plot the distribution of DMRs around TSS
unidirectional_dmr_distribution_plot = ggplot(binned_ranges_dmrs_unidirectional, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge", color = "black")
unidirectional_dmr_distribution_plot = customize_ggplot_theme(plot = unidirectional_dmr_distribution_plot, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of DMRs",
  title = NULL, fill_title = "DMR Status", axis_text_size = 14,
  fill_labels = c("Fast", "Slow"), scale_x = scale_x_continuous(breaks = seq(-20000, 20000, 2500), expand = c(0, 0), labels = scales::comma)) +
  geom_vline(xintercept = 0, linetype = "dotted")
unidirectional_dmr_distribution_plot
ggsave(plot = unidirectional_dmr_distribution_plot, "plots/unidirectional_dmr_distribution_around_tss_barplot.pdf", width = 16, height = 9)

###

# Separate gene bodies into 20 equally sized sections
gene_body_sections = tile(nvec_gene_gr, n = 20)

# Reverse the GRanges for transcripts on the - strand
gene_body_sections[strand(nvec_gene_gr) == "-"] = GRangesList(lapply(gene_body_sections[strand(nvec_gene_gr) == "-"], rev))

# Convert gene_body_sections into a flat GRanges and number the sections
gene_body_sections = unlist(unname(gene_body_sections))
gene_body_sections$gene = names(gene_body_sections)
gene_body_sections$region = paste("gene_body", 1:20, sep = "_")

# Get upstream and downstream regions for each gene
upstream_sections = tile(flank(nvec_gene_gr, width = 20000, start = T), width = 500)
upstream_sections[strand(nvec_gene_gr) == "-"] = GRangesList(lapply(upstream_sections[strand(nvec_gene_gr) == "-"], rev))
upstream_sections = unlist(unname(upstream_sections))
upstream_sections$gene = names(upstream_sections)
upstream_sections$region = paste0("TSS-", seq(20000, 500, -500))

downstream_sections = tile(flank(nvec_gene_gr, width = 20000, start = F), width = 500)
downstream_sections[strand(nvec_gene_gr) == "-"] = GRangesList(lapply(downstream_sections[strand(nvec_gene_gr) == "-"], rev))
downstream_sections = unlist(unname(downstream_sections))
downstream_sections$gene = names(downstream_sections)
downstream_sections$region = paste0("TES+", seq(500, 20000, 500))

# Combine upstream, gene body and downstream sections
combined_sections = c(upstream_sections, gene_body_sections, downstream_sections)
combined_sections$region = factor(combined_sections$region, 
  levels = c(paste0("TSS-", seq(20000, 500, -500)), paste("gene_body", 1:20, sep = "_"), paste0("TES+", seq(500, 20000, 500))))

# Find overlaps of DMRs and combined gene sections
overlaps_df = data.frame(findOverlaps(dmso_dmrs, combined_sections, ignore.strand = F))
overlaps_df$dmr_status = dmso_dmrs$dmr_status[overlaps_df$queryHits]
overlaps_df$gene = combined_sections$gene[overlaps_df$subjectHits]
overlaps_df$region = combined_sections$region[overlaps_df$subjectHits]
overlaps_df$region_class = ifelse(grepl("TSS", overlaps_df$region), "Upstream", 
  ifelse(grepl("gene", overlaps_df$region), "Gene Body", "Downstream"))
overlaps_df$region_class = factor(overlaps_df$region_class, levels = c("Upstream", "Gene Body", "Downstream"))

overlaps_summarized = summarise(group_by(overlaps_df, region_class, region, dmr_status), count = n())
regions_plot_normalized = ggplot(overlaps_summarized, aes(x = region, y = count, fill = dmr_status)) +
    geom_col(position = "dodge", color = "black")
regions_plot_normalized = customize_ggplot_theme(regions_plot_normalized, 
    title = NULL, xlab = "Region", ylab = "DMR Count")
regions_plot_normalized = regions_plot_normalized + facet_grid(~region_class, scales = "free_x", space = "free_x") + 
    theme(panel.spacing = unit(0,'lines'), strip.background = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
      strip.text = element_text(size = 12))
regions_plot_normalized
ggsave(plot = regions_plot_normalized , "plots/dmr_distributions_whole_gene_body.pdf", width = 16, height = 9)

# Convert transcripts_gr_sections into a flat GRanges
nvec_gene_gr_sections = unlist(unname(nvec_gene_gr_sections))

# Add transcript ID and the section number as metadata columns
transcripts_gr_sections$transcript_id = names(transcripts_gr_sections)
transcripts_gr_sections$region = paste("Region", 1:100)

# Convert transcripts_gr_sections back into a GRangesList
tss_grl = split(transcripts_gr_sections, transcripts_gr_sections$transcript_id)

# Load introns and exons for transcripts as a GRanges list
#tss_grl = readRDS("../auxillary_data/pc_transcripts_exons_and_introns_grl.rds")
tss_grl_expanded = expand_transcripts(grl = tss_grl, expand_upstream = seq(500, 5000, 500), expand_downstream = seq(500, 5000, 500))

### Create boxplots of expression of genes associated with DMRs in DMSO samples

# Separate fast and slow DMRs
fast_dmrs = dmso_dmrs[dmso_dmrs$dmr_status == "Fast"]
slow_dmrs = dmso_dmrs[dmso_dmrs$dmr_status == "Slow"]

# Get GRanges for genes
nvec_genes_gr = readRDS("../nematostella_genome/nvec_pc_genes_gr.rds")
nvec_genes_tss_gr = promoters(nvec_genes_gr, 0, 1)

# Find which genes only overlap fast DMRs, which genes only overlap slow DMRs, which genes overlap both and which don't overlap any DMRs
fast_dmr_genes = subsetByOverlaps(nvec_genes_gr, fast_dmrs)$gene
slow_dmr_genes = subsetByOverlaps(nvec_genes_gr, slow_dmrs)$gene
dual_dmr_genes = intersect(fast_dmr_genes, slow_dmr_genes)
fast_dmr_genes = setdiff(fast_dmr_genes, dual_dmr_genes)
slow_dmr_genes = setdiff(slow_dmr_genes, dual_dmr_genes)
non_dmr_genes = subsetByOverlaps(nvec_genes_gr, dmso_dmrs, invert = T)$gene

# Turn into a list
dmr_genes_list = list(
  "Dual DMRs" = dual_dmr_genes,
  "Fast DMRs" = fast_dmr_genes,
  "Slow DMRs" = slow_dmr_genes,
  "No DMRs" = non_dmr_genes
)

# Find for what proportion of dual DMR genes, the fast DMR is closest to the TSS
dual_dmrs = subsetByOverlaps(dmso_dmrs, nvec_genes_gr[nvec_genes_gr$gene %in% dual_dmr_genes])
dual_dmr_gene_overlaps = data.frame(findOverlaps(dual_dmrs, nvec_genes_gr))
dual_dmr_gene_overlaps$dmr_type = dual_dmrs$dmr_status[dual_dmr_gene_overlaps$queryHits]
dual_dmr_gene_overlaps$gene = nvec_genes_gr$gene[dual_dmr_gene_overlaps$subjectHits]
dual_dmr_gene_overlaps$distance = distance(dual_dmrs[dual_dmr_gene_overlaps$queryHits], nvec_genes_tss_gr[dual_dmr_gene_overlaps$subjectHits])
dual_dmr_gene_overlaps = filter(dual_dmr_gene_overlaps, gene %in% dual_dmr_genes)
dual_dmr_gene_overlaps = group_by(dual_dmr_gene_overlaps, gene)
dual_dmr_gene_overlaps_nearest = filter(dual_dmr_gene_overlaps, distance == min(distance))
prop.table(table(dual_dmr_gene_overlaps_nearest$dmr_type))

# Turn dmr_genes_list inside out and convert into a data.frame
genes_to_dmr_list = unlist(with(reshape2::melt(dmr_genes_list), split(L1, value)))
dmr_genes_df = tidyr::unnest(tibble::enframe(dmr_genes_list, value = "gene"), "gene")

# Load table with gene expression counts and convert into TPM
nvec_gene_counts = data.frame(data.table::fread("../rnaseq/nvec_complete_gene_counts_raw.tsv.gz"), row.names = 1)
nvec_gene_tpm = data.frame(apply(nvec_gene_counts, 2, function(x) x/sum(x)*1e6))

# Get mean TPM values in DMSO and GSK samples
mean_gene_tpm_dmso = rowMeans(select(nvec_gene_tpm, starts_with("DMSO")))
mean_gene_tpm_gsk = rowMeans(select(nvec_gene_tpm, starts_with("GSK")))
dmr_genes_df$dmso_value = mean_gene_tpm_dmso[dmr_genes_df$gene]
dmr_genes_df$gsk_value = mean_gene_tpm_gsk[dmr_genes_df$gene]
dmr_genes_df$name = factor(dmr_genes_df$name, levels = c("No DMRs", "Slow DMRs", "Fast DMRs", "Dual DMRs"))

# Create a boxplot of the expression of genes which overlap no DMRs, only the DMSO DMRs or the fast DMRs
dmr_expression_boxplots = ggplot(dmr_genes_df, aes(x = name, y = log10(dmso_value), fill = name)) +
  geom_boxplot()
dmr_expression_boxplots = customize_ggplot_theme(dmr_expression_boxplots, 
  title = "Expression of Genes Associated with DMRs", ylab = "Log10 TPM",
  fill_colors = c("darkgrey", "#00BFC4", "#F8766D", "#BA66B0"), 
  show_legend = F) + ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")), size = 10,
    comparisons = list(c("No DMRs", "Slow DMRs"), c("Slow DMRs", "Fast DMRs"), c("Fast DMRs", "Dual DMRs")))
dmr_expression_boxplots
ggsave(plot = dmr_expression_boxplots, "plots/dmr_associated_gene_expression_boxplots_with_dual_dmrs.pdf", width = 16, height = 9)

dmr_expression_boxplots_gsk_samples = ggplot(dmr_genes_df, aes(x = name, y = log10(gsk_value), fill = name)) +
  geom_boxplot()
dmr_expression_boxplots_gsk_samples = customize_ggplot_theme(dmr_expression_boxplots_gsk_samples, 
  title = "Expression of Genes Associated with DMRs in GSK-Treated Samples", ylab = "Log10 TPM",
  fill_colors = c("darkgrey", "#00BFC4", "#F8766D", "#BA66B0"), 
  show_legend = F) + ggpubr::stat_compare_means(method = "wilcox.test", 
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")), size = 10,
    comparisons = list(c("No DMRs", "Slow DMRs"), c("Slow DMRs", "Fast DMRs"), c("Fast DMRs", "Dual DMRs")))
dmr_expression_boxplots_gsk_samples
ggsave(plot = dmr_expression_boxplots_gsk_samples, "plots/dmr_associated_gene_expression_boxplots_with_dual_dmrs_gsk_samples.pdf", width = 16, height = 9)

### Create boxplots of expression of genes associated with DMRs across developmental timepoints

# Get the mean expression of genes at each timepoint
combined_timepoint_raw_counts = data.frame(data.table::fread("../rnaseq/combined_timepoint_gene_counts.tsv.gz"), row.names = 1)

# Normalize combined_timepoint_raw_counts across samples using DESeq2
combined_timepoint_raw_counts_dds = DESeqDataSetFromMatrix(countData = combined_timepoint_raw_counts, 
  colData = data.frame(sample = names(combined_timepoint_raw_counts)), design = ~ 1)
combined_timepoint_raw_counts_dds  = estimateSizeFactors(combined_timepoint_raw_counts_dds) 
combined_timepoint_deseq_counts = data.frame(counts(combined_timepoint_raw_counts_dds, normalized = T))

# Normalize across genes using max normalization
combined_timepoint_deseq_counts_max_normalized = data.frame(t(apply(combined_timepoint_deseq_counts, 1, function(x) x/max(x))))

# Convert combined_timepoint_deseq_counts_max_normalized to long format
combined_timepoint_deseq_counts_max_normalized_long = 
  tidyr::pivot_longer(tibble::rownames_to_column(combined_timepoint_deseq_counts_max_normalized, "gene_id"), -gene_id, names_to = "timepoint")

# Indicate which DMRs genes are associated with
combined_timepoint_deseq_counts_max_normalized_long$dmr_class = genes_to_dmr_list[combined_timepoint_deseq_counts_max_normalized_long$gene_id]

# Remove genes with missing values
combined_timepoint_deseq_counts_max_normalized_long = 
  combined_timepoint_deseq_counts_max_normalized_long[complete.cases(combined_timepoint_deseq_counts_max_normalized_long), ]

# Put timepoint and DMR class in desired order  
combined_timepoint_deseq_counts_max_normalized_long$timepoint = factor(combined_timepoint_deseq_counts_max_normalized_long$timepoint, 
  levels = names(combined_timepoint_deseq_counts_max_normalized))
combined_timepoint_deseq_counts_max_normalized_long$dmr_class = 
  factor(combined_timepoint_deseq_counts_max_normalized_long$dmr_class, levels = c("No DMRs", "Slow DMRs", "Fast DMRs", "Dual DMRs"))

# Create boxplots of DMR gene expression across timepoints
timepoint_boxplots = ggplot(combined_timepoint_deseq_counts_max_normalized_long, aes(x = timepoint, y = value, fill = dmr_class)) + 
  geom_boxplot()
timepoint_boxplots = customize_ggplot_theme(timepoint_boxplots, facet = "dmr_class", facet_nrow = 4, 
  x_labels = gsub("Hpf_", "", levels(combined_timepoint_deseq_counts_max_normalized_long$timepoint)),
  fill_colors = c("darkgrey", "#00BFC4", "#F8766D", "#BA66B0"), 
  title = "Expression of Genes Associated with DMRs Across Development",
    xlab = "Hours Post-Fertilization", ylab = "Max Normalized Gene Expression") +
  theme(strip.background = element_rect(fill = "white"))
timepoint_boxplots
ggsave(plot = timepoint_boxplots, "plots/dmr_timepoint_boxplots_dual_dmrs.pdf", width = 16, height = 9)

# Calculate TPM values from combined_timepoint_raw_counts
combined_timepoint_raw_counts_tpm = data.frame(apply(combined_timepoint_raw_counts, 2, function(x) x/sum(x)*1e6))

# Find which genes are expressed before 24 Hpf in DMSO samples using median TPM of 10 as a cutoff
expressed_before_24hpf = names(which(matrixStats::rowMedians(as.matrix(combined_timepoint_raw_counts_tpm[1:3])) > 10))

# Use a chi-squared test to test if there is a significant difference in the proportion of fast_dmr_genes and 
# slow_dmr_genes expressed before and after 24 Hpf. There is a highly significant difference (0.83 vs 0.59; p-value < 2.2e-16)
prop.test(x = c(
  length(intersect(fast_dmr_genes, expressed_before_24hpf)),
  length(intersect(slow_dmr_genes, expressed_before_24hpf))), 
  n = c(length(fast_dmr_genes), length(slow_dmr_genes)), 
  alternative = "greater")



### Old stuff

# Trim complete_annotation_gr for regions in gene bodies and turn the results into a list
complete_annotation_gr = trim_overlaps(complete_annotation_gr, gene_body_gr)
complete_annotation_grl = split(complete_annotation_gr, complete_annotation_gr$class)
complete_annotation_grl$pc_gene_body = NULL
complete_annotation_grl$introns = gene_body_gr[gene_body_gr$subclass == "intron"]
complete_annotation_grl$exon = gene_body_gr[gene_body_gr$subclass == "exon"]
complete_annotation_grl$all_repeats = trim_overlaps(nvec_repeats_gr, gene_body_gr)

# Get slow and fast DMRs inside gene bodies
slow_dmrs_gene_bodies = trim_overlaps(dmso_dmrs[dmso_dmrs$dmr_status == "Slow"], gene_body_gr)
fast_dmrs_gene_bodies = trim_overlaps(dmso_dmrs[dmso_dmrs$dmr_status == "Fast"], gene_body_gr)

# Count proportion of DMSO vs GSK DMRs overlapping each class of genomic region
slow_dmr_region_overlap_proportion = sort(sapply(complete_annotation_grl, function(x) 
  genomicTools::calculate_regions_intersections(slow_dmrs_gene_bodies, x, overlap_measure = "proportion")))
fast_dmr_region_overlap_proportion = sort(sapply(complete_annotation_grl, function(x) 
  genomicTools::calculate_regions_intersections(fast_dmrs_gene_bodies, x, overlap_measure = "proportion")))
background_overlap_proportion = sort(sapply(complete_annotation_grl, function(x) 
  genomicTools::calculate_regions_intersections(gene_body_gr, x, overlap_measure = "proportion")))

# Combine results in a data.frame
x = bind_rows(background = stack(background_overlap_proportion), slow = stack(slow_dmr_region_overlap_proportion), fast = stack(fast_dmr_region_overlap_proportion), .id = "tmr_status")
ggplot(x, aes(y = ind, x = values, fill = tmr_status)) +
  geom_col(position = "dodge")

dmr_region_overlap_proportion = data.frame(
  region = names(dmr_region_overlap_proportion),
  proportion = dmr_region_overlap_proportion, row.names = NULL)
dmr_region_overlap_proportion = filter(arrange(dmr_region_overlap_proportion, proportion), region != "ARTEFACT")
dmr_region_overlap_proportion$region = factor(dmr_region_overlap_proportion$region, levels = dmr_region_overlap_proportion$region)

# Count proportion of DMSO vs GSK DMRs overlapping each class of genomic region
dmr_region_overlap_proportion = sort(sapply(complete_annotation_class_grl, function(x) 
  length(subsetByOverlaps(dmso_dmrs, x, ignore.strand = T))/length(dmr_gsk_vs_dmso_gr)))
dmr_region_overlap_proportion = data.frame(
  region = names(dmr_region_overlap_proportion),
  proportion = dmr_region_overlap_proportion, row.names = NULL)
dmr_region_overlap_proportion = filter(arrange(dmr_region_overlap_proportion, proportion), region != "ARTEFACT")
dmr_region_overlap_proportion$region = factor(dmr_region_overlap_proportion$region, levels = dmr_region_overlap_proportion$region)

# Find overlaps between DMRs and complete_annotation_gr
dmr_annotation_overlaps = data.frame(findOverlaps(dmso_dmrs, complete_annotation_gr, ignore.strand = T))
dmr_annotation_overlaps$dmr_status = dmso_dmrs$dmr_status[dmr_annotation_overlaps$queryHits]
dmr_annotation_overlaps$annotation_class = complete_annotation_gr$class[dmr_annotation_overlaps$subjectHits]
dmr_annotation_overlaps$annotation_subclass = complete_annotation_gr$subclass[dmr_annotation_overlaps$subjectHits]

# Group 
summarize(group_by(dmr_annotation_overlaps, dmr_status, annotation_class), count = n())

# Count proportion of DMSO vs GSK DMRs overlapping each class of genomic region
dmr_region_overlap_proportion = sort(sapply(complete_annotation_class_grl, function(x) 
  length(subsetByOverlaps(dmso_dmrs, x, ignore.strand = T))/length(dmr_gsk_vs_dmso_gr)))
dmr_region_overlap_proportion = data.frame(
  region = names(dmr_region_overlap_proportion),
  proportion = dmr_region_overlap_proportion, row.names = NULL)
dmr_region_overlap_proportion = filter(arrange(dmr_region_overlap_proportion, proportion), region != "ARTEFACT")
dmr_region_overlap_proportion$region = factor(dmr_region_overlap_proportion$region, levels = dmr_region_overlap_proportion$region)

# Create x labels for the plot
x_labels = as.character(dmr_region_overlap_proportion$region)
x_labels[c(4, 6, 7, 11, 14, 15)] = c("Pseudogene", "Simple Repeat", "lncRNA", "PC Promoter", "DNA Transposon", "PC Gene Body")

# Create barplot for number of DMRs overlapping each class of region
dmr_region_overlap_proportion_plot = ggplot(dmr_region_overlap_proportion, aes(x = region, y = proportion, fill = region)) +
  geom_col(color="black", fill = colorRampPalette(RColorBrewer::brewer.pal(n = 20, name = "RdBu"))(15)) +
  scale_fill_brewer(palette="Set1")
dmr_region_overlap_proportion_plot = customize_ggplot_theme(dmr_region_overlap_proportion_plot, 
  xlab = "Region", ylab = "Proportion", x_labels = x_labels,
  title = "Proportion of DMSO vs GSK DMRs Overlapping\nDifferent Classes of Genomic Region") +
  coord_flip()
dmr_region_overlap_proportion_plot
ggsave(plot = dmr_region_overlap_proportion_plot, "plots/dmso_dmr_region_overlap_proportion_barplot.pdf", width = 9, height = 9)


###

##### Sort this stuff

# Find DMLs between sperm and GSK samples. Took 20 minutes with 1 core.
system.time({sperm_dmls = DMLtest(nematostella_complete_bsseq, smoothing = T,
  group1 = gsk_sperm_F0_samples, group2 = gsk_f0_samples, ncores = 1)})

# Merge differentially metyhlyated CpGs into regions with an absolute change grater than 0.5.
# All sperm DMRs gain methylation
sperm_dmrs = callDMR(sperm_dmls, delta = 0.5)
sperm_dmrs = makeGRangesFromDataFrame(sperm_dmrs, keep.extra.columns = T)
names(sperm_dmrs) = paste0("dmr_", seq_along(sperm_dmrs))
saveRDS(sperm_dmrs, "sperm_dmrs.rds")
sperm_dmrs = readRDS("sperm_dmrs.rds")
rtracklayer::export.bed(sperm_dmrs, "sperm_dmrs.bed")

# Find DMLs between egg and GSK samples. Took 20 minutes with 1 core.
system.time({egg_dmls = DMLtest(nematostella_complete_bsseq, smoothing = T,
  group1 = gsk_egg_F0_samples, group2 = gsk_f0_samples, ncores = 1)})

# Merge differentially metyhlyated CpGs into regions with an absolute change grater than 0.5.
# No egg DMRs found
egg_dmrs = callDMR(egg_dmls, delta = 0.5)

# Find DMLs between F1 and GSK samples. Took 20 minutes with 1 core.
system.time({f1_dmls = DMLtest(nematostella_complete_bsseq, smoothing = T,
  group1 = gsk_f1_samples, group2 = gsk_f0_samples, ncores = 1)})

# Merge differentially metyhlyated CpGs into regions with an absolute change grater than 0.5.
# All f1 DMRs gain methylation
f1_cross_dmrs = callDMR(f1_dmls, delta = 0.5)
f1_cross_dmrs = makeGRangesFromDataFrame(f1_cross_dmrs, keep.extra.columns = T)
names(f1_cross_dmrs) = paste0("dmr_", seq_along(f1_cross_dmrs))
saveRDS(f1_cross_dmrs, "f1_cross_dmrs.rds")
f1_cross_dmrs = readRDS("f1_cross_dmrs.rds")

# Check overlap of dmso_dmrs, sperm_dmrs and f1_cross_dmrs. 
# 88% of sperm DMRs overlap DMSO DMRs and 95% of F1 cross DMRs overlap sperm DMRs
sum(width(reduce(intersect(sperm_dmrs, dmso_dmrs), ignore.strand = T)))/
  sum(width(reduce(sperm_dmrs, ignore.strand = T)))
sum(width(reduce(intersect(f1_cross_dmrs, sperm_dmrs), ignore.strand = T)))/
  sum(width(reduce(f1_cross_dmrs, ignore.strand = T)))

# Create a Venn Diagram for overlap of DMSO and sperm DMRs
pdf("plots/dmso_and_sperm_dmr_venn.pdf", width = 16, height = 9)
plot(
  eulerr::euler(c(
    DMSO = length(subsetByOverlaps(dmso_dmrs, sperm_dmrs, invert = T)), 
    Sperm = length(subsetByOverlaps(sperm_dmrs, dmso_dmrs, invert = T)), 
    "DMSO&Sperm" = length(subsetByOverlaps(dmso_dmrs, sperm_dmrs)))
    ), fills = c("#7B5C90", "#3288BD"), labels = T, quantities = T)
dev.off()

# Find the regions in common to DMSO and sperm DMRs. 88% of the base pairs in sperm_dmrs overlap dmso_dmrs
common_dmrs = intersect(dmso_dmrs, sperm_dmrs)
common_dmr_genes = unique(subsetByOverlaps(nvec_gene_gr, common_dmrs)$gene)

# Find the regions in DMSO DMRs which do not overlap sperm DMRs
unique_dmso_dmrs = subsetByOverlaps(dmso_dmrs, sperm_dmrs, invert = T)
unique_dmso_dmr_genes = unique(subsetByOverlaps(nvec_gene_gr, unique_dmso_dmrs)$gene)

### Find which classes of regions DMSO DMRs overlap

# Load Granges with genomic annotation of Nematostella and split into a GRangesList based on class of region
complete_annotation_gr = readRDS("~/genomes/nematostella/repeats/nvec_complete_annotation_gr.rds")
complete_annotation_class_grl = split(complete_annotation_gr, complete_annotation_gr$class)

# Count proportion of DMSO vs GSK DMRs overlapping each class of genomic region
dmr_region_overlap_proportion = sort(sapply(complete_annotation_class_grl, function(x) 
  length(subsetByOverlaps(dmso_dmrs, x, ignore.strand = T))/length(dmr_gsk_vs_dmso_gr)))
dmr_region_overlap_proportion = data.frame(
  region = names(dmr_region_overlap_proportion),
  proportion = dmr_region_overlap_proportion, row.names = NULL)
dmr_region_overlap_proportion = filter(arrange(dmr_region_overlap_proportion, proportion), region != "ARTEFACT")
dmr_region_overlap_proportion$region = factor(dmr_region_overlap_proportion$region, levels = dmr_region_overlap_proportion$region)

# Create x labels for the plot
x_labels = as.character(dmr_region_overlap_proportion$region)
x_labels[c(4, 6, 7, 11, 14, 15)] = c("Pseudogene", "Simple Repeat", "lncRNA", "PC Promoter", "DNA Transposon", "PC Gene Body")

# Create barplot for number of DMRs overlapping each class of region
dmr_region_overlap_proportion_plot = ggplot(dmr_region_overlap_proportion, aes(x = region, y = proportion, fill = region)) +
  geom_col(color="black", fill = colorRampPalette(RColorBrewer::brewer.pal(n = 20, name = "RdBu"))(15)) +
  scale_fill_brewer(palette="Set1")
dmr_region_overlap_proportion_plot = customize_ggplot_theme(dmr_region_overlap_proportion_plot, 
  xlab = "Region", ylab = "Proportion", x_labels = x_labels,
  title = "Proportion of DMSO vs GSK DMRs Overlapping\nDifferent Classes of Genomic Region") +
  coord_flip()
dmr_region_overlap_proportion_plot
ggsave(plot = dmr_region_overlap_proportion_plot, "plots/dmso_dmr_region_overlap_proportion_barplot.pdf", width = 9, height = 9)

### Plot Methylation of DMRs in DMSO and GSK samples

# Add a column indicating if DMRs overlap sperm DMRs
dmso_dmr_methylation$overlap_sperm_dmrs = ifelse(dmso_dmrs %over% sperm_dmrs, "Overlap Sperm DMRs", "Unique DMSO DMRs")

# Convert dmso_dmr_methylation to long format 
dmso_dmr_methylation_long = tidyr::pivot_longer(dmso_dmr_methylation, cols = -overlap_sperm_dmrs, 
  values_to = "methylation", names_to = "sample")

# Filter for DMSO, GSK and sperm samples
dmso_dmr_methylation_long = filter(dmso_dmr_methylation_long, sample %in% 
    c(dmso_f0_samples, gsk_f0_samples, gsk_sperm_F0_samples, gsk_egg_F0_samples, gsk_f1_samples))

# Indicate what sample group samples are from
dmso_dmr_methylation_long$sample_group = case_when(
  grepl("DMSO", dmso_dmr_methylation_long$sample) ~ "F0 DMSO",
  grepl("GSK", dmso_dmr_methylation_long$sample) ~ "F0 GSK",
  grepl("sperm", dmso_dmr_methylation_long$sample) ~ "F0 Sperm",
  grepl("egg", dmso_dmr_methylation_long$sample) ~ "F0 Egg",
  grepl("_Cr", dmso_dmr_methylation_long$sample) ~ "F1 GSK Cross"
)

# Put groups in desired order
dmso_dmr_methylation_long$sample_group = factor(dmso_dmr_methylation_long$sample_group, 
  levels = c("F0 DMSO", "F0 GSK", "F0 Sperm", "F0 Egg", "F1 GSK Cross"))

# Make boxplots of methylation of DMRs in different sample groups
dmr_methylation_boxplots = ggplot(dmso_dmr_methylation_long, aes(x = sample_group, y = methylation, fill = overlap_sperm_dmrs)) +
  geom_boxplot()
dmr_methylation_boxplots = customize_ggplot_theme(dmr_methylation_boxplots, 
  xlab = "Sample Group", ylab = "Mean Methylation", fill_colors = c("#3288BD", "#7B5C90"))
dmr_methylation_boxplots
ggsave(plot = dmr_methylation_boxplots, "plots/dmr_methylation_boxplots.pdf", width = 16, height = 9)

### Make boxplots of expression of genes overlapping DMRs

# Get GRanges for genes
nvec_genes_gr = readRDS("../nematostella_genome/nvec_pc_genes_gr.rds")

# Find which genes overlap dmr_gsk_vs_dmso_gr and which overlap dmr_gsk_vs_dmso_gr_overlapping_dmr_F0_sperm_vs_gsk_gr
sperm_dmr_genes = subsetByOverlaps(nvec_genes_gr, sperm_dmrs)$gene
unique_dmso_dmr_genes = setdiff(subsetByOverlaps(nvec_genes_gr, dmso_dmrs)$gene, sperm_dmr_genes)
non_dmr_genes = subsetByOverlaps(nvec_genes_gr, dmso_dmrs, invert = T)$gene

# Turn into a list
dmr_genes_list = list(
  "Sperm DMRs" = sperm_dmr_genes,
  "DMSO DMRs Only" = unique_dmso_dmr_genes,
  "No DMRs" = non_dmr_genes
)

# Turn dmr_genes_list inside out and convert into a data.frame
genes_to_dmr_list = unlist(with(reshape2::melt(dmr_genes_list), split(L1, value)))
dmr_genes_df = tidyr::unnest(tibble::enframe(dmr_genes_list, value = "gene"), "gene")

# Load table with gene expression counts and convert into TPM
nvec_gene_counts = data.frame(data.table::fread("../rnaseq/nvec_complete_gene_counts_raw.tsv.gz"), row.names = 1)
nvec_gene_tpm = data.frame(apply(nvec_gene_counts, 2, function(x) x/sum(x)*1e6))

# Get mean TPM values in DMSO samples
mean_gene_tpm_dmso = rowMeans(select(nvec_gene_tpm, starts_with("DMSO")))
dmr_genes_df$value = mean_gene_tpm_dmso[dmr_genes_df$gene]
dmr_genes_df$name = factor(dmr_genes_df$name, levels = c("No DMRs", "DMSO DMRs Only", "Sperm DMRs"))

# Create a boxplot of the expression of genes which overlap no DMRs, only the DMSO DMRs or the sperm DMRs
dmr_expression_boxplots = ggplot(dmr_genes_df, aes(x = name, y = log10(value), fill = name)) +
  geom_boxplot()
dmr_expression_boxplots = customize_ggplot_theme(dmr_expression_boxplots, 
  title = "Log10 TPM of Genes Associated with DMRs", xlab = "Group", ylab = "Log10 TPM",
  fill_colors = colour_list$three_blues, show_legend = F)
dmr_expression_boxplots
ggsave(plot = dmr_expression_boxplots, "plots/dmr_associated_gene_expression_boxplots.pdf", width = 16, height = 9)

### Make boxplots of the expression of genes overlapping DMRs at different developmental timepoints

# Get the mean expression of genes at each timepoint
combined_timepoint_raw_counts = data.frame(data.table::fread("../rnaseq/combined_timepoint_gene_counts.tsv.gz"), row.names = 1)

# Normalize combined_timepoint_raw_counts across samples using DESeq2
combined_timepoint_raw_counts_dds = DESeqDataSetFromMatrix(countData = combined_timepoint_raw_counts, 
  colData = data.frame(sample = names(combined_timepoint_raw_counts)), design = ~ 1)
combined_timepoint_raw_counts_dds  = estimateSizeFactors(combined_timepoint_raw_counts_dds) 
combined_timepoint_deseq_counts = data.frame(counts(combined_timepoint_raw_counts_dds, normalized = T))

# Normalize across genes using max normalization
combined_timepoint_deseq_counts_max_normalized = data.frame(t(apply(combined_timepoint_deseq_counts, 1, function(x) x/max(x))))

# Convert combined_timepoint_deseq_counts_max_normalized to long format
combined_timepoint_deseq_counts_max_normalized_long = 
  tidyr::pivot_longer(tibble::rownames_to_column(combined_timepoint_deseq_counts_max_normalized, "gene_id"), -gene_id, names_to = "timepoint")

# Indicate which DMRs genes are associated with
combined_timepoint_deseq_counts_max_normalized_long$dmr_class = genes_to_dmr_list[combined_timepoint_deseq_counts_max_normalized_long$gene_id]

# Remove genes with missing values
combined_timepoint_deseq_counts_max_normalized_long = 
  combined_timepoint_deseq_counts_max_normalized_long[complete.cases(combined_timepoint_deseq_counts_max_normalized_long), ]

# Put timepoint and DMR class in desired order  
combined_timepoint_deseq_counts_max_normalized_long$timepoint = factor(combined_timepoint_deseq_counts_max_normalized_long$timepoint, 
  levels = names(combined_timepoint_deseq_counts_max_normalized))
combined_timepoint_deseq_counts_max_normalized_long$dmr_class = 
  factor(combined_timepoint_deseq_counts_max_normalized_long$dmr_class, levels = c("No DMRs", "DMSO DMRs Only", "Sperm DMRs"))

# Create boxplots of DMR gene expression across timepoints
timepoint_boxplots = ggplot(combined_timepoint_deseq_counts_max_normalized_long, aes(x = timepoint, y = value, fill = dmr_class)) + 
  geom_boxplot()
timepoint_boxplots = customize_ggplot_theme(timepoint_boxplots, facet = "dmr_class", facet_nrow = 3, 
  fill_colors = c("#bfab25", "#7B5C90", "#3288BD"), xlab = "Timepoint", ylab = "Max Normalized Gene Expression")
timepoint_boxplots
ggsave(plot = timepoint_boxplots, "plots/dmr_timepoint_boxplots.pdf", width = 16, height = 9)

# Calculate TPM values from combined_timepoint_raw_counts
combined_timepoint_raw_counts_tpm = data.frame(apply(combined_timepoint_raw_counts, 2, function(x) x/sum(x)*1e6))

# Find which genes are expressed before 24 Hpf in DMSO samples using median TPM of 10 as a cutoff
expressed_before_24hpf = names(which(matrixStats::rowMedians(as.matrix(combined_timepoint_raw_counts_tpm[1:3])) > 10))

# Use a chi-squared test to test if there is a significant difference in the proportion of sperm_dmr_genes and 
# unique_dmso_dmr_genes expressed before and after 24 Hpf. There is a highly significant difference (0.83 vs 0.59; p-value < 2.2e-16)
prop.test(x = c(
  length(intersect(sperm_dmr_genes, expressed_before_24hpf)),
  length(intersect(unique_dmso_dmr_genes, expressed_before_24hpf))), 
  n = c(length(sperm_dmr_genes), length(unique_dmso_dmr_genes)), 
  alternative = "greater")

### Make boxplots of histone marks in DMRs

# Get paths to bigwig files for histone marks and extract name of marks
histone_bigwigs = list.files("../histone_marks/bigwigs/", full.names = T)
histone_marks = gsub(".log2r.input.bw", "", gsub(".*Nvec_", "", histone_bigwigs))

# Sumamrize histone mark values for each DMR in dmr_gsk_vs_dmso_gr
system.time({dmso_dmr_histone_values = genomicTools::bigwig_summarize_over_regions(bw_filepaths = histone_bigwigs, 
  gr = dmso_dmrs, column_names = histone_marks, statistic = "mean0")})
data.table::fwrite(dmso_dmr_histone_values, "dmso_dmr_histone_values.tsv.gz", sep = "\t")
dmso_dmr_histone_values = data.frame(data.table::fread("dmso_dmr_histone_values.tsv.gz"), row.names = 1)

# # Add a column indicating if DMRs overlap sperm DMRs
# dmso_dmr_histone_values$overlap_sperm_dmrs = ifelse(dmso_dmrs %over% sperm_dmrs, "fast_recovery", "slow_recovery")
# 
# # Convert table to long format
# dmso_dmr_histone_values_long = tidyr::pivot_longer(dmso_dmr_histone_values, -overlap_sperm_dmrs)

# Convert table to long format
dmso_dmr_histone_values_long = tidyr::pivot_longer(tibble::rownames_to_column(dmr_gsk_vs_dmso_histone_values, "genomic_region"), -genomic_region)

# Indicate if DMRs overlapped dmr_F0_sperm_vs_gsk_gr
dmso_dmr_histone_values_long$overlap_sperm_dmrs = ifelse(dmso_dmr_histone_values_long$genomic_region %in% 
    names(dmr_gsk_vs_dmso_gr_overlapping_dmr_F0_sperm_vs_gsk_gr), "Overlapping Sperm/GSK DMRs", "Not Overlapping Sperm/GSK DMRs")

# Make boxplots of histone marks for DMRs
histone_mark_boxplots = ggplot(dmso_dmr_histone_values_long, aes(x = name, y = value, fill = overlap_sperm_dmrs)) +
  geom_boxplot()
histone_mark_boxplots = customize_ggplot_theme(histone_mark_boxplots, title = "Histone Marks at DMRs", 
  xlab = "Histone Mark", ylab = "Value", x_labels_angle = 45, fill_colors = c("#7B5C90", "#3288BD"))
histone_mark_boxplots
ggsave(plot = histone_mark_boxplots, "histone_mark_boxplots.pdf", width = 16, height = 9)

# Create a PCA plot 
dmr_pca = prcomp(as.matrix(dmr_gsk_vs_dmso_histone_values), , scale. = T, center = T)
dmr_pca_plot = pca_ggplot(dmr_pca, colour_groups = ifelse(row.names(dmr_gsk_vs_dmso_histone_values) %in% 
    names(dmr_gsk_vs_dmso_gr_overlapping_dmr_F0_sperm_vs_gsk_gr), "Overlapping Sperm/GSK DMRs", "Not Overlapping Sperm/GSK DMRs"),
  alpha = 0.5, size = 2, legend_key_size = 1) + theme(legend.position = c(0.85, 0.2))
dmr_pca_plot
ggsave(plot = dmr_pca_plot, "dmr_pca_plot.pdf", width = 16, height = 16)

### Old

### Extract methylation values for unique_dmr_gsk_vs_dmso_gr

# Load meth RSE object for Nematostella
system.time({nematostella_complete_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("nematostella_complete_meth_rse/")})
bpparam = BiocParallel::SerialParam()

# Get mean methylation for unique_dmr_gsk_vs_dmso_gr in all samples. Took 10 seconds
system.time({dmr_gsk_vs_dmso_meth = summarizeRegionMethylation(meth_rse = nematostella_complete_meth_rse, 
  genomic_regions = dmr_gsk_vs_dmso_gr, BPPARAM = bpparam)})

pheatmap::pheatmap(dmr_gsk_vs_dmso_meth[complete.cases(dmr_gsk_vs_dmso_meth), grep("egg", names(dmr_gsk_vs_dmso_meth), invert = T)], 
  kmeans_k = 200, show_rownames = F, labels_col = grep("egg", colData(nematostella_complete_meth_rse)$description, invert = T, value = T), 
  filename = "dmr_heatmap_all_samples.pdf", width = 16, height = 9)

# Get mean methylation for unique_dmr_gsk_vs_dmso_gr in all samples. Took 10 seconds
system.time({unique_dmr_meth = summarizeRegionMethylation(meth_rse = nematostella_complete_meth_rse, 
  genomic_regions = unique_dmr_gsk_vs_dmso_gr, BPPARAM = bpparam)})
unique_dmr_meth_long = tidyr::pivot_longer(unique_dmr_meth, tidyselect::everything())

# Get mean methylation for unique_dmr_F0_sperm_vs_gsk_gr in all samples. Took 10 seconds
system.time({sperm_dmr_meth = summarizeRegionMethylation(meth_rse = nematostella_complete_meth_rse, 
  genomic_regions = dmr_F0_sperm_vs_gsk_gr, BPPARAM = bpparam)})
sperm_dmr_meth_long = tidyr::pivot_longer(sperm_dmr_meth, tidyselect::everything())

# Get mean methylation for unique_dmr_F0_sperm_vs_gsk_gr in all samples. Took 10 seconds
system.time({unique_sperm_dmr_meth = summarizeRegionMethylation(meth_rse = nematostella_complete_meth_rse, 
  genomic_regions = unique_dmr_F0_sperm_vs_gsk_gr, BPPARAM = bpparam)})
unique_sperm_dmr_meth_long = tidyr::pivot_longer(unique_sperm_dmr_meth, tidyselect::everything())

combined_dmr_df = bind_rows(list(gsk_vs_dmso_unique = unique_dmr_meth_long, sperm_vs_gsk = sperm_dmr_meth_long, unique_sperm = unique_sperm_dmr_meth_long), .id = "dmr_class")
combined_dmr_df$group = setNames(colData(nematostella_complete_meth_rse)$group, colnames(nematostella_complete_meth_rse))[combined_dmr_df$name]
combined_dmr_df$name = setNames(colData(nematostella_complete_meth_rse)$description, colnames(nematostella_complete_meth_rse))[combined_dmr_df$name]
combined_dmr_df$group = factor(combined_dmr_df$group, 
  levels = c("Control F0", "Control F0 Sperm", "Control F0 Egg", "GSK F0", "GSK F0 Sperm", "GSK F0 Egg", 
  "Morpholino F0 Sperm", "Control x Control F1", "GSK x GSK F1", "Morpholino x Morpholino F1", "Control x GSK F1"))

colData(nematostella_complete_meth_rse)$group = c(rep("Control F0", 3), rep("GSK F0", 3), "GSK F0 Egg", 
  rep("GSK F0 Sperm", 2), "GSK F0 Egg", "GSK F0 Sperm", "Morpholino F0 Sperm", rep("Control F0 Sperm", 2), 
  rep("Morpholino x Morpholino F1", 2), rep("GSK x GSK F1", 2), rep("Control F0 Egg", 2), "GSK F0 Sperm",
  "Control x Control F1", rep("Control x GSK F1", 2))

boxplots = ggplot(combined_dmr_df, aes(x = name, y = value, fill = dmr_class)) +
  geom_boxplot(outlier.shape = NA)
boxplots = customize_ggplot_theme(boxplots, ylab = "DMR Methylation", x_labels_angle = 45, xlab = "Sample", 
  fill_title = "DMR Class", fill_colors = c("#7B5C90", "#3288BD", "#bfab25"),
  fill_labels = c("Unique to GSK VS DMSO", "Shared with Sperm Vs GSK", "Unique Sperm Vs GSK"),
  facet = "group", facet_nrow = 1, facet_scales = "free_x", strip_text_size = 10) + 
  theme(strip.background = element_blank())
boxplots
ggsave(plot = boxplots, "dmr_methylation_boxplots.pdf", width = 24, height = 13.5)

### Make boxplots of histone marks in DMRs

# Get paths to bigwig files for histone marks and extract name of marks
histone_bigwigs = list.files("../histone_marks/bigwigs/", full.names = T)
histone_marks = gsub(".log2r.input.bw", "", gsub(".*Nvec_", "", histone_bigwigs))

# Sumamrize histone mark values for each DMR in dmr_gsk_vs_dmso_gr
system.time({dmr_gsk_vs_dmso_histone_values = genomicTools::bigwig_summarize_over_regions(bw_filepaths = histone_bigwigs, 
  gr = dmr_gsk_vs_dmso_gr, column_names = histone_marks, statistic = "mean0")})
data.table::fwrite(dmr_gsk_vs_dmso_histone_values, "dmr_gsk_vs_dmso_histone_values_test.tsv.gz", sep = "\t")
dmr_gsk_vs_dmso_histone_values = data.frame(data.table::fread("dmr_gsk_vs_dmso_histone_values_test.tsv.gz"), row.names = 1)
dmr_gsk_vs_dmso_histone_values = tibble::column_to_rownames(dmr_gsk_vs_dmso_histone_values, "genomic_region")

# Convert table to long format
dmr_gsk_vs_dmso_histone_values_long = tidyr::pivot_longer(tibble::rownames_to_column(dmr_gsk_vs_dmso_histone_values, "genomic_region"), -genomic_region)

# Indicate if DMRs overlapped dmr_F0_sperm_vs_gsk_gr
dmr_gsk_vs_dmso_histone_values_long$dmr_overlap = ifelse(dmr_gsk_vs_dmso_histone_values_long$genomic_region %in% 
    names(dmr_gsk_vs_dmso_gr_overlapping_dmr_F0_sperm_vs_gsk_gr), "Fast Recovery", "Slow Recovery")

dmr_gsk_vs_dmso_histone_values_long$effect = histones_to_effects[dmr_gsk_vs_dmso_histone_values_long$name]

# Make boxplots of histone marks for DMRs
histone_mark_boxplots = ggplot(dmr_gsk_vs_dmso_histone_values_long, aes(x = name, y = value, fill = dmr_overlap)) +
  geom_boxplot()
histone_mark_boxplots = customize_ggplot_theme(histone_mark_boxplots, title = "Histone Marks at DMRs", 
  xlab = "Histone Mark", ylab = "Value", fill_title = "DMR Group", x_labels_angle = 45, 
  facet = "effect", facet_scales = "free_x")
histone_mark_boxplots
ggsave(plot = histone_mark_boxplots, "plots/histone_mark_boxplots_test.pdf", width = 16, height = 9)

# Create a PCA plot 
dmr_pca = prcomp(as.matrix(dmr_gsk_vs_dmso_histone_values), , scale. = T, center = T)
dmr_pca_plot = pca_ggplot(dmr_pca, colour_groups = ifelse(row.names(dmr_gsk_vs_dmso_histone_values) %in% 
    names(dmr_gsk_vs_dmso_gr_overlapping_dmr_F0_sperm_vs_gsk_gr), "Overlapping Sperm/GSK DMRs", "Not Overlapping Sperm/GSK DMRs"),
  alpha = 0.5, size = 2, legend_key_size = 1) + theme(legend.position = c(0.85, 0.2))
dmr_pca_plot
ggsave(plot = dmr_pca_plot, "dmr_pca_plot.pdf", width = 16, height = 16)




### DELETE

### Compare GSK vs DMSO samples

# Find differentially methylated CpGs. Took 16 minutes with 20 cores
system.time({dml_test_gsk_vs_dmso = DMLtest(nematostella_complete_bsseq, smoothing = T,
  group1 = gsk_f0_samples, group2 = dmso_f0_samples, ncores = 20)})

### Rerun this
system.time({dml_test_dmso_vs_gsk = DMLtest(nematostella_complete_bsseq, smoothing = T,
  group1 = dmso_f0_samples, group2 = gsk_f0_samples, ncores = 4)})

# Merge differentially metyhlyated CpGs into regions with a change grater than 0.5
dml_test_gsk_vs_dmso_0.5 = callDMR(dml_test_gsk_vs_dmso, delta = 0.5)

# Convert into a GRanges, sort, name regions and save
dmr_gsk_vs_dmso_gr = makeGRangesFromDataFrame(dml_test_gsk_vs_dmso_0.5, keep.extra.columns = T)
dmr_gsk_vs_dmso_gr = sort(dmr_gsk_vs_dmso_gr)
names(dmr_gsk_vs_dmso_gr) = paste0("dmr_", seq_along(dmr_gsk_vs_dmso_gr))
saveRDS(dmr_gsk_vs_dmso_gr, "dmr_gsk_vs_dmso_gr.rds")
dmr_gsk_vs_dmso_gr = readRDS("dmr_gsk_vs_dmso_gr.rds")

### Compare recovered sperm vs DMSO samples

# Find differentially methylated CpGs. Took 16 minutes with 20 cores
system.time({dml_test_F0_sperm_vs_gsk = DMLtest(nematostella_complete_bsseq, smoothing = T,
  group1 = gsk_sperm_F0_samples, group2 = gsk_f0_samples, ncores = 20)})

# Merge differentially metyhlyated CpGs into regions with a change grater than 0.5
dml_test_F0_sperm_vs_gsk_0.5 = callDMR(dml_test_F0_sperm_vs_gsk, delta = 0.5)

# Convert into a GRanges, sort, name regions and save
dmr_F0_sperm_vs_gsk_gr = makeGRangesFromDataFrame(dml_test_F0_sperm_vs_gsk_0.5, keep.extra.columns = T)
dmr_F0_sperm_vs_gsk_gr = sort(dmr_F0_sperm_vs_gsk_gr)
names(dmr_F0_sperm_vs_gsk_gr) = paste0("dmr_", seq_along(dmr_F0_sperm_vs_gsk_gr))
saveRDS(dmr_F0_sperm_vs_gsk_gr, "dmr_F0_sperm_vs_gsk_gr.rds")
dmr_F0_sperm_vs_gsk_gr = readRDS("dmr_F0_sperm_vs_gsk_gr.rds")

# Indicate if DMRs in dmr_gsk_vs_dmso_gr overlap those in dmr_F0_sperm_vs_gsk_gr or not
dmr_gsk_vs_dmso_gr_overlapping_dmr_F0_sperm_vs_gsk_gr = subsetByOverlaps(dmr_gsk_vs_dmso_gr, dmr_F0_sperm_vs_gsk_gr)

# Find DMRs in GSK vs DMSO comparison which do not overlap DMRs from sperm vs GSK comparison
unique_dmr_gsk_vs_dmso_gr = subsetByOverlaps(dmr_gsk_vs_dmso_gr, dmr_F0_sperm_vs_gsk_gr, invert = T)
unique_dmr_F0_sperm_vs_gsk_gr = subsetByOverlaps(dmr_F0_sperm_vs_gsk_gr, dmr_gsk_vs_dmso_gr, invert = T)

# Create a Venn Diagram for overlap of DMRs
length(subsetByOverlaps(dmr_F0_sperm_vs_gsk_gr, dmr_gsk_vs_dmso_gr, minoverlap = 100))
library(VennDiagram)
pdf("dmr_venn.pdf", width = 9, height = 9)
venn = venneuler(c("GSK vs DMSO" = length(dmr_gsk_vs_dmso_gr), "Sperm Vs GSK" = length(dmr_F0_sperm_vs_gsk_gr), 
  "GSK vs DMSO&Sperm Vs GSK" = length(subsetByOverlaps(dmr_F0_sperm_vs_gsk_gr, dmr_gsk_vs_dmso_gr, ignore.strand = T)))) 
venn$labels = rev(c("Sperm Vs GSK\nn = 8,020", "GSK Vs DMSO\nn = 21,470"))
plot(venn)
dev.off()

sperm_dmrs_genes = unique(subsetByOverlaps(nvec_gene_gr, sperm_dmrs)$gene)
dmso_dmr_genes_only = setdiff(unique(subsetByOverlaps(nvec_gene_gr, dmso_dmrs)$gene), sperm_dmrs_genes)

# Turn into a list
dmr_genes_list = list(
  "Fast Recovery" = sperm_dmrs_genes,
  "Slow Recovery" = dmso_dmr_genes_only
  #no_dmrs = nvec_genes_not_overlapping_dmr_gsk_vs_dmso
)

# Turn dmr_genes_list inside out
genes_to_dmr_list = unlist(with(reshape2::melt(dmr_genes_list), split(L1, value)))

dmr_genes_df = tidyr::unnest(tibble::enframe(dmr_genes_list, value = "gene"), "gene")

# Load table with gene expression counts
nvec_gene_counts = data.frame(data.table::fread("../rnaseq/nvec_complete_gene_counts_raw.tsv.gz"), row.names = 1)
nvec_gene_tpm = data.frame(apply(nvec_gene_counts, 2, function(x) x/sum(x)*1e6))

# Get mean DMSO values
mean_gene_tpm_dmso = rowMeans(select(nvec_gene_tpm, starts_with("DMSO")))
dmr_genes_df$value = mean_gene_tpm_dmso[dmr_genes_df$gene]
dmr_genes_df$name = factor(dmr_genes_df$name, levels = c("gsk_dmso_dmrs_only", "sperm_vs_gsk_dmrs"))

# 
dmr_expression_boxplots = ggplot(dmr_genes_df, aes(x = name, y = log10(value), fill = name)) +
  geom_boxplot()
dmr_expression_boxplots = customize_ggplot_theme(dmr_expression_boxplots, xlab = "Group", ylab = "Log10 TPM",
  fill_colors = colour_list$three_blues)
dmr_expression_boxplots
ggsave(plot = dmr_expression_boxplots, "dmr_expression_boxplots.pdf", width = 16, height = 9)

### Add the timepoint with max expression for each gene
combined_timepoint_raw_counts = data.frame(data.table::fread("../rnaseq/combined_timepoint_gene_counts.tsv.gz"), row.names = 1)

# Convert to TPM
combined_timepoint_raw_counts_tpm = data.frame(apply(combined_timepoint_raw_counts, 2, function(x) x/sum(x)*1e6))

# Find which genes are expressed before 24 Hpf
expressed_before_24hpf = names(which(apply(combined_timepoint_raw_counts_tpm[1:3] > 2, 1, any)))
fisher_test_vectors(nvec_genes_overlapping_dmr_F0_sperm_vs_gsk, expressed_before_24hpf, row.names(combined_timepoint_raw_counts_tpm))
fisher_test_vectors(nvec_genes_overlapping_dmr_gsk_vs_dmso_gr_only, expressed_before_24hpf, row.names(combined_timepoint_raw_counts_tpm))
fisher_test_vectors(nvec_genes_not_overlapping_dmr_gsk_vs_dmso, expressed_before_24hpf, row.names(combined_timepoint_raw_counts_tpm))
chi_test = prop.test(x = c(
  length(intersect(nvec_genes_overlapping_dmr_F0_sperm_vs_gsk, expressed_before_24hpf)),
  length(intersect(nvec_genes_overlapping_dmr_gsk_vs_dmso_gr_only, expressed_before_24hpf))), 
  n = c(length(nvec_genes_overlapping_dmr_F0_sperm_vs_gsk), length(nvec_genes_overlapping_dmr_gsk_vs_dmso_gr_only)), 
  alternative = "greater")

# Normalize combined_timepoint_raw_counts using DESeq2
combined_timepoint_raw_counts_dds = DESeqDataSetFromMatrix(countData = combined_timepoint_raw_counts, 
  colData = data.frame(sample = names(combined_timepoint_raw_counts)), design = ~ 1)
combined_timepoint_raw_counts_dds  = estimateSizeFactors(combined_timepoint_raw_counts_dds) 
combined_timepoint_deseq_counts = data.frame(counts(combined_timepoint_raw_counts_dds, normalized = T))

# Max normalize each gene
combined_timepoint_deseq_counts_max_normalized = data.frame(t(apply(combined_timepoint_deseq_counts, 1, function(x) x/max(x))))

# Convert to long format
combined_timepoint_deseq_counts_max_normalized_long = 
  tidyr::pivot_longer(tibble::rownames_to_column(combined_timepoint_deseq_counts_max_normalized, "gene_id"), -gene_id)

# Indicate which DMRs genes are associated with
combined_timepoint_deseq_counts_max_normalized_long$dmr_class = genes_to_dmr_list[combined_timepoint_deseq_counts_max_normalized_long$gene_id]
combined_timepoint_deseq_counts_max_normalized_long = combined_timepoint_deseq_counts_max_normalized_long[complete.cases(combined_timepoint_deseq_counts_max_normalized_long), ]
combined_timepoint_deseq_counts_max_normalized_long$name = factor(combined_timepoint_deseq_counts_max_normalized_long$name, 
  levels = names(combined_timepoint_deseq_counts_max_normalized))

timepoint_boxplots = ggplot(combined_timepoint_deseq_counts_max_normalized_long, aes(x = name, y = value, fill = dmr_class)) + 
  geom_boxplot()
timepoint_boxplots = customize_ggplot_theme(timepoint_boxplots, title = "Expression of Genes Associated with DMRs Across Development",
  xlab = "Hours Post-Fertilization", ylab = "Max Normalized Gene Expression",  fill_title = "DMR Group") + 
  scale_x_discrete(labels = function(x) gsub("Hpf_", "", x))
timepoint_boxplots
ggsave(plot = timepoint_boxplots, "plots/dmr_timepoint_boxplots_test.pdf", width = 16, height = 9)

max_timepoint_per_gene = readRDS("../rnaseq/max_timepoint_per_gene.rds")
dmr_genes_df$max_timepoint = max_timepoint_per_gene[dmr_genes_df$gene]

timepoint_plot = ggplot(dmr_genes_df, aes(x = max_timepoint)) + 
  geom_bar() + facet_wrap("~ name", scales = "free_y", nrow = 3)
timepoint_plot
ggsave(plot = timepoint_plot, "timepoint_plot.pdf", width = 16, height = 9)

### Find DMRs between individual gsk_sperm_F0_samples/F1 samples and DMSO samples

# Find DMLs between DMSO and GSK samples. Took 8 minutes with 4 cores.
system.time({NvIT16K4_sperm_Oct31_dmls = DMLtest(nematostella_complete_bsseq, smoothing = T,
  group1 = "NvIT16K4_sperm_Oct31", group2 = dmso_f0_samples, ncores = 1)})

# Merge differentially metyhlyated CpGs into regions with a change greater than 0.1.
NvIT16K4_sperm_Oct31_dmrs = callDMR(NvIT16K4_sperm_Oct31_dmls, delta = 0.1)
NvIT16K4_sperm_Oct31_dmrs = makeGRangesFromDataFrame(NvIT16K4_sperm_Oct31_dmrs, keep.extra.columns = T)
NvIT16K4_sperm_Oct31_dmrs = sort(NvIT16K4_sperm_Oct31_dmrs)
names(NvIT16K4_sperm_Oct31_dmrs) = paste0("dmr_", seq_along(NvIT16K4_sperm_Oct31_dmrs))
saveRDS(NvIT16K4_sperm_Oct31_dmrs, "NvIT16K4_sperm_Oct31_dmrs.rds")
NvIT16K4_sperm_Oct31_dmrs = readRDS("NvIT16K4_sperm_Oct31_dmrs.rds")

# Find DMLs between DMSO and GSK samples. Took 8 minutes with 4 cores.
system.time({NvIT17B14_sperm_dmls = DMLtest(nematostella_complete_bsseq, smoothing = T,
  group1 = "NvIT17B14_sperm", group2 = dmso_f0_samples, ncores = 1)})

# Merge differentially metyhlyated CpGs into regions with a change greater than 0.1.
NvIT17B14_sperm_dmrs = callDMR(NvIT17B14_sperm_dmls, delta = 0.1)
NvIT17B14_sperm_dmrs = makeGRangesFromDataFrame(NvIT17B14_sperm_dmrs, keep.extra.columns = T)
NvIT17B14_sperm_dmrs = sort(NvIT17B14_sperm_dmrs)
names(NvIT17B14_sperm_dmrs) = paste0("dmr_", seq_along(NvIT17B14_sperm_dmrs))
saveRDS(NvIT17B14_sperm_dmrs, "NvIT17B14_sperm_dmrs.rds")
NvIT17B14_sperm_dmrs = readRDS("NvIT17B14_sperm_dmrs.rds")

# Find DMLs between DMSO and GSK samples. Took 8 minutes with 4 cores.
system.time({NvIT17B19_sperm_dmls = DMLtest(nematostella_complete_bsseq, smoothing = T,
  group1 = "NvIT17B19_sperm", group2 = dmso_f0_samples, ncores = 1)})

# Merge differentially metyhlyated CpGs into regions with a change greater than 0.1.
NvIT17B19_sperm_dmrs = callDMR(NvIT17B19_sperm_dmls, delta = 0.1)
NvIT17B19_sperm_dmrs = makeGRangesFromDataFrame(NvIT17B19_sperm_dmrs, keep.extra.columns = T)
NvIT17B19_sperm_dmrs = sort(NvIT17B19_sperm_dmrs)
names(NvIT17B19_sperm_dmrs) = paste0("dmr_", seq_along(NvIT17B19_sperm_dmrs))
saveRDS(NvIT17B19_sperm_dmrs, "NvIT17B19_sperm_dmrs.rds")
NvIT17B19_sperm_dmrs = readRDS("NvIT17B19_sperm_dmrs.rds")

# Find DMLs between DMSO and GSK samples. Took 8 minutes with 4 cores.
system.time({NvIT17C6_sperm_dmls = DMLtest(nematostella_complete_bsseq, smoothing = T,
  group1 = "NvIT17C6_sperm", group2 = dmso_f0_samples, ncores = 1)})

# Merge differentially metyhlyated CpGs into regions with a change greater than 0.1.
NvIT17C6_sperm_dmrs = callDMR(NvIT17C6_sperm_dmls, delta = 0.1)
NvIT17C6_sperm_dmrs = makeGRangesFromDataFrame(NvIT17C6_sperm_dmrs, keep.extra.columns = T)
NvIT17C6_sperm_dmrs = sort(NvIT17C6_sperm_dmrs)
names(NvIT17C6_sperm_dmrs) = paste0("dmr_", seq_along(NvIT17C6_sperm_dmrs))
saveRDS(NvIT17C6_sperm_dmrs, "NvIT17C6_sperm_dmrs.rds")
NvIT17C6_sperm_dmrs = readRDS("NvIT17C6_sperm_dmrs.rds")

# Combine all sperm DMRs together
all_sperm_dmrs = c(NvIT16K4_sperm_Oct31_dmrs, NvIT17B14_sperm_dmrs, NvIT17B19_sperm_dmrs, NvIT17C6_sperm_dmrs)

# Filter for those which gain methylation, with at least 10 CpGs and of at least 400 bp
all_sperm_dmrs_increase = all_sperm_dmrs[all_sperm_dmrs$diff.Methy > 0 & all_sperm_dmrs$nCG >=10 & width(all_sperm_dmrs) >= 400]

all_sperm_dmrs_increase_cov = summarizeRegionMethylation(nematostella_complete_meth_rse, assay = 2, 
  genomic_regions = all_sperm_dmrs_increase, BPPARAM = bpparam)
all_sperm_dmrs_increase = all_sperm_dmrs_increase[all_sperm_dmrs_increase_cov[["NvIT16K4_sperm_Oct31"]] > 10]
all_sperm_dmrs_increase_meth = summarizeRegionMethylation(nematostella_complete_meth_rse, assay = 1, 
  genomic_regions = all_sperm_dmrs_increase, BPPARAM = bpparam)

###

NvIT16K4_sperm_Oct31_dmrs_increase = NvIT16K4_sperm_Oct31_dmrs[NvIT16K4_sperm_Oct31_dmrs$diff.Methy > 0]
NvIT16K4_sperm_Oct31_dmrs_increase = NvIT16K4_sperm_Oct31_dmrs_increase[width(NvIT16K4_sperm_Oct31_dmrs_increase) > 400]

bpparam = BiocParallel::SerialParam()
NvIT16K4_sperm_Oct31_dmrs_increase_cov = summarizeRegionMethylation(nematostella_complete_meth_rse, assay = 2, 
  genomic_regions = NvIT16K4_sperm_Oct31_dmrs_increase, BPPARAM = bpparam)
NvIT16K4_sperm_Oct31_dmrs_increase = NvIT16K4_sperm_Oct31_dmrs_increase[NvIT16K4_sperm_Oct31_dmrs_increase_cov[["NvIT16K4_sperm_Oct31"]] > 10]
NvIT16K4_sperm_Oct31_dmrs_increase_meth = summarizeRegionMethylation(nematostella_complete_meth_rse, assay = 1, 
  genomic_regions = NvIT16K4_sperm_Oct31_dmrs_increase, BPPARAM = bpparam)

# Impute missing values
set.seed(123)
NvIT16K4_sperm_Oct31_dmrs_increase_meth_imputed = impute::impute.knn(as.matrix(NvIT16K4_sperm_Oct31_dmrs_increase_meth[, timepoint_samples]))$data

NvIT16K4_sperm_Oct31_dmrs_increase = NvIT16K4_sperm_Oct31_dmrs_increase[names(which(rowMeans(NvIT16K4_sperm_Oct31_dmrs_increase_meth[, dmso_f0_samples], na.rm = T) < 0.1))]
NvIT16K4_sperm_Oct31_dmrs_increase = NvIT16K4_sperm_Oct31_dmrs_increase[NvIT16K4_sperm_Oct31_dmrs_increase$nCG >= 10]

plotR::heatmap_without_clustering(as.matrix(NvIT16K4_sperm_Oct31_dmrs_increase_meth[, timepoint_samples]), mask_diagonal = F)

Heatmap(matrix = as.matrix(NvIT16K4_sperm_Oct31_dmrs_increase_meth_imputed[, ]), col = c("#4B878BFF", "white", "#D01C1FFF"), border = T,
  column_title = "Clustering of Population DMRs", column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  clustering_method_rows = "ward.D2", row_km = 0, show_row_names = F, row_title = NULL, row_dend_width = unit(25, "mm"), gap = unit(5, "mm"), 
  cluster_columns = F, column_names_rot = 45, column_names_centered = T, column_names_gp = gpar(fontsize = 10), 
  name = "DMR\nMethylation", heatmap_legend_param = list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14), legend_height = unit(6, "cm")))

subsetByOverlaps(NvIT16K4_sperm_Oct31_dmrs_increase, population_dmrs_reduced)

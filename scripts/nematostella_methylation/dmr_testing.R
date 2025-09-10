# Test differential methylation between different conditions

# Load required packages
library(DSS)
library(methodical)
library(ComplexHeatmap)
library(dplyr)
library(DESeq2)
source("../auxillary_scripts/plotting_functions.R")
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
wilcoxon_test_results$significance = sig_sym(wilcoxon_test_results$q_value)

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
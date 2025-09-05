# Identify DMRs between different populations of Nematostella

# Load required packages
library(dplyr)
library(methrix)
library(methodical)
library(BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19)
library(DSS)
library(ComplexHeatmap)
source("../auxillary_scripts/bigwig_summarize_over_regions.R")

### Create BSSeq and Meth RSE objects from CGmap files

# Create lists of FL, MA and Vienna sample files
fl_samples = list.files("cgmap_files", pattern = "FL", full.names = T)
ma_samples = list.files("cgmap_files", pattern = "MA", full.names = T)
vienna_samples = list.files("cgmap_files", pattern = "Vienna", full.names = T)

# Combine samples into a single list
combined_samples = c(fl_samples, ma_samples, vienna_samples)

# Filter CGmap files for NC contigs
filter_nc = function(cgmap_file){
  nc_file = gsub(".CGmap_corrected.gz", "_nc_contigs.CGmap_corrected.gz", cgmap_file)
  system(paste("zcat", cgmap_file, "| grep NC_ | gzip >", nc_file ))
}
system.time({parallelR::parallelize(object = combined_samples, ncores = 9, filter_nc)})

# Get paths to NC contig files and name them
nc_files = list.files("cgmap_files/", pattern = "nc_contigs")
names(nc_files) = gsub("_nc.*", "", gsub("mC_Nvectensis_", "", basename(nc_files)))

# Create colData for samples
sample_coldata = data.frame(
  group = gsub(".*_", "", names(nc_files)),
  timepoint = gsub("_.*", "", names(nc_files)), 
  row.names = names(nc_files)
)

# Get reference CpGs from Nematostella genome. 
source("../auxillary_scripts/extract_CPGs_modified.R")
system.time({nc_cpg_sites = extract_CPGs_modified("BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19",
  seqnames = grep("NC", seqnames(BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19), value = T))})

# Create a list of NC contigs
nc_contigs = grep("NC", seqnames(BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19), value = T)
  
# Create a methrix object for the NC contigs. 
system.time({nvec_population_methrix = read_bedgraphs(files = nc_files, 
  chr_idx = 1, start_idx = 3, beta_idx = 6, cov_idx = 8, contigs = nc_contigs, zero_based = F,
  stranded = TRUE, collapse_strands = T, ref_cpgs = nc_cpg_sites, coldata = sample_coldata, 
  h5 = T, h5_dir = "nvec_population_methrix")})
nvec_population_methrix = methrix::load_HDF5_methrix("nvec_population_methrix")

# Convert methrix to a BSSeq object. 
system.time({nvec_population_bsseq = methrix::methrix2bsseq(nvec_population_methrix)})
system.time({HDF5Array::saveHDF5SummarizedExperiment(nvec_population_bsseq, "nvec_population_bsseq")})

# Convert methrix to a meth RSE object. 
system.time({nvec_population_meth_rse = methodical::methrixToRSE(nvec_population_methrix)})
system.time({HDF5Array::saveHDF5SummarizedExperiment(nvec_population_meth_rse, "nvec_population_meth_rse")})

# Combine nvec_population_bsseq with DMSO samples from nematostella_nvec_population_bsseq
nvec_population_bsseq = HDF5Array::loadHDF5SummarizedExperiment("nvec_population_bsseq")
dmso_bsseq  = HDF5Array::loadHDF5SummarizedExperiment("../nematostella_methylation/nematostella_complete_bsseq//")
dmso_bsseq = dmso_bsseq[, grep("DMSO", colnames(dmso_bsseq))]
colData(dmso_bsseq) = DataFrame(group = c(rep("DMSO", 3), rep("GSK", 3)), timepoint = rep("gastrula", 6), row.names = colnames(dmso_bsseq))
system.time({nvec_population_bsseq = combine(nvec_population_bsseq, dmso_bsseq)})

### Idenfify DMRs

# Define groups to test DMRs
dmso_samples = c("DMSO_rep1", "DMSO_rep2", "DMSO_rep3")
fl_samples = grep("FL", colnames(nvec_population_bsseq), value = T)
vienna_samples = grep("Vienna", colnames(nvec_population_bsseq), value = T)
ma_samples = grep("MA", colnames(nvec_population_bsseq), value = T)

# Find DMLs between FL and DMSO samples. Took 20 minutes with 1 core
system.time({fl_dmls = DMLtest(nvec_population_bsseq, smoothing = T,
  group1 = fl_samples, group2 = dmso_samples, ncores = 1)})
fl_dmrs = callDMR(fl_dmls, delta = 0.5)
fl_dmrs = makeGRangesFromDataFrame(fl_dmrs, keep.extra.columns = T)
fl_dmrs = sort(fl_dmrs)
names(fl_dmrs) = paste0("dmr_", seq_along(fl_dmrs))
saveRDS(fl_dmrs, "fl_dmrs.rds")

# Find DMLs between MA and DMSO samples. Took 20 minutes with 1 core
system.time({ma_dmls = DMLtest(nvec_population_bsseq, smoothing = T,
  group1 = ma_samples, group2 = dmso_samples, ncores = 1)})
ma_dmrs = callDMR(ma_dmls, delta = 0.5)
ma_dmrs = makeGRangesFromDataFrame(ma_dmrs, keep.extra.columns = T)
ma_dmrs = sort(ma_dmrs)
names(ma_dmrs) = paste0("dmr_", seq_along(ma_dmrs))
saveRDS(ma_dmrs, "ma_dmrs.rds")

# Find DMLs between Vienna and DMSO samples. Took 20 minutes with 1 core
system.time({vienna_dmls = DMLtest(nvec_population_bsseq, smoothing = T,
  group1 = vienna_samples, group2 = dmso_samples, ncores = 1)})
vienna_dmrs = callDMR(vienna_dmls, delta = 0.5)
vienna_dmrs = makeGRangesFromDataFrame(vienna_dmrs, keep.extra.columns = T)
vienna_dmrs = sort(vienna_dmrs)
names(vienna_dmrs) = paste0("dmr_", seq_along(vienna_dmrs))
saveRDS(vienna_dmrs, "vienna_dmrs.rds")

# Load DMRs for three groups
fl_dmrs = readRDS("fl_dmrs.rds")
ma_dmrs = readRDS("ma_dmrs.rds")
vienna_dmrs = readRDS("vienna_dmrs.rds")

# Combine the three groups of DMRs and filter for DMRs of at least 400 bp
population_dmrs = c(fl_dmrs, ma_dmrs, vienna_dmrs)
population_dmrs = population_dmrs[width(population_dmrs) >= 400]

# Combine all DMRs within 1000 bp
population_dmrs_reduced = reduce(population_dmrs, min.gapwidth = 1000)
population_dmrs_reduced = saveRDS(population_dmrs_reduced, "population_dmrs_reduced.rds")
population_dmrs_reduced = readRDS("population_dmrs_reduced.rds")

# Combine meth RSEs for population and DMSO samples
nvec_population_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("nvec_population_meth_rse")
dmso_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../nematostella_methylation/nematostella_complete_meth_rse//")
dmso_meth_rse = dmso_meth_rse[, grep("DMSO", colnames(dmso_meth_rse))]
dmso_meth_rse = subsetByOverlaps(dmso_meth_rse, rowRanges(nvec_population_meth_rse))
colData(dmso_meth_rse) = DataFrame(group = rep("DMSO", 3), timepoint = rep("gastrula", 3), row.names = colnames(dmso_meth_rse))
nvec_population_meth_rse = cbind(nvec_population_meth_rse, dmso_meth_rse)

# Get methylation values for DMRs
bpparam = BiocParallel::SerialParam()
system.time({population_dmrs_reduced_methylation = 
  methodical::summarizeRegionMethylation(meth_rse = nvec_population_meth_rse, genomic_regions = population_dmrs_reduced, BPPARAM = bpparam)})

# Impute missing values
set.seed(123)
population_dmrs_reduced_methylation = impute::impute.knn(as.matrix(population_dmrs_reduced_methylation))$data

# Create a heatmap of population DMR methylation
set.seed(123)
pdf(file = "popoulation_dmr_methylation_heatmap.pdf", width = 16, height = 16, bg = "white")
system.time({dmr_heatmap = draw(Heatmap(matrix = as.matrix(population_dmrs_reduced_methylation), col = c("#4B878BFF", "white", "#D01C1FFF"), border = T,
  column_title = "Clustering of DMRs", column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  clustering_method_rows = "ward.D2", row_km = 0, show_row_names = F, row_title = NULL, row_dend_width = unit(25, "mm"), gap = unit(5, "mm"), 
  cluster_columns = T, clustering_method_columns = "ward.D2", column_dend_height = unit(25, "mm"), column_names_rot = 0, column_names_centered = T, column_names_gp = gpar(fontsize = 10), 
  name = "DMR\nMethylation", heatmap_legend_param = list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14), legend_height = unit(6, "cm"))))})
dev.off()

# Check overlap of population_dmrs_reduced with DMSO DMRs. 53% of population DMRs overlap slow DMRs
dmso_dmrs = readRDS("../nematostella_methylation/dmso_dmrs.rds")
length(subsetByOverlaps(population_dmrs_reduced, dmso_dmrs[dmso_dmrs$dmr_status == "Fast"]))/length(population_dmrs_reduced)
length(subsetByOverlaps(population_dmrs_reduced, dmso_dmrs[dmso_dmrs$dmr_status == "Slow"]))/length(population_dmrs_reduced)

# Get methylation of slow and fast DMRs across populations
system.time({dmso_dmr_methylation_populations = 
  methodical::summarizeRegionMethylation(meth_rse = nvec_population_meth_rse, genomic_regions = dmso_dmrs, BPPARAM = bpparam)})
set.seed(123)

# Impute missing values
dmso_dmr_methylation_populations = impute::impute.knn(as.matrix(dmso_dmr_methylation_populations))$data
dmso_dmr_methylation_populations = dmso_dmr_methylation_populations[intersect(row.names(dmso_dmr_methylation_populations), names(dmso_dmrs[dmso_dmrs$dmr_status == "Slow"])), ]

# Create a heatmap of DMSO DMR methylation across populations
set.seed(123)
pdf(file = "dmso_dmr_methylation_across_populations_heatmap.pdf", width = 16, height = 16, bg = "white")
system.time({slow_dmr_heatmap = draw(Heatmap(matrix = as.matrix(dmso_dmr_methylation_populations), col = c("#4B878BFF", "white", "#D01C1FFF"), border = T,
  column_title = "Clustering of Slow DMRs", column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  clustering_method_rows = "ward.D2", row_km = 2, show_row_names = F, row_title = NULL, row_dend_width = unit(25, "mm"), gap = unit(5, "mm"), 
  cluster_columns = T, clustering_method_columns = "ward.D2", column_dend_height = unit(25, "mm"), column_names_rot = 0, column_names_centered = T, column_names_gp = gpar(fontsize = 10), 
  name = "DMR\nMethylation", heatmap_legend_param = list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14), legend_height = unit(6, "cm"))))})
dev.off()

# Extract the DMRs from the two clusters for the DMRs that vary and those that are stable across populations
slow_dmr_clusters = row_order(slow_dmr_heatmap)
variable_dmr_names = row.names(dmso_dmr_methylation_populations)[slow_dmr_clusters [[1]]]
stable_dmr_names = row.names(dmso_dmr_methylation_populations)[slow_dmr_clusters [[2]]]

# Combine variable and stable DMRs
variable_dmrs = dmso_dmrs[names(dmso_dmrs) %in% variable_dmr_names]
variable_dmrs$dmr_status = "variable"
stable_dmrs = dmso_dmrs[names(dmso_dmrs) %in% stable_dmr_names]
stable_dmrs$dmr_status = "stable"
stable_and_variable_dmrs = c(stable_dmrs, variable_dmrs)
saveRDS(stable_and_variable_dmrs, "stable_and_variable_dmrs.rds")

### Make boxplots of of the histone mark values for stable and variable DMRs

# Get paths to bigwig files for histone marks and extract name of marks
histone_bigwigs = list.files("../histone_marks/navarrete_bigwigs/", full.names = T)
histone_marks = gsub(".log2r.input.bw", "", gsub(".*Nvec_", "", histone_bigwigs))

# Summarize histone mark values for each DMR in dmr_gsk_vs_dmso_gr. Took 12 seconds
system.time({slow_dmr_histone_values = bigwig_summarize_over_regions(bw_filepaths = histone_bigwigs, 
  gr = stable_and_variable_dmrs, column_names = histone_marks, statistic = "mean0")})

# Add DMR status to slow_dmr_histone_values
slow_dmr_histone_values$dmr_status = stable_and_variable_dmrs$dmr_status

# Test if there is a difference in medians for each histone mark between the slow and fast DMRs using Mood's median test
mood_test_results = lapply(histone_marks, function(x)
  data.frame(p_value = coin::pvalue(coin::median_test(slow_dmr_histone_values[[x]] ~ factor(slow_dmr_histone_values$dmr_status)))))

# Convert results to a data.frame and correct p-values
mood_test_results = bind_rows(setNames(mood_test_results, histone_marks), .id = "histone_mark")
mood_test_results$q_value = p.adjust(mood_test_results$p_value, method = "BH")
mood_test_results$significance = plotR::sig_sym(mood_test_results$q_value)

# Convert table to long format
slow_dmr_histone_values_long = tidyr::pivot_longer(slow_dmr_histone_values, -dmr_status, names_to = "histone_mark")

# Make a vector matching histone marks to their effect and add effects to slow_dmr_histone_values_long
histones_to_effects  = setNames(c("Activating Marks", "Repressive Marks", "Gene Body Marks", "Activating Marks", "Activating Marks", "Gene Body Marks", 
  "Gene Body Marks", "Gene Body Marks", "Activating Marks", "Repressive Marks", "Repressive Marks", "Activating Marks"), histone_marks)
slow_dmr_histone_values_long$effect = histones_to_effects[slow_dmr_histone_values_long$histone_mark]
mood_test_results$effect = histones_to_effects[mood_test_results$histone_mark]

# Make boxplots of histone marks for DMRs
histone_mark_boxplots = ggplot() +
  geom_boxplot(data = slow_dmr_histone_values_long, mapping = aes(x = histone_mark, y = value, fill = dmr_status)) +
  geom_text(data = mood_test_results, mapping = aes(label = significance, x = histone_mark), y = 3, size = 10)
histone_mark_boxplots = customize_ggplot_theme(histone_mark_boxplots, title = "Histone Marks at Stable Vs Variable DMRs", 
  xlab = "Histone Mark", ylab = "Value", x_labels_angle = 45, facet = "effect", facet_scales = "free_x", fill_labels = c("Stable", "Variable"))
histone_mark_boxplots
ggsave(plot = histone_mark_boxplots, "stable_vs_variable_dmr_histone_mark_boxplots.pdf", width = 16, height = 9)

### Calculate methylation of population DMRs in 

bpparam = BiocParallel::SerialParam()

nematostella_complete_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../nematostella_methylation/nematostella_complete_meth_rse")

# Get mean methylation for unique_dmr_gsk_vs_dmso_gr in all samples. Took 10 seconds
system.time({population_dmrs_meth_timepoints = summarizeRegionMethylation(meth_rse = nematostella_complete_meth_rse, 
  genomic_regions = population_dmrs_reduced, BPPARAM = bpparam)})

# Select samples representing timepoints of interest
dmso_f0_samples = c("DMSO_rep1", "DMSO_rep2", "DMSO_rep3")
gsk_f0_samples = c("GSK_rep1", "GSK_rep2", "GSK_rep3")
gsk_sperm_F0_samples = c("NvIT17B14_sperm", "NvIT17B19_sperm", "NvIT17C6_sperm", "NvIT16K4_sperm_Oct31")
gsk_f1_samples = c("Nv_Cr1", "Nv_Cr2")
timepoint_samples = c(dmso_f0_samples, gsk_f0_samples, gsk_sperm_F0_samples, gsk_f1_samples)

# Create a vector with the names of timepoints and name this vector with corresponding sample names
timepoints = c(rep("Untreated\nGastrula", 3), rep("Treated F0\nGastrula", 3), "Treated F0\n≤ 3 Months", "Treated F0\n≥ 5 Months", rep("Treated F0\nSperm", 4), 
  rep("Treated F1\nGastrula", 2))
timepoints = setNames(timepoints, timepoint_samples)

# Create a long-format data.frame with DMSO DMR methylation for the selected samples
dmso_dmr_methylation_timepoints_df_long = tidyr::pivot_longer(tibble::rownames_to_column(population_dmrs_meth_timepoints[timepoint_samples], "dmr_id"), -dmr_id)
dmso_dmr_methylation_timepoints_df_long$timepoint = timepoints[dmso_dmr_methylation_timepoints_df_long$name]

# Put timepoint levels in correct order
dmso_dmr_methylation_timepoints_df_long$timepoint = factor(dmso_dmr_methylation_timepoints_df_long$timepoint, 
  levels = unique(dmso_dmr_methylation_timepoints_df_long$timepoint))

# Create boxplot for methylation recovery in DMSO samples
dmr_meth_recovery_plot = ggplot(dmso_dmr_methylation_timepoints_df_long, aes(x = timepoint, y = value)) +
  geom_boxplot(outlier.shape = NA, fill = "#D2C465")
dmr_meth_recovery_plot = customize_ggplot_theme(dmr_meth_recovery_plot, title = "Recovery of Methylation in Population DMRs", 
  xlab = "\nDevelopmental Timepoint", ylab = "Mean Methylation")
dmr_meth_recovery_plot

set.seed(123)
population_dmrs_meth_timepoints_imputed = impute::impute.knn(as.matrix(population_dmrs_meth_timepoints))$data

# Create a heatmap of population DMR methylation
set.seed(123)
pdf(file = "popoulation_dmr_methylation_timepoints_heatmap.pdf", width = 24, height = 13.5, bg = "white")
system.time({dmr_heatmap = draw(Heatmap(matrix = as.matrix(population_dmrs_meth_timepoints_imputed[, timepoint_samples]), col = c("#4B878BFF", "white", "#D01C1FFF"), border = T,
  column_title = "Clustering of Population DMRs", column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  clustering_method_rows = "ward.D2", row_km = 2, show_row_names = F, row_title = NULL, row_dend_width = unit(25, "mm"), gap = unit(5, "mm"), 
  cluster_columns = F, column_names_rot = 0, column_names_centered = T, column_names_gp = gpar(fontsize = 10), 
  name = "DMR\nMethylation", heatmap_legend_param = list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14), legend_height = unit(6, "cm"))))})
dev.off()

# Extract the DMRs from the two clusters for the DMRs that vary and those that are stable across populations
population_clusters = row_order(dmr_heatmap)
unmethylated_dmrs = row.names(population_dmrs_meth_timepoints_imputed)[population_clusters[[1]]]
methylated_dmrs = row.names(population_dmrs_meth_timepoints_imputed)[population_clusters[[2]]]

methylation_gain_regions = names(which(population_dmrs_meth_timepoints_imputed[unmethylated_dmrs, "NvIT16K4_sperm_Oct31"] - 
  rowMeans(population_dmrs_meth_timepoints_imputed[unmethylated_dmrs, dmso_f0_samples]) > 0.1))
population_dmrs_meth_timepoints_imputed[methylation_gain_regions, "Nv_Cr2"] - 
  rowMeans(population_dmrs_meth_timepoints_imputed[methylation_gain_regions, dmso_f0_samples]) > 0.1

system.time({Heatmap(matrix = as.matrix(population_dmrs_meth_timepoints_imputed[methylation_gain_regions, timepoint_samples]), col = c("#4B878BFF", "white", "#D01C1FFF"), border = T,
  column_title = "Clustering of Population DMRs", column_title_gp = gpar(fontsize = 20, fontface = "bold"), 
  cluster_rows = F, row_km = 0, show_row_names = F, row_title = NULL, row_dend_width = unit(25, "mm"), gap = unit(5, "mm"), 
  cluster_columns = F, column_names_rot = 0, column_names_centered = T, column_names_gp = gpar(fontsize = 10), 
  name = "DMR\nMethylation", heatmap_legend_param = list(title_gp = gpar(fontsize = 16), labels_gp = gpar(fontsize = 14), legend_height = unit(6, "cm")))})

plotR::heatmap_without_clustering(as.matrix(population_dmrs_meth_timepoints_imputed[methylation_gain_regions, timepoint_samples]), mask_diagonal = F)

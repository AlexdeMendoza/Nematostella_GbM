# Get methylation of different classes of regions in untreated samples

# Load required packages
library(methodical)
library(dplyr)
library(plotR)

# Load nematostella_complete_meth_rse and filter for DMSO and GSK treated samples
nematostella_complete_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("nematostella_complete_meth_rse")
nematostella_complete_meth_rse_dmso_gsk = 
  nematostella_complete_meth_rse[, c(grep("^DMSO|^GSK", colnames(nematostella_complete_meth_rse), value = T), "NvIT17B14_sperm", "NvIT17B19_sperm")]

# Load Granges with genomic annotation of Nematostella and split into a GRangesList based on class of region
complete_annotation_gr = readRDS("~/genomes/nematostella/repeats/nvec_complete_annotation_gr.rds")
complete_annotation_class_grl = split(complete_annotation_gr, complete_annotation_gr$class)

# Calculate the mean methylation in DMSO samples for each class of region. Took 30 seconds
system.time({region_class_methylation_in_dmso_and_gsk_samples = lapply(complete_annotation_class_grl, function(x) 
  colMeans(assay(subsetByOverlaps(nematostella_complete_meth_rse_dmso_gsk, x)), na.rm = T))})

# Convert result into a data.frame
region_class_methylation_in_dmso_and_gsk_samples = data.frame(do.call(rbind, region_class_methylation_in_dmso_and_gsk_samples))

# Get mean region methylation in DMSO, GSK and sperm samples and combine results
mean_region_methylation_dmso = data.frame(
  region = row.names(region_class_methylation_in_dmso_and_gsk_samples),
  mean_methylation = rowMeans(select(region_class_methylation_in_dmso_and_gsk_samples, starts_with("DMSO"))),
  treatment = "DMSO", row.names =NULL)
mean_region_methylation_gsk = data.frame(
  region = row.names(region_class_methylation_in_dmso_and_gsk_samples),
  mean_methylation = rowMeans(select(region_class_methylation_in_dmso_and_gsk_samples, starts_with("GSK"))),
  treatment = "GSK", row.names = NULL)
mean_region_methylation_sperm = data.frame(
  region = row.names(region_class_methylation_in_dmso_and_gsk_samples),
  mean_methylation = rowMeans(select(region_class_methylation_in_dmso_and_gsk_samples, contains("sperm"))),
  treatment = "Sperm", row.names = NULL)
mean_region_methylation_combined = rbind(mean_region_methylation_dmso, mean_region_methylation_gsk, mean_region_methylation_sperm)

# Put regions in the same order as complete_annotation_gr and remove ARTEFACT region
mean_region_methylation_combined$region = factor(mean_region_methylation_combined$region, levels = levels(complete_annotation_gr$class))
mean_region_methylation_combined = filter(mean_region_methylation_combined, region != "ARTEFACT")

# Create x-axis labels for the plot
x_labels = levels(complete_annotation_gr$class)
x_labels[c(1, 2, 3, 7, 8, 12)] = c("PC Promoter", "PC Gene Body", "lncRNA", "Pseudogene", "DNA Transposon", "Simple Repeat")

# Create a barplot with the mean methylation of regions in DMSO and GSK samples
region_methylation_barplot_gsk_vs_dmso = ggplot(filter(mean_region_methylation_combined, treatment != "Sperm"), aes(x = region, y = mean_methylation, fill = treatment)) +
  geom_col(position = "dodge", color = "black")
region_methylation_barplot_gsk_vs_dmso = customize_ggplot_theme(region_methylation_barplot_gsk_vs_dmso, 
  title = "Mean Methylation of Genomic Regions", xlab = "Genomic Region", ylab = "Mean Methylation", 
  x_labels = x_labels, x_labels_angle = 45, 
  fill_colors = c("#7B5C90", "#bfab25"))
region_methylation_barplot_gsk_vs_dmso
ggsave(plot = region_methylation_barplot_gsk_vs_dmso, filename = "mean_methylation_genomic_regions_dmso_and_gsk_samples.pdf", width = 16, height = 9)

# Create a barplot with the mean methylation of regions in DMSO, GSK and sperm samples
region_methylation_barplot_all = ggplot(mean_region_methylation_combined, aes(x = region, y = mean_methylation, fill = treatment)) +
  geom_col(position = "dodge", color = "black")
region_methylation_barplot_all = customize_ggplot_theme(region_methylation_barplot_all, 
  title = "Mean Methylation of Genomic Regions", xlab = "Genomic Region", ylab = "Mean Methylation",
  x_labels = x_labels, x_labels_angle = 45, 
  fill_colors = c("#7B5C90", "#bfab25", "#3288BD"))
region_methylation_barplot_all
ggsave(plot = region_methylation_barplot_all, filename = "mean_methylation_genomic_regions_dmso_gsk_and_sperm_samples.pdf", width = 16, height = 9)

# Make a data.frame giving the proportion of recovery in sperm compared to DMSO
region_sperm_recovery_df = data.frame(
  region = mean_region_methylation_sperm$region,
  recovery = mean_region_methylation_sperm$mean_methylation/mean_region_methylation_dmso$mean_methylation,
  dmso_mean_methylation = mean_region_methylation_dmso$mean_methylation)
region_sperm_recovery_df$region = factor(region_sperm_recovery_df$region, levels = levels(complete_annotation_gr$class))
region_sperm_recovery_df = filter(region_sperm_recovery_df, region != "ARTEFACT")

# Create a barplot showing the recovery of methylation in sperm
region_sperm_recovery_barplot = ggplot(region_sperm_recovery_df, aes(x = region, y = recovery)) +
  geom_col(color = "black", fill = "#3288BD")
region_sperm_recovery_barplot = region_methylation_barplot_all = customize_ggplot_theme(region_sperm_recovery_barplot, 
  title = "Recovery of Methylation in F0 Sperm", xlab = "Genomic Region", ylab = "Proportion of Methylation Recovered", 
  x_labels = x_labels, x_labels_angle = 45)
region_sperm_recovery_barplot
ggsave(plot = region_sperm_recovery_barplot, "region_sperm_recovery_barplot.pdf", width = 16, height = 9)
  
### Create barplot for subclasses of repetitive elements

# Get GRanges for repeats and split into a GRangesList based on subclass of repeat
repeats_gr = readRDS("~/genomes/nematostella/repeats/nvec_repeats_gr.rds")
repeats_grl = split(repeats_gr, repeats_gr$subclass)

# Make a list matching repeat subclasses to class
repeat_subclass_to_class = lapply(split(repeats_gr$class, repeats_gr$subclass), function(x) unique(x)[1])

# Calculate the mean methylation in DMSO and GSK samples for each subclass of repeat. Took 3 minutes
system.time({repeat_subclass_methylation_in_dmso_and_gsk_samples = lapply(repeats_grl, function(x) 
  colMeans(assay(subsetByOverlaps(nematostella_complete_meth_rse_dmso_gsk, x)), na.rm = T))})

# Convert result into a data.frame
repeat_subclass_methylation_in_dmso_and_gsk_samples = data.frame(do.call(rbind, repeat_subclass_methylation_in_dmso_and_gsk_samples))

# Get mean region methylation in DMSO and GSK samples and combine results
mean_repeat_subclass_methylation_dmso = data.frame(
  region = row.names(repeat_subclass_methylation_in_dmso_and_gsk_samples),
  mean_methylation = rowMeans(select(repeat_subclass_methylation_in_dmso_and_gsk_samples, starts_with("DMSO"))),
  treatment = "DMSO", row.names =NULL)
mean_repeat_subclass_methylation_gsk = data.frame(
  region = row.names(repeat_subclass_methylation_in_dmso_and_gsk_samples),
  mean_methylation = rowMeans(select(repeat_subclass_methylation_in_dmso_and_gsk_samples, starts_with("GSK"))),
  treatment = "GSK", row.names = NULL)
mean_repeat_subclass_methylation_sperm = data.frame(
  region = row.names(repeat_subclass_methylation_in_dmso_and_gsk_samples),
  mean_methylation = rowMeans(select(repeat_subclass_methylation_in_dmso_and_gsk_samples, contains("sperm"))),
  treatment = "Sperm", row.names = NULL)
mean_repeat_subclass_methylation_combined = rbind(mean_repeat_subclass_methylation_dmso, mean_repeat_subclass_methylation_gsk, mean_repeat_subclass_methylation_sperm)

# Add a column indicating class of repeat
mean_repeat_subclass_methylation_combined$class = unlist(repeat_subclass_to_class[mean_repeat_subclass_methylation_combined$region])

# Find the 25 most methylated subclasses in DMSO samples and filter for these subclasses
most_methylated_subclasses = head(arrange(mean_repeat_subclass_methylation_dmso, desc(mean_methylation))$region, 25)
mean_repeat_subclass_methylation_combined = filter(mean_repeat_subclass_methylation_combined, region %in% most_methylated_subclasses)

# Create a barplot with the mean methylation of repeat subclasses in DMSO and GSK samples
repeat_subclass_methylation_barplot_gsk_vs_dmso = ggplot(filter(mean_repeat_subclass_methylation_combined, treatment != "Sperm"), aes(x = region, y = mean_methylation, fill = treatment)) +
  geom_col(position = "dodge", color = "black")
repeat_subclass_methylation_barplot_gsk_vs_dmso = customize_ggplot_theme(repeat_subclass_methylation_barplot_gsk_vs_dmso, 
  title = "25 Most Methylated Repeat Subclasses", xlab = "Repeat Subclass", ylab = "Mean Methylation", 
  fill_colors = c("#7B5C90", "#bfab25"), x_labels_angle = 45) + facet_grid(~class, scales = "free_x", space = "free_x", 
  labeller = labeller(class = c(DNA = "DNA Transposon", LINE = "LINE", LTR = "LTR", SINE = "SINE"))) + 
    theme(panel.spacing = unit(0,'lines'), strip.background = element_blank(), strip.text = element_text(size = 16))
repeat_subclass_methylation_barplot_gsk_vs_dmso
ggsave(plot = repeat_subclass_methylation_barplot_gsk_vs_dmso, filename = "mean_methylation_repeat_subclass_dmso_and_gsk_samples.pdf", width = 16, height = 9)

# Create a barplot with the mean methylation of repeat subclasses in DMSO, GSK and sperm samples
repeat_subclass_methylation_barplot_all = ggplot(mean_repeat_subclass_methylation_combined, aes(x = region, y = mean_methylation, fill = treatment)) +
  geom_col(position = "dodge", color = "black")
repeat_subclass_methylation_barplot_all = customize_ggplot_theme(repeat_subclass_methylation_barplot_all, 
  title = "25 Most Methylated Repeat Subclasses", xlab = "Repeat Subclass", ylab = "Mean Methylation", 
  fill_colors = c("#7B5C90", "#bfab25", "#3288BD"), x_labels_angle = 45) + facet_grid(~class, scales = "free_x", space = "free_x", 
  labeller = labeller(class = c(DNA = "DNA Transposon", LINE = "LINE", LTR = "LTR", SINE = "SINE"))) + 
    theme(panel.spacing = unit(0,'lines'), strip.background = element_blank(), strip.text = element_text(size = 16))
repeat_subclass_methylation_barplot_all
ggsave(plot = repeat_subclass_methylation_barplot_all, filename = "mean_methylation_repeat_subclass_dmso_gsk_and_sperm_samples.pdf", width = 16, height = 9)

# Make a data.frame giving the proportion of recovery in sperm compared to DMSO
repeat_sperm_recovery_df = data.frame(
  region = mean_repeat_subclass_methylation_sperm$region,
  recovery = mean_repeat_subclass_methylation_sperm$mean_methylation/mean_repeat_subclass_methylation_dmso$mean_methylation,
  dmso_mean_methylation = mean_repeat_subclass_methylation_dmso$mean_methylation)
repeat_sperm_recovery_df = filter(repeat_sperm_recovery_df, region != "ARTEFACT")

# Add a column indicating class of repeat
repeat_sperm_recovery_df$class = unlist(repeat_subclass_to_class[repeat_sperm_recovery_df$region])

# Filter for 25 repeats with highest methylation in DMSO
repeat_sperm_recovery_df_subset = head(arrange(repeat_sperm_recovery_df, desc(dmso_mean_methylation)), 25)

# Create a barplot showing the recovery of methylation in sperm
repeat_sperm_recovery_barplot = ggplot(repeat_sperm_recovery_df_subset, aes(x = region, y = recovery)) +
  geom_col(color = "black", fill = "#3288BD")
repeat_sperm_recovery_barplot = customize_ggplot_theme(repeat_sperm_recovery_barplot, title = "Recovery of Methylation in F0 Sperm", 
  xlab = "Repeat Subclass", ylab = "Proportion of Methylation Recovered", x_labels_angle = 45) + 
  facet_grid(~class, scales = "free_x", space = "free_x", 
    labeller = labeller(class = c(DNA = "DNA Transposon", LINE = "LINE", LTR = "LTR", SINE = "SINE"))) + 
  theme(panel.spacing = unit(0,'lines'), strip.background = element_blank(), strip.text = element_text(size = 16))
repeat_sperm_recovery_barplot
ggsave(plot = repeat_sperm_recovery_barplot, "repeat_sperm_recovery_barplot.pdf", width = 16, height = 9)

# Load median Kimura score for each repeat sublclass
repeat_subclass_median_kimura_score = readRDS("~/genomes/nematostella/repeats/repeat_subclass_median_kimura_score.rds")

# Add median Kimura score as a column to repeat_sperm_recovery_df
repeat_sperm_recovery_df$median_kimura_score = repeat_subclass_median_kimura_score[repeat_sperm_recovery_df$region]

# Create scatter plots of Kimura score with mean methylation in DMSO
dmso_vs_kimura_scatter = plot_features(repeat_sperm_recovery_df$median_kimura_score, repeat_sperm_recovery_df$dmso_mean_methylation, 
  title = "Repeat Subclass Kimura Distance Vs Methylation\nin DMSO Samples", method = "s", position = "topright",  
  xlab = "Median Kimura Distance", ylab = "Mean Methylation", 
  label_size = 8, correlation_label = "Spearman Correlation = ")
ggsave(plot = dmso_vs_kimura_scatter, "dmso_vs_kimura_scatter.pdf", width = 9, height = 9)

# Create scatter plots of Kimura score with methylation recovery in sperm
sperm_vs_kimura_scatter = plot_features(repeat_sperm_recovery_df$median_kimura_score, repeat_sperm_recovery_df$recovery, 
  title = "Repeat Subclass Kimura Distance Vs\nMethylation Recovery in Sperm Samples", method = "s", position = "topright",  
  xlab = "Median Kimura Distance", ylab = "Proportion of Methylation Recovered", 
  label_size = 8, correlation_label = "Spearman Correlation = ")
ggsave(plot = sperm_vs_kimura_scatter, "sperm_vs_kimura_scatter.pdf", width = 9, height = 9)
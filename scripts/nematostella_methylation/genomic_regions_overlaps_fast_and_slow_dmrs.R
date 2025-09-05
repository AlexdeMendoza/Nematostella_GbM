### Compare genomic overlaps for fast and slow recovery DMRs

# Load required packages
library(GenomicRanges)
library(dplyr)
library(plotR)

# Split dmso_gr into slow and fast recovery DMRs
dmso_dmrs = readRDS("dmso_dmrs.rds")
slow_dmrs = dmso_dmrs[dmso_dmrs$dmr_status == "Slow"]
fast_dmrs = dmso_dmrs[dmso_dmrs$dmr_status == "Fast"]

### Find which superclasses of regions DMSO DMRs overlap

# Load Granges with genomic annotation of Nematostella and split into a GRangesList based on superclass of region
complete_annotation_gr = readRDS("~/nematostella_project/nematostella_genome/repeats/nvec_complete_annotation_gr.rds")

# Split complete_annotation_gr into a GRangesList based on superclass of region
complete_annotation_grl_superclass = split(complete_annotation_gr, complete_annotation_gr$superclass)

# Calculate the proportions of fast and slow DNRs overlapping each class of region and test the differences using Fisher's exact test
fast_vs_slow_dmr_proportions_superclass = data.frame(
  regions = names(complete_annotation_grl_superclass),
  slow_dmr = sapply(complete_annotation_grl_superclass, function(x) length(subsetByOverlaps(slow_dmrs, x)))/length(slow_dmrs),
  fast_dmr = sapply(complete_annotation_grl_superclass, function(x) length(subsetByOverlaps(fast_dmrs, x)))/length(fast_dmrs),
  p_value = sapply(complete_annotation_grl_superclass, function(x) fisher.test(rbind(
    table(factor(slow_dmrs %over% x, levels = c(FALSE, TRUE))),
    table(factor(fast_dmrs %over% x, levels = c(FALSE, TRUE)))
  ))$p.value),
  row.names = NULL
)

# Correct p-values and add significance symbol
fast_vs_slow_dmr_proportions_superclass$q_value = p.adjust(fast_vs_slow_dmr_proportions_superclass$p_value, method = "BH")
fast_vs_slow_dmr_proportions_superclass$significance = plotR::sig_sym(fast_vs_slow_dmr_proportions_superclass$q_value)
fast_vs_slow_dmr_proportions_superclass$significance_position = pmax(fast_vs_slow_dmr_proportions_superclass$slow_dmr, fast_vs_slow_dmr_proportions_superclass$fast_dmr) + 0.02

# Sort regions by those with the most significant difference
fast_vs_slow_dmr_proportions_superclass = filter(arrange(fast_vs_slow_dmr_proportions_superclass, q_value), q_value < 0.05)
fast_vs_slow_dmr_proportions_superclass$regions = factor(fast_vs_slow_dmr_proportions_superclass$regions, fast_vs_slow_dmr_proportions_superclass$regions)

# Convert fast_vs_slow_dmr_proportions_superclass to long format
fast_vs_slow_dmr_proportions_superclass_long = tidyr::pivot_longer(select(fast_vs_slow_dmr_proportions_superclass, regions, slow_dmr, fast_dmr), -regions)

# Remove tRNAs, snRNAs and rRNAs
fast_vs_slow_dmr_proportions_superclass_long = filter(fast_vs_slow_dmr_proportions_superclass_long, !regions %in% (c("tRNA", "snRNA", "rRNA")))

# Create a barplot comparing the enrichment of different classes of genomic regions among slow and fast recovery DMRs
region_superclass_enrichment_barplot = ggplot() +
  geom_col(data = fast_vs_slow_dmr_proportions_superclass_long, mapping = aes(y = value, x = regions, fill = name), color = "black", position = "dodge") +
  # geom_linerange(data = fast_vs_slow_dmr_proportions_superclass_long, mapping = aes(xmax = value, y = regions, group = name), xmin = 0, 
  #   linewidth = 1.5, position = position_dodge(width = 0.5)) +
  # geom_point(data = fast_vs_slow_dmr_proportions_superclass_long, mapping = aes(x = value, y = regions, fill = name), 
  #   size = 8, shape = 21, color = "black", position = position_dodge(width = 0.5)) +
  geom_text(data = fast_vs_slow_dmr_proportions_superclass, mapping = aes(label = significance, y = significance_position, x = regions), size = 10)
region_superclass_enrichment_barplot = customize_ggplot_theme(region_superclass_enrichment_barplot, 
  title = "Proportion of DMRs Overlapping Genomic Regions", xlab = "Genomic Regions", ylab = "Proportion", 
  #x_labels = c("PC Gene Body", "Repetitive DNA", "PC Promoter", "lncRNA", "CpG Islands", "Pseudogene"),
  fill_title = "DMR Group", fill_labels = c("Fast Recovery", "Slow Recovery")) # + scale_x_continuous(expand = expansion(mult = c(0, 0.2), add = 0))
region_superclass_enrichment_barplot
ggsave(plot = region_superclass_enrichment_barplot, "plots/dmr_region_superclass_enrichment_barplot.pdf", width = 16, height = 9)

### Find which classes of regions intragenic DMSO DMRs overlap

# Subset slow and fast DMRs for regions within the gene body
gene_body_gr = complete_annotation_gr[complete_annotation_gr$superclass == "Gene Body"]
slow_dmrs_gene_body = genomicTools::bedtools_intersect(slow_dmrs, gene_body_gr)
fast_dmrs_gene_body = genomicTools::bedtools_intersect(fast_dmrs, gene_body_gr)

# Split complete_annotation_gr into a GRangesList based on class of region
complete_annotation_grl_class = split(complete_annotation_gr, complete_annotation_gr$class)

# Calculate the proportions of fast and slow DNRs overlapping each class of region and test the differences using Fisher's exact test
fast_vs_slow_dmr_proportions_class = data.frame(
  regions = names(complete_annotation_grl_class),
  slow_dmr = sapply(complete_annotation_grl_class, function(x) length(subsetByOverlaps(slow_dmrs_gene_body, x))/length(slow_dmrs_gene_body)),
  fast_dmr = sapply(complete_annotation_grl_class, function(x) length(subsetByOverlaps(fast_dmrs_gene_body, x))/length(fast_dmrs_gene_body)),
  p_value = sapply(complete_annotation_grl_class, function(x) fisher.test(rbind(
    table(factor(slow_dmrs_gene_body %over% x, levels = c(FALSE, TRUE))),
    table(factor(fast_dmrs_gene_body %over% x, levels = c(FALSE, TRUE)))
  ))$p.value),
  row.names = NULL
)

# Correct p-values and add significance symbol
fast_vs_slow_dmr_proportions_class$q_value = p.adjust(fast_vs_slow_dmr_proportions_class$p_value, method = "BH")
fast_vs_slow_dmr_proportions_class$significance = plotR::sig_sym(fast_vs_slow_dmr_proportions_class$q_value)
fast_vs_slow_dmr_proportions_class$significance_position = pmax(fast_vs_slow_dmr_proportions_class$slow_dmr, fast_vs_slow_dmr_proportions_class$fast_dmr) + 0.02

# Filter for significant regions and sort regions by those with the most significant difference
fast_vs_slow_dmr_proportions_class = filter(arrange(fast_vs_slow_dmr_proportions_class, q_value), q_value < 0.05)
fast_vs_slow_dmr_proportions_class$regions = factor(fast_vs_slow_dmr_proportions_class$regions, fast_vs_slow_dmr_proportions_class$regions)

# Convert fast_vs_slow_dmr_proportions_class to long format
fast_vs_slow_dmr_proportions_class_long = tidyr::pivot_longer(select(fast_vs_slow_dmr_proportions_class, regions, slow_dmr, fast_dmr), -regions)

# Create a barplot comparing the enrichment of different classes of genomic regions among slow and fast recovery DMRs
region_class_enrichment_barplot = ggplot() +
  geom_col(data = fast_vs_slow_dmr_proportions_class_long, mapping = aes(y = value, x = regions, fill = name), color = "black", position = "dodge") +
  # geom_linerange(data = fast_vs_slow_dmr_proportions_class_long, mapping = aes(xmax = value, y = regions, group = name), xmin = 0, 
  #   linewidth = 1.5, position = position_dodge(width = 0.5)) +
  # geom_point(data = fast_vs_slow_dmr_proportions_class_long, mapping = aes(x = value, y = regions, fill = name), 
  #   size = 8, shape = 21, color = "black", position = position_dodge(width = 0.5)) +
  geom_text(data = fast_vs_slow_dmr_proportions_class, mapping = aes(label = significance, y = significance_position, x = regions), size = 10)
region_class_enrichment_barplot = customize_ggplot_theme(region_class_enrichment_barplot, 
  title = "Proportion of Intragenic DMRs Overlapping Genomic Regions", xlab = "Genomic Regions", ylab = "Proportion", 
  x_labels_angle = 45,
  fill_title = "DMR Group", fill_labels = c("Fast Recovery", "Slow Recovery")) # + scale_x_continuous(expand = expansion(mult = c(0, 0.2), add = 0))
region_class_enrichment_barplot
ggsave(plot = region_class_enrichment_barplot, "plots/dmr_region_class_enrichment_barplot.pdf", width = 16, height = 9)

### Find which subclasses of repeats intragenic DMSO DMRs overlap

# Create a GRanges just with repeats and create a vector matching repeat subfamilies to their family
repeat_annotation_gr = complete_annotation_gr[complete_annotation_gr$superclass == "Repetitive DNA"]

# Update the subclass of elements with LTR class and unknown subclass to LTR
repeat_annotation_gr[repeat_annotation_gr$subclass == "Unknown" & repeat_annotation_gr$class == "LTR"]$subclass = "LTR"

# Update the class of elements with tRNA subclass and SINE class to tRNA 
repeat_annotation_gr[repeat_annotation_gr$subclass == "tRNA" & repeat_annotation_gr$class == "SINE"]$class = "tRNA"

# Remove any regions with class as ARTEFACT or Unknown
repeat_annotation_gr = repeat_annotation_gr[!repeat_annotation_gr$class %in% c("ARTEFACT", "Unknown")]

# Create a vector matching repeats subclasses to their class
repeat_subclass_to_class = sapply(split(repeat_annotation_gr$class, repeat_annotation_gr$subclass), unique)

# Split complete_annotation_gr into a GRangesList based on subclass of region
repeat_annotation_grl_subclass = split(repeat_annotation_gr , repeat_annotation_gr$subclass)

# Calculate the proportions of fast and slow DNRs overlapping each subclass of region and test the differences using Fisher's exact test
fast_vs_slow_dmr_proportions_subclass = data.frame(
  regions = names(repeat_annotation_grl_subclass),
  slow_dmr = sapply(repeat_annotation_grl_subclass, function(x) length(subsetByOverlaps(slow_dmrs, x)))/length(slow_dmrs),
  fast_dmr = sapply(repeat_annotation_grl_subclass, function(x) length(subsetByOverlaps(fast_dmrs, x)))/length(fast_dmrs),
  p_value = sapply(repeat_annotation_grl_subclass, function(x) fisher.test(rbind(
    table(factor(slow_dmrs %over% x, levels = c(FALSE, TRUE))),
    table(factor(fast_dmrs %over% x, levels = c(FALSE, TRUE)))
  ))$p.value),
  row.names = NULL
)

# Correct p-values and add significance symbol
fast_vs_slow_dmr_proportions_subclass$q_value = p.adjust(fast_vs_slow_dmr_proportions_subclass$p_value, method = "BH")
fast_vs_slow_dmr_proportions_subclass$significance = plotR::sig_sym(fast_vs_slow_dmr_proportions_subclass$q_value)
fast_vs_slow_dmr_proportions_subclass$significance_position = pmax(fast_vs_slow_dmr_proportions_subclass$slow_dmr, fast_vs_slow_dmr_proportions_subclass$fast_dmr) + 0.005

# Sort regions by those with the biggest overlap with fast DMRs
# fast_vs_slow_dmr_proportions_subclass = arrange(fast_vs_slow_dmr_proportions_subclass, desc((slow_dmr + fast_dmr)/2))
fast_vs_slow_dmr_proportions_subclass = filter(dplyr::arrange(fast_vs_slow_dmr_proportions_subclass, q_value), q_value < 0.05)
fast_vs_slow_dmr_proportions_subclass$regions = factor(fast_vs_slow_dmr_proportions_subclass$regions, fast_vs_slow_dmr_proportions_subclass$regions)

# Convert fast_vs_slow_dmr_proportions_subclass to long format
fast_vs_slow_dmr_proportions_subclass_long = tidyr::pivot_longer(select(fast_vs_slow_dmr_proportions_subclass, regions, slow_dmr, fast_dmr), -regions)

# Add repeat family
fast_vs_slow_dmr_proportions_subclass_long$repeat_family = repeat_subclass_to_class[as.character(fast_vs_slow_dmr_proportions_subclass_long$regions)]
fast_vs_slow_dmr_proportions_subclass$repeat_family = repeat_subclass_to_class[as.character(fast_vs_slow_dmr_proportions_subclass$regions)]

# Filter for top 25 most enriched repeat subfamilies
fast_vs_slow_dmr_proportions_subclass_long_top25 = head(fast_vs_slow_dmr_proportions_subclass_long, 50)
fast_vs_slow_dmr_proportions_subclass_top25 = head(fast_vs_slow_dmr_proportions_subclass, 25)
fast_vs_slow_dmr_proportions_subclass_top25$regions = factor(fast_vs_slow_dmr_proportions_subclass_top25$regions, 
  levels = levels(fast_vs_slow_dmr_proportions_subclass_long_top25$regions))

# Create a barplot comparing the enrichment of different subclasses of genomic regions among slow and fast recovery DMRs
region_subclass_enrichment_barplot = ggplot() +
  geom_col(data = fast_vs_slow_dmr_proportions_subclass_long_top25, mapping = aes(y = value, x = regions, fill = name), color = "black", position = "dodge") +
  # geom_linerange(data = fast_vs_slow_dmr_proportions_subclass_long, mapping = aes(xmax = value, y = regions, group = name), xmin = 0, 
  #   linewidth = 1.5, position = position_dodge(width = 0.5)) +
  # geom_point(data = fast_vs_slow_dmr_proportions_subclass_long, mapping = aes(x = value, y = regions, fill = name), 
  #   size = 8, shape = 21, color = "black", position = position_dodge(width = 0.5)) +
  geom_text(data = fast_vs_slow_dmr_proportions_subclass_top25, mapping = aes(label = significance, y = significance_position, x = regions), size = 10)
region_subclass_enrichment_barplot = customize_ggplot_theme(region_subclass_enrichment_barplot, 
  title = "Proportion of DMRs Overlapping Intragenic Repeat Subfamilies", xlab = "Repeat Subfamily", ylab = "Proportion", x_labels_angle = 45,
  fill_title = "DMR Group", fill_labels = c("Fast Recovery", "Slow Recovery")) +
  facet_grid(~repeat_family, scales = "free_x", space = "free_x") + 
  theme(panel.spacing = unit(0,'lines'), strip.background = element_blank(),
    strip.text = element_text(size = 16))
region_subclass_enrichment_barplot
ggsave(plot = region_subclass_enrichment_barplot, "plots/dmr_region_subclass_enrichment_barplot.pdf", width = 16, height = 9)
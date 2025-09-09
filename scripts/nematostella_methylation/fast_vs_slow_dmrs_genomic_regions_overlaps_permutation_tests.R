### Compare genomic overlaps for fast and slow recovery DMRs

# Load required packages
library(GenomicRanges)
library(dplyr)
source("../auxillary_scripts/plotting_functions.R")
library(doParallel)

# Load DMSO DMRs
dmso_dmrs = readRDS("dmso_dmrs.rds")

cl = makeCluster(10)
registerDoParallel(cl)

# Define a function for testing the significance of the difference in overlap between slow and fast DMRs with another group of regions using a permutation test
slow_vs_fast_dmr_premutation_test = function(dmrs, query_gr, repetitions = 10000){
  
  # Split DMRs into slow and fast recovery DMRs
  slow_dmrs = dmrs[dmrs$dmr_status == "Slow"]
  fast_dmrs = dmrs[dmrs$dmr_status == "Fast"]
  
  # Calculate proportions of slow_dmrs and fast_dmrs overlapping query_gr
  slow_dmrs_overlap = genomicTools::calculate_regions_intersections(slow_dmrs, query_gr, overlap_measure = "proportion")
  fast_dmrs_overlap = genomicTools::calculate_regions_intersections(fast_dmrs, query_gr, overlap_measure = "proportion")
  fast_slow_difference = fast_dmrs_overlap - slow_dmrs_overlap
  
  # Calculate shuffled differences
  system.time({shuffled_differences = foreach(x = seq.int(repetitions), .combine = c, .export = c("dmso_dmrs", "slow_dmrs", "fast_dmrs")) %dopar% {
  
    # Shuffle labels of DMRs
    shuffled_dmrs = dmso_dmrs
    shuffled_dmrs$dmr_status = sample(shuffled_dmrs$dmr_status)
    
    # Separate shuffled DMRs into slow and fast DMRs
    shuffled_slow_dmrs = shuffled_dmrs[shuffled_dmrs$dmr_status == "Slow"]
    shuffled_fast_dmrs = shuffled_dmrs[shuffled_dmrs$dmr_status == "Fast"]
    
    # Calculate proportios of shuffled_slow_dmrs and shuffled_fast_dmrs overlapping query_gr
    shuffled_slow_dmrs_overlap = genomicTools::calculate_regions_intersections(shuffled_slow_dmrs, query_gr, overlap_measure = "proportion")
    shuffled_fast_dmrs_overlap = genomicTools::calculate_regions_intersections(shuffled_fast_dmrs, query_gr, overlap_measure = "proportion")
    
    # Return difference between fastt and slow DMRs
    shuffled_fast_dmrs_overlap - shuffled_slow_dmrs_overlap
    
  }})
  
  # Calculate p-value
  p_value = (sum(abs(shuffled_differences) >= abs(fast_slow_difference))+1)/(repetitions + 1)
  
  # Return slow_dmrs_overlap, fast_dmrs_overlap and p-value
  return(data.frame(slow_dmrs = slow_dmrs_overlap, fast_dmrs = fast_dmrs_overlap, p_value = p_value))
  
}

### Find which superclasses of regions DMSO DMRs overlap

# Load Granges with genomic annotation of Nematostella and split into a GRangesList based on superclass of region
complete_annotation_gr = readRDS("~/nematostella_project/nematostella_genome/repeats/nvec_complete_annotation_gr.rds")

# Split complete_annotation_gr into a GRangesList based on superclass of region
complete_annotation_grl_superclass = split(complete_annotation_gr, complete_annotation_gr$superclass)

# Test the difference in the overlap of slow and fast DMRs with region superclasses. Took 1 hour
set.seed(123)
system.time({fast_vs_slow_dmr_proportions_superclass = foreach(region_name = names(complete_annotation_grl_superclass), .combine = rbind) %do% {
  
  print(region_name)
  
  # Subset complete_annotation_grl_superclass for selected regions
  region_gr = complete_annotation_grl_superclass[[region_name]]
  
  # Test difference in overlap of slow and fast DMRs with regions
  result = slow_vs_fast_dmr_premutation_test(dmrs = dmso_dmrs, query_gr = region_gr, repetitions = 10000) 
  row.names(result) = region_name
  result = tibble::rownames_to_column(result, "regions")
  result
  
}})

# Save results
write.table(fast_vs_slow_dmr_proportions_superclass, "fast_vs_slow_dmr_proportions_superclass.tsv", sep = "\t", quote = F, row.names = F)
fast_vs_slow_dmr_proportions_superclass = data.table::fread("fast_vs_slow_dmr_proportions_superclass.tsv")

# Correct p-values and add significance symbol
fast_vs_slow_dmr_proportions_superclass$q_value = p.adjust(fast_vs_slow_dmr_proportions_superclass$p_value, method = "BH")
fast_vs_slow_dmr_proportions_superclass$significance = sig_sym(fast_vs_slow_dmr_proportions_superclass$q_value)
fast_vs_slow_dmr_proportions_superclass$significance_position = pmax(fast_vs_slow_dmr_proportions_superclass$slow_dmr, fast_vs_slow_dmr_proportions_superclass$fast_dmr) + 0.02

# Sort regions by those with the most significant difference
fast_vs_slow_dmr_proportions_superclass = filter(arrange(fast_vs_slow_dmr_proportions_superclass, q_value), q_value < 0.05)
fast_vs_slow_dmr_proportions_superclass$regions = factor(fast_vs_slow_dmr_proportions_superclass$regions, fast_vs_slow_dmr_proportions_superclass$regions)

# Convert fast_vs_slow_dmr_proportions_superclass to long format
fast_vs_slow_dmr_proportions_superclass_long = tidyr::pivot_longer(select(fast_vs_slow_dmr_proportions_superclass, regions, slow_dmrs, fast_dmrs), -regions)

# Create a barplot comparing the enrichment of different classes of genomic regions among slow and fast recovery DMRs
region_superclass_enrichment_barplot = ggplot() +
  geom_col(data = fast_vs_slow_dmr_proportions_superclass_long, mapping = aes(y = value, x = regions, fill = name), color = "black", position = "dodge") +
  geom_text(data = fast_vs_slow_dmr_proportions_superclass, mapping = aes(label = significance, y = significance_position, x = regions), size = 10)
region_superclass_enrichment_barplot = customize_ggplot_theme(region_superclass_enrichment_barplot, 
  title = "Proportion of DMRs Overlapping Genomic Regions", xlab = "Genomic Regions", ylab = "Proportion of Base Pairs", 
  fill_title = "DMR Group", fill_labels = c("Fast Recovery", "Slow Recovery")) 
region_superclass_enrichment_barplot
#ggsave(plot = region_superclass_enrichment_barplot, "plots/dmr_region_superclass_enrichment_barplot_permutation_test.pdf", width = 16, height = 9)

### Find which classes of regions intragenic DMSO DMRs overlap. Took 10 minutes

# Subset slow and fast DMRs for regions within the gene body
gene_body_gr = complete_annotation_gr[complete_annotation_gr$superclass == "Gene Body"]
gene_body_dmrs = genomicTools::bedtools_intersect(dmso_dmrs, gene_body_gr)

# Split complete_annotation_gr into a GRangesList based on class of region
complete_annotation_grl_class = split(complete_annotation_gr, complete_annotation_gr$class)

# Test the difference in the overlap of slow and fast DMRs with region classes
set.seed(123)
system.time({fast_vs_slow_dmr_proportions_class = foreach(region_name = names(complete_annotation_grl_class), .combine = rbind) %do% {
  
  print(region_name)
  
  # Subset complete_annotation_grl_class for selected regions
  region_gr = complete_annotation_grl_class[[region_name]]
  
  # Test difference in overlap of slow and fast DMRs with regions
  result = slow_vs_fast_dmr_premutation_test(dmrs = gene_body_dmrs, query_gr = region_gr, repetitions = 10000) 
  row.names(result) = region_name
  result = tibble::rownames_to_column(result, "regions")
  result
  
}})

# Save results
write.table(fast_vs_slow_dmr_proportions_class, "fast_vs_slow_dmr_proportions_class.tsv", sep = "\t", quote = F, row.names = F)
fast_vs_slow_dmr_proportions_class = data.table::fread("fast_vs_slow_dmr_proportions_class.tsv")

# Correct p-values and add significance symbol
fast_vs_slow_dmr_proportions_class$q_value = p.adjust(fast_vs_slow_dmr_proportions_class$p_value, method = "BH")
fast_vs_slow_dmr_proportions_class$significance = sig_sym(fast_vs_slow_dmr_proportions_class$q_value)
fast_vs_slow_dmr_proportions_class$significance_position = pmax(fast_vs_slow_dmr_proportions_class$slow_dmr, fast_vs_slow_dmr_proportions_class$fast_dmr) + 0.02

# Sort regions by those with the most significant difference
fast_vs_slow_dmr_proportions_class = filter(arrange(fast_vs_slow_dmr_proportions_class, q_value), q_value < 0.05)
fast_vs_slow_dmr_proportions_class$regions = factor(fast_vs_slow_dmr_proportions_class$regions, fast_vs_slow_dmr_proportions_class$regions)

# Convert fast_vs_slow_dmr_proportions_class to long format
fast_vs_slow_dmr_proportions_class_long = tidyr::pivot_longer(select(fast_vs_slow_dmr_proportions_class, regions, slow_dmrs, fast_dmrs), -regions)

# Remove tRNAs, snRNAs and rRNAs
fast_vs_slow_dmr_proportions_class_long = filter(fast_vs_slow_dmr_proportions_class_long, !regions %in% (c("tRNA", "snRNA", "rRNA")))

# Create a barplot comparing the enrichment of different classes of genomic regions among slow and fast recovery DMRs
region_class_enrichment_barplot = ggplot() +
  geom_col(data = fast_vs_slow_dmr_proportions_class_long, mapping = aes(y = value, x = regions, fill = name), color = "black", position = "dodge") +
  geom_text(data = fast_vs_slow_dmr_proportions_class, mapping = aes(label = significance, y = significance_position, x = regions), size = 10)
region_class_enrichment_barplot = customize_ggplot_theme(region_class_enrichment_barplot, 
  title = "Proportion of Intragenic DMRs Overlapping Genomic Regions", xlab = "Genomic Regions", ylab = "Proportion of Base Pairs", 
  x_labels_angle = 45, fill_title = "DMR Group", fill_labels = c("Fast Recovery", "Slow Recovery")) 
region_class_enrichment_barplot
#ggsave(plot = region_class_enrichment_barplot, "plots/dmr_region_class_enrichment_barplot_permutation_test.pdf", width = 16, height = 9)

### Find which subfamilies of repeats intragenic DMSO DMRs overlap

# Create a GRanges just with repeats and create a vector matching repeat subfamilies to their family
repeat_annotation_gr = complete_annotation_gr[complete_annotation_gr$superclass== "Repetitive DNA"]

# Update the subfamily of elements with LTR family and unknown subfamily to LTR
repeat_annotation_gr[repeat_annotation_gr$subclass == "Unknown" & repeat_annotation_gr$class == "LTR"]$subclass = "LTR"

# Update the family of elements with tRNA subfamily and SINE family to tRNA 
repeat_annotation_gr[repeat_annotation_gr$subclass == "tRNA" & repeat_annotation_gr$class == "SINE"]$class = "tRNA"

# Remove any repeats with family as ARTEFACT, Unknown or Simple Repeat
repeat_annotation_gr = repeat_annotation_gr[!repeat_annotation_gr$class %in% c("ARTEFACT", "Unknown", "Simple Repeat")]

# Create a vector matching repeats subfamilies to their family
repeat_subfamily_to_family = sapply(split(repeat_annotation_gr$class, repeat_annotation_gr$subclass), unique)

# Split complete_annotation_gr into a GRangesList based on subfamily of repeat
repeat_annotation_grl_subfamily = split(repeat_annotation_gr , repeat_annotation_gr$subclass)

# Test the difference in the overlap of slow and fast DMRs with repeat subfamilies. Took 7 hours with 10000 permutations and 10 cores
set.seed(123)
system.time({fast_vs_slow_dmr_repeat_subfamily_proportions = foreach(repeat_name = names(repeat_annotation_grl_subfamily), .combine = rbind) %do% {
  
  print(repeat_name)
  
  # Subset complete_annotation_grl_subfamily for selected repeats
  repeat_gr = complete_annotation_grl_subfamily[[repeat_name]]
  
  # Test difference in overlap of slow and fast DMRs with repeats
  result = slow_vs_fast_dmr_premutation_test(dmrs = gene_body_dmrs, query_gr = repeat_gr, repetitions = 10000) 
  row.names(result) = repeat_name
  result = tibble::rownames_to_column(result, "repeats")
  result
  
}})

# Save results
write.table(fast_vs_slow_dmr_repeat_subfamily_proportions, "fast_vs_slow_dmr_repeat_subfamily_proportions.tsv", sep = "\t", quote = F, row.names = F)
fast_vs_slow_dmr_repeat_subfamily_proportions = data.table::fread("fast_vs_slow_dmr_repeat_subfamily_proportions.tsv")

# Correct p-values and add significance symbol
fast_vs_slow_dmr_repeat_subfamily_proportions$q_value = p.adjust(fast_vs_slow_dmr_repeat_subfamily_proportions$p_value, method = "BH")
fast_vs_slow_dmr_repeat_subfamily_proportions$significance = sig_sym(fast_vs_slow_dmr_repeat_subfamily_proportions$q_value)
fast_vs_slow_dmr_repeat_subfamily_proportions$significance_position = pmax(fast_vs_slow_dmr_repeat_subfamily_proportions$slow_dmr, fast_vs_slow_dmr_repeat_subfamily_proportions$fast_dmr) + 0.0005

# Sort repeats by those with the most significant difference
fast_vs_slow_dmr_repeat_subfamily_proportions = filter(arrange(fast_vs_slow_dmr_repeat_subfamily_proportions, q_value), q_value < 0.05)
fast_vs_slow_dmr_repeat_subfamily_proportions$repeats = factor(fast_vs_slow_dmr_repeat_subfamily_proportions$repeats, fast_vs_slow_dmr_repeat_subfamily_proportions$repeats)

# Convert fast_vs_slow_dmr_repeat_subfamily_proportions to long format
fast_vs_slow_dmr_repeat_subfamily_proportions_long = tidyr::pivot_longer(select(fast_vs_slow_dmr_repeat_subfamily_proportions, repeats, slow_dmrs, fast_dmrs), -repeats)

# Remove tRNAs, snRNAs and rRNAs
fast_vs_slow_dmr_repeat_subfamily_proportions_long = filter(fast_vs_slow_dmr_repeat_subfamily_proportions_long, !repeats %in% (c("tRNA", "snRNA", "rRNA")))

# Add repeat family
fast_vs_slow_dmr_repeat_subfamily_proportions_long$repeat_family = repeat_subfamily_to_family[as.character(fast_vs_slow_dmr_repeat_subfamily_proportions_long$repeats)]
fast_vs_slow_dmr_repeat_subfamily_proportions$repeat_family = repeat_subfamily_to_family[as.character(fast_vs_slow_dmr_repeat_subfamily_proportions$repeats)]

# Create a barplot comparing the enrichment of different subfamilies of repeats among slow and fast recovery DMRs
repeat_subfamily_enrichment_barplot = ggplot() +
  geom_col(data = fast_vs_slow_dmr_repeat_subfamily_proportions_long, mapping = aes(y = value, x = repeats, fill = name), color = "black", position = "dodge") +
  geom_text(data = fast_vs_slow_dmr_repeat_subfamily_proportions, mapping = aes(label = significance, y = significance_position, x = repeats), size = 10)
repeat_subfamily_enrichment_barplot = customize_ggplot_theme(repeat_subfamily_enrichment_barplot, 
  title = "Proportion of Intragenic DMRs Overlapping Genomic Repeat Subfamilies", xlab = "Repeat Subfamily", ylab = "Proportion of Base Pairs", 
  x_labels_angle = 45, fill_title = "DMR Group", fill_labels = c("Fast Recovery", "Slow Recovery")) +
  facet_grid(~repeat_family, scales = "free_x", space = "free_x") + 
  theme(panel.spacing = unit(0,'lines'), strip.background = element_blank(),
    strip.text = element_text(size = 16))
repeat_subfamily_enrichment_barplot
#ggsave(plot = repeat_subfamily_enrichment_barplot, "plots/dmr_repeat_subfamily_enrichment_barplot_permutation_test.pdf", width = 16, height = 9)

### Test enrichment of population DMRs for different genomic regions

# Load required packages
library(dplyr)
library(BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19)
library(regioneR)
library(plotR)

# Load population DMRs
population_dmrs_reduced = readRDS("population_dmrs_reduced.rds")

# Load DMSO DMRs and split into a list
dmso_dmrs = readRDS("../nematostella_methylation/dmso_dmrs.rds")
dmso_dmrs = split(dmso_dmrs, dmso_dmrs$dmr_status)
names(dmso_dmrs) = c("Fast DMRs", "Slow DMRs")

# Create a GRanges from the Nematostella genome
genome = GRanges(seqinfo(BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19))

# Load GRanges with genomic annotation of Nematostella and split into GRangesLists based on superclass, class and subclass of regions
complete_annotation_gr = readRDS("~/nematostella_project/nematostella_genome/repeats/nvec_complete_annotation_gr.rds")
complete_annotation_grl_superclass = split(complete_annotation_gr, complete_annotation_gr$superclass)
complete_annotation_grl_class = split(complete_annotation_gr, complete_annotation_gr$class)
complete_annotation_grl_subclass = split(complete_annotation_gr, complete_annotation_gr$subclass)

# Add DMSO DMRs to complete_annotation_grl_superclass
complete_annotation_grl_superclass = c(complete_annotation_grl_superclass, dmso_dmrs)

# Create a function to test the overlap of DMRs with different classes of genomic regions
dmr_overlap_permutation_test = function(dmrs, regions_grl, n, nthreads = 5){
  
  # Make sure names of regions_grl are in order
  regions_grl = regions_grl[sort(names(regions_grl))]
  
  # Calculate size of the overlaps between dmrs and annotation regions
  dmr_overlaps = sapply(regions_grl, function(x) sum(width(GenomicRanges::intersect(x, dmrs, ignore.strand = T))))
  
  # Flatten regions_grl
  regions_gr = unlist(GRangesList(lapply(regions_grl, unname)))
  
  # Create a GRangeLsit with randomized dmrs
  system.time({random_dmrs = GRangesList(parallel::mclapply(seq.int(n), function(x)
    randomizeRegions(A = dmrs, genome = genome, allow.overlaps = F, per.chromosome = F), mc.cores = nthreads))})
  
  # Name GRanges with permuation number and convert to a flat GRanges
  names(random_dmrs) = paste0("permutation_", seq_along(random_dmrs))
  random_dmrs = unlist(random_dmrs)
  names(random_dmrs) = gsub("\\..*", "", names(random_dmrs))
  
  # Find overlaps between random_dmrs and regions_gr
  random_overlaps = data.frame(findOverlaps(random_dmrs, regions_gr, ignore.strand = T))
  
  # Add permutation number and region type
  random_overlaps$permutation = names(random_dmrs)[random_overlaps$queryHits]
  random_overlaps$region_type = names(regions_gr)[random_overlaps$subjectHits]
  
  # Find size of the intersections
  random_overlaps$intersection = width(pintersect(random_dmrs[random_overlaps$queryHits], regions_gr[random_overlaps$subjectHits]))
  
  # Set region_type as a factor
  random_overlaps$region_type = factor(random_overlaps$region_type, levels = names(regions_grl))
  
  # Summarize the intersections
  random_overlaps_summary = data.frame(dplyr::summarise(
    dplyr::group_by(random_overlaps, permutation, region_type, .drop = FALSE), 
      intersection = sum(intersection)))
  
  # Convert into a wide data.frame and ensure columns are in the same order as regions_grl
  random_overlaps_summary = tidyr::pivot_wider(random_overlaps_summary, names_from = "region_type", values_from = "intersection")
  random_overlaps_summary = data.frame(tibble::column_to_rownames(random_overlaps_summary, "permutation"))
  names(random_overlaps_summary) = names(regions_grl)
  random_overlaps_summary = random_overlaps_summary[, names(regions_grl)]
  
  # Calculate p-values
  p_value = (rowSums(apply(random_overlaps_summary, 1, function(x) x >= dmr_overlaps)) + 1)/(n+1)
  
  # Create a data.frame with the results and return
  results_df = data.frame(
    region_type = names(dmr_overlaps),
    dmr_overlaps = dmr_overlaps,
    random_overlaps = colMeans(random_overlaps_summary, na.rm =T),
    enrichment = dmr_overlaps/colMeans(random_overlaps_summary, na.rm = T),
    p_value = p_value, row.names = NULL
  )
  return(results_df)
  
}

# Test the enrichment of populatiion DMRs in region superclasses
set.seed(123)
system.time({population_dmr_region_superclass_enrichment = 
  dmr_overlap_permutation_test(dmrs = population_dmrs_reduced, regions_grl = complete_annotation_grl_superclass, n = 3000)})
population_dmr_region_superclass_enrichment$q_value = p.adjust(population_dmr_region_superclass_enrichment$p_value, method = "fdr")
population_dmr_region_superclass_enrichment$significance = plotR::sig_sym(population_dmr_region_superclass_enrichment$q_value)

# Sort regions by those with the most significant difference
population_dmr_region_superclass_enrichment = arrange(population_dmr_region_superclass_enrichment, q_value)
population_dmr_region_superclass_enrichment$region_type = factor(population_dmr_region_superclass_enrichment$region_type, population_dmr_region_superclass_enrichment$region_type)

# Convert population_dmr_region_superclass_enrichment to long format
population_dmr_region_superclass_enrichment_long = tidyr::pivot_longer(select(population_dmr_region_superclass_enrichment, region_type, enrichment), -region_type)

# Create a barplot comparing the enrichment of different superclasses of genomic regions among stable and variable recovery DMRs
region_superclass_enrichment_barplot = ggplot() +
  geom_col(data = population_dmr_region_superclass_enrichment_long, mapping = aes(y = value, x = region_type), color = "black", position = "dodge") +
  geom_text(data = population_dmr_region_superclass_enrichment, mapping = aes(label = significance, y = enrichment + 0.05, x = region_type), size = 10)
region_superclass_enrichment_barplot = customize_ggplot_theme(region_superclass_enrichment_barplot, 
  title = "Proportion of DMRs Overlapping Genomic Regions", xlab = "Genomic Region", ylab = "Enrichment", x_labels_angle = 45)
region_superclass_enrichment_barplot
ggsave(plot = region_superclass_enrichment_barplot, "region_superclass_enrichment_barplot.pdf", width = 16, height = 9)

# Test the enrichment of populatiion DMRs in region classes
set.seed(123)
system.time({population_dmr_region_class_enrichment = 
  dmr_overlap_permutation_test(dmrs = population_dmrs_reduced, regions_grl = complete_annotation_grl_class, n = 3000)})
population_dmr_region_class_enrichment$q_value = p.adjust(population_dmr_region_class_enrichment$p_value, method = "fdr")
population_dmr_region_class_enrichment$significance = plotR::sig_sym(population_dmr_region_class_enrichment$q_value)

# Sort regions by those with the most significant difference
population_dmr_region_class_enrichment = arrange(population_dmr_region_class_enrichment, q_value)
population_dmr_region_class_enrichment$region_type = factor(population_dmr_region_class_enrichment$region_type, population_dmr_region_class_enrichment$region_type)

# Convert population_dmr_region_class_enrichment to long format
population_dmr_region_class_enrichment_long = tidyr::pivot_longer(select(population_dmr_region_class_enrichment, region_type, enrichment), -region_type)

# Create a barplot comparing the enrichment of different classes of genomic regions among stable and variable recovery DMRs
region_class_enrichment_barplot = ggplot() +
  geom_col(data = population_dmr_region_class_enrichment_long, mapping = aes(y = value, x = region_type), color = "black", position = "dodge") +
  geom_text(data = population_dmr_region_class_enrichment, mapping = aes(label = significance, y = enrichment + 0.05, x = region_type), size = 10)
region_class_enrichment_barplot = customize_ggplot_theme(region_class_enrichment_barplot, 
  title = "Proportion of DMRs Overlapping Genomic Regions", xlab = "Genomic Region", ylab = "Enrichment", x_labels_angle = 45)
region_class_enrichment_barplot
ggsave(plot = region_class_enrichment_barplot, "region_class_enrichment_barplot.pdf", width = 16, height = 9)

# Test the enrichment of populatiion DMRs in region subclasses
set.seed(123)
system.time({population_dmr_region_subclass_enrichment = 
  dmr_overlap_permutation_test(dmrs = population_dmrs_reduced, regions_grl = complete_annotation_grl_subclass, n = 1000)})
population_dmr_region_subclass_enrichment$q_value = p.adjust(population_dmr_region_subclass_enrichment$p_value, method = "fdr")
population_dmr_region_subclass_enrichment$significance = plotR::sig_sym(population_dmr_region_subclass_enrichment$q_value)

# Sort regions by those with the most significant difference
population_dmr_region_subclass_enrichment_sig = filter(arrange(population_dmr_region_subclass_enrichment, q_value), q_value < 0.2)
population_dmr_region_subclass_enrichment_sig$region_type = factor(population_dmr_region_subclass_enrichment_sig$region_type, population_dmr_region_subclass_enrichment_sig$region_type)

# Convert population_dmr_region_subclass_enrichment_sig to long format
population_dmr_region_subclass_enrichment_sig_long = tidyr::pivot_longer(select(population_dmr_region_subclass_enrichment_sig, region_type, enrichment), -region_type)

# Create a barplot comparing the enrichment of different classes of genomic regions among stable and variable recovery DMRs
region_subclass_enrichment_barplot = ggplot() +
  geom_col(data = population_dmr_region_subclass_enrichment_sig_long, mapping = aes(y = value, x = region_type), color = "black", position = "dodge") +
  geom_text(data = population_dmr_region_subclass_enrichment_sig, mapping = aes(label = significance, y = enrichment + 0.05, x = region_type), size = 10)
region_subclass_enrichment_barplot = customize_ggplot_theme(region_subclass_enrichment_barplot, 
  title = "Proportion of DMRs Overlapping Genomic Regions", xlab = "Genomic Region", ylab = "Enrichment", x_labels_angle = 45)
region_subclass_enrichment_barplot
ggsave(plot = region_subclass_enrichment_barplot, "region_subclass_enrichment_barplot.pdf", width = 16, height = 9)

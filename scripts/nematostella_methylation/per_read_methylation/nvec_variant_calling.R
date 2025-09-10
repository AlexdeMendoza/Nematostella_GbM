# Create a pedigree phased VCF for Nv_Cr5_3_4mpf based on which variants are homozygous in NvIT16K4_sperm_Oct31 but heterozygous in Nv_Cr5_3_4mpf
# Nv_Cr5_3_4mpf is an individual offspring from the cross of a demethylated male and an untreated female; NvIT16K4_sperm_Oct31 is sperm from the demethylated male

# Load required packages
library(VariantAnnotation)
library(dplyr)

# Read in phased variants for Nv_Cr5_3_4mpf and NvIT16K4_sperm_Oct31
Nv_Cr5_3_4mpf_clair3_and_NvIT16K4_sperm_Oct31 = readVcf("Nv_Cr5_3_4mpf_clair3_and_NvIT16K4_sperm_Oct31_merged_phased_vcf.gz")

# Get genotypes for Nv_Cr5_3_4mpf_clair3_and_NvIT16K4_sperm_Oct31
genotypes = data.frame(variant_name = names(Nv_Cr5_3_4mpf_clair3_and_NvIT16K4_sperm_Oct31), 
  geno(Nv_Cr5_3_4mpf_clair3_and_NvIT16K4_sperm_Oct31)$GT, row.names = NULL)

# Find sites that are homozygous in the sperm (0/0 or 1/1), but heterozygous in the offspring
genotypes = filter(genotypes, 
  (NvIT16K4_sperm_Oct31 == "0/0" & Nv_Cr5_3 %in% c("0/1", "0|1", "1|0")) | (NvIT16K4_sperm_Oct31 == "1/1" & Nv_Cr5_3 %in% c("0/1", "0|1", "1|0")))

# Get names of the variants associated with these sites
selected_variants = genotypes$variant_name

# If sperm is homozygous reference, alternate allele must have come from the mother so set genotype as 0|1
# If the sperm is homozygous for the alternate allele, alternate allele must have come from the father so set genotype as 1|0
genotypes$Nv_Cr5_3[genotypes$NvIT16K4_sperm_Oct31 == "0/0"] = "0|1"
genotypes$Nv_Cr5_3[genotypes$NvIT16K4_sperm_Oct31 == "1/1"] = "1|0"

# Convert genotypes back to a matrix
genotypes_matrix = as.matrix(dplyr::select(genotypes, -variant_name))
row.names(genotypes_matrix) = selected_variants

# Filter Nv_Cr5_3_4mpf_clair3_and_NvIT16K4_sperm_Oct31 for the selected variants in Nv_Cr5_3
system.time({Nv_Cr5_3_4mpf_filtered = 
  Nv_Cr5_3_4mpf_clair3_and_NvIT16K4_sperm_Oct31[names(Nv_Cr5_3_4mpf_clair3_and_NvIT16K4_sperm_Oct31) %in% selected_variants, "Nv_Cr5_3"]})

# Update the genotypes of Nv_Cr5_3_4mpf_clair3_and_NvIT16K4_sperm_Oct31_filtered 
geno(Nv_Cr5_3_4mpf_filtered) = genotypes_matrix[, "Nv_Cr5_3", drop = F]

# Output VCF
system.time({writeVcf(Nv_Cr5_3_4mpf_filtered, filename = "Nv_Cr5_3_4mpf_clair3_output/pedigree_phased_merge_output.vcf.gz", index = T)})

# Extract GRanges from Nv_Cr5_3_4mpf_filtered and add genotype of Nv_Cr5_3
Nv_Cr5_3_4mpf_filtered_gr = rowRanges(Nv_Cr5_3_4mpf_filtered)
Nv_Cr5_3_4mpf_filtered_gr$Nv_Cr5_3_4mpf_genotype = genotypes$Nv_Cr5_3

# Convert REF and ALT to characters
Nv_Cr5_3_4mpf_filtered_gr$REF = as.character(Nv_Cr5_3_4mpf_filtered_gr$REF)
Nv_Cr5_3_4mpf_filtered_gr$ALT = as.character(unlist(Nv_Cr5_3_4mpf_filtered_gr$ALT))

# Assign REF and ALT to paternal or maternal alleles depending on the genotype
Nv_Cr5_3_4mpf_filtered_gr$paternal_allele = ifelse(Nv_Cr5_3_4mpf_filtered_gr$Nv_Cr5_3_4mpf_genotype == "0|1", Nv_Cr5_3_4mpf_filtered_gr$REF, Nv_Cr5_3_4mpf_filtered_gr$ALT)
Nv_Cr5_3_4mpf_filtered_gr$maternal_allele = ifelse(Nv_Cr5_3_4mpf_filtered_gr$Nv_Cr5_3_4mpf_genotype == "0|1", Nv_Cr5_3_4mpf_filtered_gr$ALT, Nv_Cr5_3_4mpf_filtered_gr$REF)

# Find overlaps between variants and exons
nvec_exons_gr = readRDS("../../nematostella_genome/nvec_pc_transcripts_promoters_exons_and_introns_gr.rds")
nvec_exons_gr = nvec_exons_gr[nvec_exons_gr$region ==  "exon"]
nvec_transcript_gr = readRDS("../../nematostella_genome/nvec_pc_transcripts_gr.rds")
variant_exon_overlaps = data.frame(findOverlaps(Nv_Cr5_3_4mpf_filtered_gr, nvec_exons_gr))

# Create a GRanges for variants overlapping exons
variant_exon_overlaps_gr = granges(Nv_Cr5_3_4mpf_filtered_gr)[variant_exon_overlaps$queryHits]
names(variant_exon_overlaps_gr) = NULL

# Add gene, transcript, Nv_Cr5_3 genotype, reference allele, paternal allele, maternal allele and call quality to variant_exon_overlaps_gr
variant_exon_overlaps_gr$gene = nvec_exons_gr$gene[variant_exon_overlaps$subjectHits]
variant_exon_overlaps_gr$transcript = nvec_exons_gr$transcript_id[variant_exon_overlaps$subjectHits]
variant_exon_overlaps_gr$Nv_Cr5_3_4mpf_genotype = Nv_Cr5_3_4mpf_filtered_gr$Nv_Cr5_3_4mpf_genotype[variant_exon_overlaps$queryHits]
variant_exon_overlaps_gr$reference_allele = Nv_Cr5_3_4mpf_filtered_gr$REF[variant_exon_overlaps$queryHits]
variant_exon_overlaps_gr$paternal_allele = Nv_Cr5_3_4mpf_filtered_gr$paternal_allele[variant_exon_overlaps$queryHits]
variant_exon_overlaps_gr$maternal_allele = Nv_Cr5_3_4mpf_filtered_gr$maternal_allele[variant_exon_overlaps$queryHits]
variant_exon_overlaps_gr$call_quality = Nv_Cr5_3_4mpf_filtered_gr$QUAL[variant_exon_overlaps$queryHits]

# Sort variant_exon_overlaps_gr and save
variant_exon_overlaps_gr = sort(variant_exon_overlaps_gr)
saveRDS(variant_exon_overlaps_gr, "variant_exon_overlaps_gr.rds")
variant_exon_overlaps_gr =  readRDS("variant_exon_overlaps_gr.rds")
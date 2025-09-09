# Create CGmap files files from modBAM files

# Load required packages
library(dplyr)
source("../auxillary_scripts/modkit_functions.R")

# Get paths to modBAM files
modbams = gsub(" .*", "", readLines("nanopore_data.txt"))

# Create names for CGmap files
cgmap_names = gsub("_sorted.bam", ".CGmap.gz", modbams)

# Create a GRanges with methylation pileup for for Nv_Cr5_3_4mpf
system.time({Nv_Cr5_3_4mpf_pileup_gr = pileup_modcalls(modBAM = "/data/SBCS-ademendoza/07-Lan/Nanopore/Nv_Cr5_3_4mpf/Nv_Cr5_3_4mpf_sorted.bam",
  reference_fasta = "~/genomes/nematostella/Nematostella_DToL_lambda_pUC_mitochondria_originalIDs.fasta", motif = "CG",
  ignore_mod = "h", combine_strands = F, nthreads = 50)})

# Convert Nv_Cr5_3_4mpf_pileup_gr into a data.frame and format it so that it can be written as a CGmap file
Nv_Cr5_3_4mpf_pileup_df = data.frame(Nv_Cr5_3_4mpf_pileup_gr)
Nv_Cr5_3_4mpf_pileup_df = transmute(Nv_Cr5_3_4mpf_pileup_df, seqnames = seqnames, nucleotide = "C", start = start,
  context = "CG", dinucleotide = "CG", methylated_proportion = percent_modified/100,
  methylated_counts = count_modified, total_counts = valid_coverage)
data.table::fwrite(Nv_Cr5_3_4mpf_pileup_df, "cgmap_files/cgmap_files_nanopore/Nv_Cr5_3_4mpf.CGmap.gz",
  sep = "\t", col.names = F, quote = F)

# Create a GRanges with methylation pileup for for NvIT16K4_sperm_Oct31
system.time({NvIT16K4_sperm_Oct31_pileup_gr = pileup_modcalls(modBAM = "/data/SBCS-ademendoza/07-Lan/Nanopore/NvIT16K4_sperm_Oct31/NvIT16K4_sperm_Oct31_sorted.bam",
  reference_fasta = "~/genomes/nematostella/Nematostella_DToL_lambda_pUC_mitochondria_originalIDs.fasta", motif = "CG",
  ignore_mod = "h", combine_strands = F, nthreads = 50)})

# Convert NvIT16K4_sperm_Oct31_gr into a data.frame and format it so that it can be written as a CGmap file
NvIT16K4_sperm_Oct31_pileup_df = data.frame(NvIT16K4_sperm_Oct31_pileup_gr)
NvIT16K4_sperm_Oct31_pileup_df = transmute(NvIT16K4_sperm_Oct31_pileup_df, seqnames = seqnames, nucleotide = "C", start = start,
  context = "CG", dinucleotide = "CG", methylated_proportion = percent_modified/100,
  methylated_counts = count_modified, total_counts = valid_coverage)
data.table::fwrite(NvIT16K4_sperm_Oct31_pileup_df, "cgmap_files/cgmap_files_nanopore/NvIT16K4_sperm_Oct31.CGmap.gz",
  sep = "\t", col.names = F, quote = F)

# Create a GRanges with methylation pileup for for Nv_Cr5
system.time({Nv_Cr5_pileup_gr = pileup_modcalls(modBAM = "/data/SBCS-ademendoza/07-Lan/Nanopore/Nv_Cr5/Nv_Cr5_sorted.bam",
  reference_fasta = "~/genomes/nematostella/Nematostella_DToL_lambda_pUC_mitochondria_originalIDs.fasta", motif = "CG",
  ignore_mod = "h", combine_strands = F, nthreads = 50)})

# Convert Nv_Cr5_pileup_gr into a data.frame and format it so that it can be written as a CGmap file
Nv_Cr5_pileup_df = data.frame(Nv_Cr5_pileup_gr)
Nv_Cr5_pileup_df = transmute(Nv_Cr5_pileup_df, seqnames = seqnames, nucleotide = "C", start = start,
  context = "CG", dinucleotide = "CG", methylated_proportion = percent_modified/100,
  methylated_counts = count_modified, total_counts = valid_coverage)
data.table::fwrite(Nv_Cr5_pileup_df, "cgmap_files/cgmap_files_nanopore//Nv_Cr5.CGmap.gz",
  sep = "\t", col.names = F, quote = F)

# Create a GRanges with methylation pileup for for Nv_Cr18
system.time({Nv_Cr18_pileup_gr = pileup_modcalls(modBAM = "/data/SBCS-ademendoza/07-Lan/Nanopore/Nv_Cr18/Nv_Cr18_sorted.bam",
  reference_fasta = "~/genomes/nematostella/Nematostella_DToL_lambda_pUC_mitochondria_originalIDs.fasta", motif = "CG",
  ignore_mod = "h", combine_strands = F, nthreads = 50)})

# Convert Nv_Cr18_pileup_gr into a data.frame and format it so that it can be written as a CGmap file
Nv_Cr18_pileup_df = data.frame(Nv_Cr18_pileup_gr)
Nv_Cr18_pileup_df = transmute(Nv_Cr18_pileup_df, seqnames = seqnames, nucleotide = "C", start = start,
  context = "CG", dinucleotide = "CG", methylated_proportion = percent_modified/100,
  methylated_counts = count_modified, total_counts = valid_coverage)
data.table::fwrite(Nv_Cr18_pileup_df, "cgmap_files/cgmap_files_nanopore/Nv_Cr18.CGmap.gz",
  sep = "\t", col.names = F, quote = F)

# Create a GRanges with methylation pileup for for 3_mpf_15X (combined skimming samples). Took 9 minutes with 1 core
system.time({three_mpf_15X_pileup_gr = pileup_modcalls(modBAM = "/data/SBCS-ademendoza/07-Lan/Nanopore/RBK_combined/GSK_F0_juvenile_combined/3mpf_15X_sorted.bam", 
  reference_fasta = "~/genomes/nematostella/Nematostella_DToL_lambda_pUC_mitochondria_originalIDs.fasta", motif = "CG", 
  ignore_mod = "h", combine_strands = F, nthreads = 1)})

# Convert three_mpf_15X_pileup_gr into a data.frame and format it so that it can be written as a CGmap file 
three_mpf_15X_pileup_df = data.frame(three_mpf_15X_pileup_gr)
three_mpf_15X_pileup_df = transmute(three_mpf_15X_pileup_df, seqnames = seqnames, nucleotide = "C", start = start,
  context = "CG", dinucleotide = "CG", methylated_proportion = percent_modified/100, 
  methylated_counts = count_modified, total_counts = valid_coverage)
data.table::fwrite(three_mpf_15X_pileup_df, "cgmap_files/cgmap_files_nanopore/three_mpf_15X.CGmap.gz", 
  sep = "\t", col.names = F, quote = F)

# Create a GRanges with methylation pileup for for 5_mpf_11X (combined skimming samples). Took 7 minutes
system.time({five_mpf_11X_pileup_gr = pileup_modcalls(modBAM = "/data/SBCS-ademendoza/07-Lan/Nanopore/RBK_combined/GSK_F0_juvenile_combined/5mpf_11X_sorted.bam", 
  reference_fasta = "~/genomes/nematostella/Nematostella_DToL_lambda_pUC_mitochondria_originalIDs.fasta", motif = "CG", 
  ignore_mod = "h", combine_strands = F, nthreads = 1)})

# Convert five_mpf_11X_pileup_gr into a data.frame and format it so that it can be written as a CGmap file 
five_mpf_11X_pileup_df = data.frame(five_mpf_11X_pileup_gr)
five_mpf_11X_pileup_df = transmute(five_mpf_11X_pileup_df, seqnames = seqnames, nucleotide = "C", start = start,
  context = "CG", dinucleotide = "CG", methylated_proportion = percent_modified/100, 
  methylated_counts = count_modified, total_counts = valid_coverage)
data.table::fwrite(five_mpf_11X_pileup_df, "cgmap_files/cgmap_files_nanopore/five_mpf_11X.CGmap.gz", 
  sep = "\t", col.names = F, quote = F)
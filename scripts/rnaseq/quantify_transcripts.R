# Quantify transcripts for Nematostella FASTQ files produced in house and downloaded from PRJNA189768 and PRJNA418421

# Load required packages
library(dplyr)
library(doParallel)
library(DESeq2)
library(methodical)
source("../auxillary_scripts/alignments_functions.R")

# Make a cluster
cl = makeCluster(6)
registerDoParallel(cl, cores = 6)

### Quantify transcripts for Nematostella FASTQ files generated in June 2024

# Get Paths to RNA-seq files
rnaseq_files = list.files(list.files("/data/home/btx717/storage/methylation_raw_data/NovoGene_NvecInhibitorATAC_RNA_EMseq_ControlPreplicatum_Cyanophora_2024_June/01.RawData", 
  full.names = T, pattern = "RNA"), full.names = T, pattern = "fq")
names(rnaseq_files) = gsub("/data/home/btx717/storage/methylation_raw_data/NovoGene_NvecInhibitorATAC_RNA_EMseq_ControlPreplicatum_Cyanophora_2024_June/01.RawData/", "", rnaseq_files)

# Create a data.frame with paths to files and their names
rnaseq_file_df = data.frame(filepaths = rnaseq_files, row.names = NULL)
rnaseq_file_df$sample = basename(dirname(rnaseq_file_df$filepaths))
rnaseq_file_df$file = gsub(".fq.gz", "", stringr::str_extract(rnaseq_file_df$filepaths, "\\d+\\.fq\\.gz"))

# For each sample, combine the two forward and two reverse read files into single files for forward and reverse reads
system.time(foreach(sample_name = unique(rnaseq_file_df$sample), .packages = "dplyr") %dopar% {
  
  sample_df = filter(rnaseq_file_df, sample == sample_name)
  forward_files = paste(filter(sample_df, file == "1")$filepaths, collapse = " ")
  reverse_files = paste(filter(sample_df, file == "2")$filepaths, collapse = " ")
  system(paste("zcat", forward_files, ">", paste0("rnaseq_june_2024/", sample_name, ".f1.fq")))
  system(paste("zcat", reverse_files, ">", paste0("rnaseq_june_2024/", sample_name, ".r2.fq")))
  
})

# Get paths to FASTQ files
gsk_fastq_files = list.files("rnaseq_june_2024", full.names = T)

# Get names of all forward FASTQ files
gsk_forward_fastqs = grep("f1", gsk_fastq_files, value = T)
gsk_reverse_fastqs = grep("r2", gsk_fastq_files, value = T)

# Name files with their SRR accession 
sample_names = gsub("_RNAseq", "", gsub("\\..*", "", basename(gsk_forward_fastqs)))

# Quantify transcripts with Kallisto. 
system.time(methodical::kallisto_quantify(path_to_kallisto = "~/programs/kallisto/kallisto", 
  kallisto_index = "../nematostella_genome/GCF__932526225.1_jaNemVect1.1_rna.kallisto.idx", 
  forward_fastqs = gsk_forward_fastqs, reverse_fastqs = gsk_reverse_fastqs, 
  sample_names = sample_names, 
  output_directory = "rnaseq_june_2024/kallisto_quantification", n_cores = 30))

### Quantify transcripts for Nematostella FASTQ files generated on the 20/03/2025

# Get Paths to RNA-seq files from 20-03-2025
rnaseq_files = list.files(list.files("/data/home/btx717/storage/methylation_raw_data/NovoGene_2025_March_Nvec_F1_and_Eggs/01.RawData", 
  full.names = T, pattern = "RNA"), full.names = T, pattern = "fq")
names(rnaseq_files) = gsub("/data/home/btx717/storage/methylation_raw_data/NovoGene_NvecInhibitorATAC_RNA_EMseq_ControlPreplicatum_Cyanophora_2024_June/01.RawData/", "", rnaseq_files)

# Create a data.frame with paths to files and their names
rnaseq_file_df = data.frame(filepaths = rnaseq_files, row.names = NULL)
rnaseq_file_df$sample = basename(dirname(rnaseq_file_df$filepaths))
rnaseq_file_df$file = gsub(".fq.gz", "", stringr::str_extract(rnaseq_file_df$filepaths, "\\d+\\.fq\\.gz"))

# For each sample, combine the two forward and two reverse read files into single files for forward and reverse reads. 
system.time(foreach(sample_name = unique(rnaseq_file_df$sample), .packages = "dplyr") %dopar% {
  
  sample_df = filter(rnaseq_file_df, sample == sample_name)
  forward_files = paste(filter(sample_df, file == "1")$filepaths, collapse = " ")
  reverse_files = paste(filter(sample_df, file == "2")$filepaths, collapse = " ")
  system(paste("zcat", forward_files, ">", paste0("rnaseq_march_2025/", sample_name, ".f1.fq")))
  system(paste("zcat", reverse_files, ">", paste0("rnaseq_march_2025/", sample_name, ".r2.fq")))
  
})

# Get paths to FASTQ files
combined_fastq_files = list.files("rnaseq_march_2025", full.names = T)

# Get names of all forward FASTQ files
forward_fastqs = grep("\\.f1", combined_fastq_files, value = T)
reverse_fastqs = grep("\\.r2", combined_fastq_files, value = T)

# Name files with their SRR accession 
sample_names = gsub("_RNAseq", "", gsub("\\..*", "", basename(forward_fastqs)))

# Quantify transcripts with Kallisto. 
system.time(methodical::kallisto_quantify(path_to_kallisto = "~/programs/kallisto/kallisto", 
  kallisto_index = "../nematostella_genome/GCF_932526225.1_jaNemVect1.1_rna_kallisto.idx", 
  forward_fastqs = forward_fastqs, reverse_fastqs = reverse_fastqs, 
  sample_names = sample_names, 
  output_directory = "rnaseq_march_2025/kallisto_quantification", n_cores = 30))

### Get FASTQ files for PRJNA189768

# Read in metadata table for PRJNA189768
PRJNA189768_metadata = data.table::fread("PRJNA189768_metadata.txt", data.table = F)

# Add column with number of hours post fertilization
PRJNA189768_metadata$timepoint = stringr::str_extract(PRJNA189768_metadata$experiment_desc, "\\d+\\s+(hours|days)")
PRJNA189768_metadata$number = as.numeric(gsub(" .*", "", PRJNA189768_metadata$timepoint))
PRJNA189768_metadata$timepoint = paste0(ifelse(grepl("days", PRJNA189768_metadata$timepoint), 
  PRJNA189768_metadata$number * 24, PRJNA189768_metadata$number), "Hpf")

# Add column with the replicate number and a column giving samples their final name
PRJNA189768_metadata$replicate = stringr::str_extract(PRJNA189768_metadata$experiment_desc, "replicate \\d+")
PRJNA189768_metadata$replicate = ifelse(PRJNA189768_metadata$replicate == "replicate 1", "r1", "r2")
PRJNA189768_metadata$sample_name = paste("C", PRJNA189768_metadata$timepoint, PRJNA189768_metadata$replicate, sep = "_")

# Prefetch files from PRJNA189768
sra_prefetch(path_to_sratk = "~/programs/sratoolkit.3.2.0-ubuntu64/bin", 
  srr_accessions = PRJNA189768_metadata$run_accession, 
  output_directory = "PRJNA189768_srr_files", parallel_files = 12)

# Extract FASTQ files from PRJNA189768_srr_files. 
system.time(sra_fastq_dump(path_to_sratk = "~/programs/sratoolkit.3.2.0-ubuntu64/bin", 
  srr_directory_list = list.files("PRJNA189768_srr_files", full.names = T), output_directory = "PRJNA189768_fastq_files", 
  parallel_files = 12))

# Remove PRJNA189768_srr_files
unlink("PRJNA189768_srr_files", recursive = T)

### Get FASTQ files for PRJNA418421

# Read in metadata table for PRJNA418421
PRJNA418421_metadata = data.table::fread("PRJNA418421_metadata.txt", data.table = F)

# Add a column with final sample names
PRJNA418421_metadata$sample_name = paste("R", 
  gsub("RNAseq of Nematostella vectensis embryo ", "", PRJNA418421_metadata$experiment_desc), sep = "_")
PRJNA418421_metadata$sample_name = gsub(" replicate ", "_r", PRJNA418421_metadata$sample_name)

# Prefetch files from PRJNA41842
system.time(sra_prefetch(path_to_sratk = "~/programs/sratoolkit.3.2.0-ubuntu64/bin", 
  srr_accessions = PRJNA418421_metadata$run_accession, 
  output_directory = "PRJNA418421_srr_files", parallel_files = 16))

# # Extract FASTQ files from PRJNA418421_srr_files. 
system.time(sra_fastq_dump(path_to_sratk = "~/programs/sratoolkit.3.2.0-ubuntu64/bin", 
  srr_directory_list = list.files("PRJNA418421_srr_files", full.names = T), output_directory = "PRJNA418421_fastq_files", 
  parallel_files = 16))

# Remove PRJNA189768_srr_files
unlink("PRJNA418421_srr_files", recursive = T)

### Quantify transcripts with Kallisto for PRJNA189768 and PRJNA418421

# Combine metadata for PRJNA418421 and PRJNA189768
combined_metadata = dplyr::bind_rows(PRJNA418421_metadata, PRJNA189768_metadata)

# Create a vector matching run accessions to sample names
accession_to_sample_name_dict = setNames(combined_metadata$sample_name, combined_metadata$run_accession)

# Get paths to FASTQ files
fastq_files = c(list.files("PRJNA189768_fastq_files", full.names = T), 
  list.files("PRJNA418421_fastq_files", full.names = T))

# Name files with their SRR accession 
names(fastq_files) = gsub(".fastq.gz", "", basename(fastq_files))

# Quantify transcripts with Kallisto. 
system.time(methodical::kallisto_quantify(path_to_kallisto = "~/programs/kallisto/kallisto", 
  kallisto_index = "../nematostella_genome/GCF__932526225.1_jaNemVect1.1_rna.kallisto.idx", 
  forward_fastqs = fastq_files, reverse_fastqs = NULL, 
  sample_names = accession_to_sample_name_dict[names(fastq_files)], 
  output_directory = "external_rnaseq_kallisto_quantification", n_cores = 30, messages_file = "", number_bootstraps = 100))

# Load transcript counts for external data
external_transcript_counts = round(data.frame(data.table::fread("external_rnaseq_kallisto_quantification/kallisto_transcript_counts_merged.tsv.gz", data.table = F), row.names = 1))

# Create a vector with the timepoints for each sample and a list with the samples for each timepoint
sample_timepoints = gsub(".*_", "", gsub("Hpf.*", "", names(external_transcript_counts)))
sample_timepoints_list = split(names(external_transcript_counts), sample_timepoints)
sample_timepoints_list = sample_timepoints_list[as.character(sort(as.numeric(names(sample_timepoints_list))))]
names(sample_timepoints_list) = paste0("Hpf_", names(sample_timepoints_list))

# Get the mean gene counts for the sample at each timepoint and round to integers
combined_timepoint_counts = data.frame(lapply(sample_timepoints_list, function(x) 
  round(rowMeans(external_transcript_counts[, x]))))

# Load a list matching genes to transcripts
genes_to_transcript_list = readRDS("../nematostella_genome/genes_to_transcript_list.rds")

# Combine transcript counts for genes
system.time({combined_timepoint_gene_counts = methodical::sumTranscriptValuesForGenes(combined_timepoint_counts, genes_to_transcript_list)})
data.table::fwrite(tibble::rownames_to_column(combined_timepoint_gene_counts, "gene_id"), 
  "combined_timepoint_gene_counts.tsv.gz", sep = "\t")

# Normalize combined_timepoint_gene_counts using DESeq2
combined_timepoint_gene_counts_dds = DESeqDataSetFromMatrix(countData = combined_timepoint_gene_counts, 
  colData = data.frame(sample = names(combined_timepoint_gene_counts)), design = ~ 1)
combined_timepoint_gene_counts_dds  = estimateSizeFactors(combined_timepoint_gene_counts_dds) 
combined_timepoint_gene_counts_normalized = data.frame(counts(combined_timepoint_gene_counts_dds, normalized = T))

# Max normalize each gene
combined_timepoint_gene_counts_max_normalized = data.frame(t(apply(combined_timepoint_gene_counts_normalized, 1, function(x) x/max(x))))

# Find the timepoint of maximum expression for each gene
max_timepoint_per_gene = apply(combined_timepoint_gene_counts_normalized, 1, function(x) names(which.max(x)))
max_timepoint_per_gene = factor(max_timepoint_per_gene, levels = names(combined_timepoint_gene_counts_normalized))
saveRDS(max_timepoint_per_gene, "max_timepoint_per_gene.rds")

### Create a table with normalized counts for our Nematostella samples

# Load tables with Kallisto counts
june_2024_transcript_counts = round(data.frame(data.table::fread("rnaseq_june_2024/kallisto_quantification/kallisto_transcript_counts_merged.tsv.gz", data.table = F), row.names = 1))
march_2025_transcript_counts = round(data.frame(data.table::fread("rnaseq_march_2025/kallisto_quantification/kallisto_transcript_counts_merged.tsv.gz", data.table = F), row.names = 1))[c(3, 4, 1, 2)]

# Combine tables and update names
complete_transcript_counts = cbind(june_2024_transcript_counts, march_2025_transcript_counts)
names(complete_transcript_counts) = gsub("Nvec_24hpf_", "", names(complete_transcript_counts))
names(complete_transcript_counts) = gsub("_25hpf", "", names(complete_transcript_counts))
data.table::fwrite(tibble::rownames_to_column(complete_transcript_counts, "transcript_id"), 
  "nvec_complete_transcript_counts_raw.tsv.gz", sep = "\t")

# Load a list matching genes to transcripts
genes_to_transcript_list = readRDS("../nematostella_genome/genes_to_transcript_list.rds")

# Combine transcript counts for the same gene. 
system.time({nvec_complete_gene_counts_raw = methodical::sumTranscriptValuesForGenes(complete_transcript_counts, genes_to_transcript_list)})
data.table::fwrite(tibble::rownames_to_column(nvec_complete_gene_counts_raw, "gene_id"), 
  "nvec_complete_gene_counts_raw.tsv.gz", sep = "\t")

# Normalize nvec_gene_counts using DESeq2
complete_transcript_counts_dds = DESeqDataSetFromMatrix(countData = complete_transcript_counts, 
  colData = data.frame(row.names = names(complete_transcript_counts)), design = ~ 1)
complete_transcript_counts_dds  = estimateSizeFactors(complete_transcript_counts_dds) 
complete_transcript_counts_normalized = data.frame(counts(complete_transcript_counts_dds, normalized = T))

# Save table
data.table::fwrite(tibble::rownames_to_column(complete_transcript_counts_normalized, "transcript_id"), 
  "nvec_complete_transcript_counts_normalized.tsv.gz", sep = "\t")
nvec_complete_transcript_counts_normalized = data.frame(data.table::fread("nvec_complete_transcript_counts_normalized.tsv.gz"), row.names = 1)

# Perform differential expression

# Subset nvec_complete_gene_counts_normalized for DMSO and GSK samples
nvec_complete_gene_counts_raw = data.frame(data.table::fread("nvec_complete_gene_counts_raw.tsv.gz"), row.names = 1)
dmso_and_gsk_gene_counts = nvec_complete_gene_counts_raw[grep("DMSO|GSK", names(nvec_complete_gene_counts_raw), value = T)]

# Create a metadata table for dmso_and_gsk_gene_counts
de_coldata = data.frame(sample = names(dmso_and_gsk_gene_counts), 
  treatment = c(rep("DMSO", 3), rep("GSK", 3)))

# Test differential expression with DESeq2
library(DESeq2)
nvec_gene_counts_dds = DESeqDataSetFromMatrix(countData = dmso_and_gsk_gene_counts, colData = de_coldata, design = ~ treatment)
nvec_gene_counts_dds  = estimateSizeFactors(nvec_gene_counts_dds) 
system.time({nvec_gene_counts_dds = DESeq(nvec_gene_counts_dds)})
deseq_results_gsk_vs_dmso = data.frame(results(nvec_gene_counts_dds, contrast = c("treatment", "GSK", "DMSO"), alpha = 0.05))
deseq_results_gsk_vs_dmso_sig = filter(deseq_results_gsk_vs_dmso, padj < 0.05)
data.table::fwrite(tibble::rownames_to_column(deseq_results_gsk_vs_dmso_sig, "gene_id"),
  "deseq_gsk_vs_dmso_significant_results.tsv.gz", sep = "\t")
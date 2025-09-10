# Create BSSeq and meth RSE objects for Nematostella samples

# Load required packages
library(dplyr)
library(methrix)
library(BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19)

# Read file with description of EMseq and Nanopore samples
sample_desciption = filter(data.table::fread("nvec_sample_description.tsv", sep = "\t"), `EM-seq` | `ONT`)
names(sample_desciption) = gsub(" ", "_", names(sample_desciption))
sample_desciption$Sample_Name = gsub("-", "", sample_desciption$Sample_Name)

# Update names of samples in sample_desciption which were incorrectly named
sample_desciption$Sample_Name[19:24] = c(paste0("GSK_rep", 1:3), paste0("DMSO_rep", 1:3))

# Create colData for EMseq samples
sample_coldata = data.frame(
  description = sample_desciption$Description, 
  timepoint = sample_desciption$Sampling_time,
  method = ifelse(sample_desciption$`EM-seq`, "EMseq", "Nanopore"),
  row.names = sample_desciption$Sample_Name
)

# Create metadata for skimming samples and combine with sample_coldata
skimming_samples_sample_coldata = data.frame(
  description = c("GSK F0 pooled skimming samples 3.6%", "GSK F0 pooled skimming samples 3.1%"),
  timepoint = c("F0 5mpf", "F0 3mpf"),
  method = "Nanopore",
  row.names = c("skimming_five_mpf_11X", "skimming_three_mpf_15X")
)
sample_coldata = rbind(sample_coldata, skimming_samples_sample_coldata)

# Get paths to EMseq BSseeker CGmap files
cgmap_files_emseq = list.files("cgmap_files/cgmap_files_emseq", full.names = T, pattern = "CGmap.gz")

# Get paths to Nanopore CGmap files
cgmap_files_nanopore = list.files("cgmap_files/cgmap_files_nanopore", full.names = T, pattern = "CGmap.gz")

# Get paths to Nanopore skimming CGmap files and combine with cgmap_files_nanopore
cgmap_files_nanopore_skimming = list.files("cgmap_files/cgmap_files_nanopore_skimming", full.names = T, pattern = "CGmap.gz")
cgmap_files_nanopore = c(cgmap_files_nanopore, cgmap_files_nanopore_skimming)

# Filter CGmap files for NC contigs
filter_nc = function(cgmap_file){
  nc_file = gsub(".CGmap.gz", "_nc_contigs.CGmap.gz", cgmap_file)
  system(paste("zcat", cgmap_file, "| grep NC_ | gzip >", nc_file ))
}
lapply(c(cgmap_files_emseq, cgmap_files_nanopore), filter_nc)

# Get paths to EMseq files for NC contigs and name them
cgmap_files_emseq_nc = list.files("cgmap_files/cgmap_files_emseq", full.names = T, pattern = "nc_contigs")
names(cgmap_files_emseq_nc) = tools::file_path_sans_ext(gsub("EMseq_", "", gsub("_nc_contigs.CGmap", "", basename(cgmap_files_emseq_nc))))

# Get paths to Nanopore files for NC contigs and name them
cgmap_files_nanopore_nc = list.files(c("cgmap_files/cgmap_files_nanopore", "cgmap_files/cgmap_files_nanopore_skimming"), full.names = T, pattern = "nc_contigs")
names(cgmap_files_nanopore_nc) = tools::file_path_sans_ext(gsub("_nc_contigs.CGmap", "", basename(cgmap_files_nanopore_nc)))

# Combine emseq and nanopore files
combined_files = c(cgmap_files_emseq_nc, cgmap_files_nanopore_nc)

# Put sample_annotation in the same order as combine_files
sample_coldata = sample_coldata[names(combined_files), ]

# Get reference CpGs from Nematostella genome. 
source("../auxillary_scripts/extract_CPGs_modified.R")
system.time({nc_cpg_sites = extract_CPGs_modified("BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19", 
  seqnames = grep("NC", seqnames(BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19), value = T))})
system.time({chr_cpg_sites = extract_CPGs_modified("BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19", 
   seqnames = grep("chr", seqnames(BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19), value = T))})

# Create a subset of contigs
nc_contigs = grep("NC", seqnames(BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19), value = T)
chr_contigs = grep("chr", seqnames(BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19), value = T)

# Create a methrix object for the NC contigs. 
system.time({nc_contigs_methrix = read_bedgraphs(files = combined_files, 
    chr_idx = 1, start_idx = 3, beta_idx = 6, cov_idx = 8, contigs = nc_contigs, zero_based = F,
    stranded = TRUE, collapse_strands = T, ref_cpgs = nc_cpg_sites, coldata = sample_coldata, 
    h5 = T, h5_dir = "nc_contigs_methrix")})
nc_contigs_methrix = methrix::load_HDF5_methrix("nc_contigs_methrix")

# Create a methrix object for the chr contigs. 
system.time({chr_contigs_methrix = read_bedgraphs(files = combined_files, 
  chr_idx = 1, start_idx = 3, beta_idx = 6, cov_idx = 8, contigs = chr_contigs, zero_based = F,
  stranded = TRUE, collapse_strands = T, ref_cpgs = chr_cpg_sites, coldata = sample_coldata, 
  h5 = T, h5_dir = "chr_contigs_methrix")})
chr_contigs_methrix = methrix::load_HDF5_methrix("chr_contigs_methrix")

# Combine methrix objects and save. 
complete_methrix = c(nc_contigs_methrix, chr_contigs_methrix)
system.time({save_HDF5_methrix(complete_methrix, dir = "nematostella_complete_methrix")})

# Convert methrix to a BSSeq object. 
system.time({nematostella_complete_bsseq = methrix::methrix2bsseq(complete_methrix)})
system.time({HDF5Array::saveHDF5SummarizedExperiment(nematostella_complete_bsseq, "nematostella_complete_bsseq")})

# Convert methrix to a meth RSE object. 
system.time({nematostella_complete_meth_rse = methodical::methrixToRSE(complete_methrix)})
system.time({HDF5Array::saveHDF5SummarizedExperiment(nematostella_complete_meth_rse, "nematostella_complete_meth_rse")})

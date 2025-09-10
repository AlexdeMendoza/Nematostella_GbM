# Create a BSgenome object from Nematostella_DToL_lambda_pUC_mitochondria_originalIDs.fasta and
# create a Kallisto index from the transcripts FASTA file

# Load required packages
library(methodical)
library(Biostrings)

# Read in Nematostella_DToL_lambda_pUC_mitochondria_originalIDs.fasta and shorten chromosome names to just the sequence name
nematostella_fasta = readDNAStringSet("Nematostella_DToL_lambda_pUC_mitochondria_originalIDs.fasta")
names(nematostella_fasta) = gsub(" .*", "", names(nematostella_fasta))

# Use seqkit to split nematostella_standard_chroms.fasta into separate FASTA files
system("~/programs/seqkit split -i Nematostella_DToL_lambda_pUC_mitochondria_originalIDs.fasta -O nematostella_fastas")

# Create a BSgenome package using nematostella_bsgenome_seed.txt and nematostella_fastas/
system.time({BSgenome::forgeBSgenomeDataPkg(x = "nematostella_bsgenome_seed.txt", 
  seqs_srcdir = "nematostella_fastas", destdir = ".")})

# Create a Kallisto index from GCF_932526225.1_jaNemVect1.1_rna.fna.gz. 
# GCF_932526225.1_jaNemVect1.1_rna.fna.gz downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/932/526/225/GCF_932526225.1_jaNemVect1.1/GCF_932526225.1_jaNemVect1.1_rna.fna.gz
system.time({methodical::kallistoIndex(path_to_kallisto = "~/programs/kallisto/kallisto",  
  transcripts_fasta = "GCF_932526225.1_jaNemVect1.1_rna.fna.gz", 
  index_name = "GCF_932526225.1_jaNemVect1.1_rna_kallisto.idx")})
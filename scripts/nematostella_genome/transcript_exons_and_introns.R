# Create GRangesList objects with exons and introns for transcripts from GCF_932526225.1_jaNemVect1.1_genomic.gff.gz

# Load required packages
library(GenomicFeatures)
library(dplyr)

# Read in GCF_932526225.1_jaNemVect1.1_genomic.gff.gz.gz as a GRanges 
system.time({nvec_gff_gr = rtracklayer::import.gff3("GCF_932526225.1_jaNemVect1.1_genomic.gff.gz")})

# Drop unnecessary metadata columns 
mcols(nvec_gff_gr) = mcols(nvec_gff_gr)[c("type", "ID", "gene", "gene_biotype", "transcript_id", "product")]

# Make GRanges object for genes 
gene_gr = nvec_gff_gr[nvec_gff_gr$type == "gene"]
names(gene_gr) = gene_gr$gene

# Remove transcript_id and product metadata columns
mcols(gene_gr)[c("transcript_id", "product")] = NULL

# Subset for protein-coding genes and save and export as a BED file
protein_gene_gr = gene_gr[gene_gr$gene_biotype == "protein_coding"]
saveRDS(protein_gene_gr, "nvec_pc_genes_gr.rds")
rtracklayer::export.bed(protein_gene_gr, "nvec_pc_genes.bed")

# Make GRanges object for transcripts. Note that transcript_id and not ID should be used to name transcripts
transcripts_gr = nvec_gff_gr[nvec_gff_gr$type %in% c("mRNA", "transcript")]
names(transcripts_gr) = transcripts_gr$transcript_id

# Remove gene_biotype metadata column
transcripts_gr$gene_biotype = NULL

# Make a list matching genes to their gene type
genes_to_type_list = setNames(gene_gr$gene_biotype, gene_gr$gene)

# Make a data.frame matching genes to transcript IDs
gene_to_transcript_df = data.frame(gene = transcripts_gr$gene, transcript = transcripts_gr$transcript_id)
gene_to_transcript_df$gene_type = genes_to_type_list[gene_to_transcript_df$gene]

# Get a list of protein coding transcripts
protein_transcripts = filter(gene_to_transcript_df, gene_type == "protein_coding")$transcript

# Create GRanges objects for protein-coding transcripts
pc_transcripts_gr = transcripts_gr[transcripts_gr$transcript_id %in% protein_transcripts]
saveRDS(pc_transcripts_gr, "nvec_pc_transcripts_gr.rds")

# Get TSS for protein-coding transcripts and save
pc_transcripts_tss_gr = resize(pc_transcripts_gr, width = 1, fix = "start")
saveRDS(pc_transcripts_tss_gr, "nvec_pc_transcripts_tss_gr.rds")

# Create a list matching gene names to transcripts and save
genes_to_transcript_list = split(pc_transcripts_gr$transcript_id, pc_transcripts_gr$gene)
saveRDS(genes_to_transcript_list, "genes_to_transcript_list.rds")

### Create a GRanges with promoters, exons and introns for protein-coding transcripts

# Create a GRanges for promoters, defined as the 1 kb upstream of the TSS
promoters_gr = promoters(pc_transcripts_tss_gr, upstream = 1000)
mcols(promoters_gr) = data.frame(transcript_id = promoters_gr$transcript_id, region = "promoter", exon_rank = 0)
names(promoters_gr) = promoters_gr$transcript_id

# Create a TxDb object from GCF_932526225.1_jaNemVect1.1_genomic.gff.gz 
system.time({nvec_txdb = makeTxDbFromGFF("GCF_932526225.1_jaNemVect1.1_genomic.gff.gz")})

# Fetch exons and introns for each protein-coding transcript in nvec_txdb 
exons_gr = unlist(exonsBy(nvec_txdb, "tx", use.names = T)[protein_transcripts])
exons_gr$region = "exon"
introns_gr = unlist(intronsByTranscript(nvec_txdb, use.names = T)[protein_transcripts])
introns_gr$region = "intron"

# Combine exons and introns into a single GRanges object
promoters_exons_and_introns_gr = c(promoters_gr, exons_gr, introns_gr)

# Add transcript name as a column
promoters_exons_and_introns_gr$transcript_id = names(promoters_exons_and_introns_gr)

# Convert exons_and_introns into a data.frame
promoters_exons_and_introns_df = data.frame(promoters_exons_and_introns_gr)

# Group data.frame by transcript_id
promoters_exons_and_introns_df = group_by(promoters_exons_and_introns_df, transcript_id)

# Arrange promoters_exons_and_introns_df by transcript_id and then by start. 
system.time({promoters_exons_and_introns_df = arrange(promoters_exons_and_introns_df, transcript_id, start)})

# Rename exon_rank to just rank
promoters_exons_and_introns_df = dplyr::rename(promoters_exons_and_introns_df, "rank" = "exon_rank")

# Add rank to introns based on the preceding exon and whether they are on the "+" or negative "-" strand
promoters_exons_and_introns_df$rank[which(is.na(promoters_exons_and_introns_df$rank) & promoters_exons_and_introns_df$strand == "+")] =  
  promoters_exons_and_introns_df$rank[which(is.na(promoters_exons_and_introns_df$rank) & promoters_exons_and_introns_df$strand == "+") - 1]
promoters_exons_and_introns_df$rank[which(is.na(promoters_exons_and_introns_df$rank) & promoters_exons_and_introns_df$strand == "-")] =  
  promoters_exons_and_introns_df$rank[which(is.na(promoters_exons_and_introns_df$rank) & promoters_exons_and_introns_df$strand == "-") + 1]

# Add gene ID
promoters_exons_and_introns_df$gene = nvec_transcripts_gr$gene[match(promoters_exons_and_introns_df$transcript_id, nvec_transcripts_gr$transcript_id)]

# Remove unneccessary columns
promoters_exons_and_introns_df = dplyr::select(promoters_exons_and_introns_df, seqnames, start, end, width, strand, gene, transcript_id = transcript_id, region, rank)

# Recreate a GRanges from the data.frame
promoters_exons_and_introns_gr = makeGRangesFromDataFrame(promoters_exons_and_introns_df, keep.extra.columns = T)

# Name regions with transcript ID, exon/intron status and rank
names(promoters_exons_and_introns_gr) = paste(promoters_exons_and_introns_gr$transcript_id, promoters_exons_and_introns_gr$region, promoters_exons_and_introns_gr$rank, sep = "_")
saveRDS(promoters_exons_and_introns_gr, "nvec_pc_transcripts_promoters_exons_and_introns_gr.rds")

### Create a GRanges with the promoter, exons and introns for the longestg transcript for each gene

# Convert pc_transcripts_gr into a data.frame and group by gene
pc_transcripts_df = group_by(data.frame(pc_transcripts_gr), gene)

# Find the longest transcript for each gene, selecting the first transcript by alphabetical order for ties
pc_transcripts_df = filter(pc_transcripts_df, width == max(width))
longest_transcripts = filter(pc_transcripts_df, transcript_id == min(transcript_id))$transcript_id

# Filter promoters_exons_and_introns_gr for the longest transcript for each gene
longest_transcripts_promoters_exons_and_introns_gr = promoters_exons_and_introns_gr[promoters_exons_and_introns_gr$transcript_id %in% longest_transcripts]
saveRDS(longest_transcripts_promoters_exons_and_introns_gr, "longest_transcripts_promoters_exons_and_introns_gr.rds")

### Create GRanges objects for other types of regions besides genes

# Define types of non-coding regions of interest
selected_regions = c("pseudogene", "snRNA", "lnc_RNA", "rRNA", "tRNA")

# Create a GRanges object with these regions. Pseudogenes and tRNAs do not have any transcript_id
non_coding_regions_gr = unlist(GRangesList(lapply(selected_regions, function(x) 
  nvec_gff_gr[nvec_gff_gr$type == x])))

# Delete gene_biotype metadata column as it is missing in most regions
non_coding_regions_gr$gene_biotype = NULL

# Sort regions and save
non_coding_regions_gr = sort(non_coding_regions_gr)
saveRDS(non_coding_regions_gr, "nvec_non_coding_regions_gr.rds")

### Create an SAF file with genes as features for use with featureCounts with the exons from GCF_932526225.1_jaNemVect1.1_genomic.gff.gz

# Make GRanges object for exons and remove all metadata columns except gene
exons_gr = nvec_gff_gr[nvec_gff_gr$type == "exon"]
mcols(exons_gr) = data.frame("gene" = exons_gr$gene)

# Convert exons_gr into a data.frame, put columns in correct order and rename them
genes_saf_df = dplyr::select(data.frame(exons_gr), 
  "Geneid" = gene, "Chr" = seqnames, "Start" = start, "End" = end, "Strand" = strand)
data.table::fwrite(genes_saf_df, "nvec_genes.saf.gz", sep = "\t")
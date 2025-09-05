# # Create bigWigs with chromosome names in UCSC format for mouse WGBS and H3K79me2/3 data

# Load required packages
library(rtracklayer)
library(plotR)

### Process WGBS bigWig files

# Get names of WGBS bedGraph files
wgbs_bg_files = list.files("mouse_wgbs/bedgraphs", full.names = T)

# Convert chromosome names to UCSC names
system.time({convert_chrom_names_to_ucsc(chrom_to_ucsc_path = "~/programs/ucsc_tools/chromToUcsc", chrom_alias_file = "../mm10.chromAlias.tsv", 
  input_files = wgbs_bg_files, output_dir = "mouse_wgbs/bedgraphs", ncores = 10)})

# Get paths to bedGraphs with UCSC chromosome names
ucsc_wgbs_bg_files = list.files("mouse_wgbs/bedgraphs", pattern = "UCSC", full.names = T)

# Convert bedGraphs to bigWigs
system.time({convert_bedgraphs_to_bigwigs(bedgraph_to_bigwig_path = "~/programs/ucsc_tools/bedGraphToBigWig", ucsc_wgbs_bg_files, 
  chrom_sizes_file = "mm10.chrom.sizes", output_dir = "mouse_wgbs/bigwigs", ncores = 10)})

### Process H3K79me2/3 bedGraph files

# Get Paths to H3K79 bigWigs
h3k79_bigwigs = list.files(path = "h3k79/bigwigs/", full.names = T)

# Convert bigWigs to bedGRaphs
convert_bigwigs_to_bedgraphs(bigwig_to_bedgraph_path = "~/programs/ucsc_tools/bigWigToBedGraph", 
  bw_files = h3k79_bigwigs, compress_output_bedgraphs = F, output_dir = "h3k79/bedgraphs", ncores = 6)

# Get paths to all bedGraphs
h3k79_bedgraphs = list.files(path = "h3k79/bedgraphs/", full.names = T)

# Convert chromosome names to UCSC names
system.time({convert_chrom_names_to_ucsc(chrom_to_ucsc_path = "~/programs/ucsc_tools/chromToUcsc", chrom_alias_file = "mm10.chromAlias.tsv", 
  input_files = h3k79_bedgraphs, output_dir = "h3k79/bedgraphs/", ncores = 1)})

# Get paths to bedGraphs with UCSC chromosome names
ucsc_h3k79_bedgraphs = list.files("h3k79/bedgraphs/", pattern = "UCSC", full.names = T)

# Convert bedGraphs to bigWigs
system.time({convert_bedgraphs_to_bigwigs(bedgraph_to_bigwig_path = "~/programs/ucsc_tools/bedGraphToBigWig", ucsc_h3k79_bedgraphs, 
  chrom_sizes_file = "mm10.chrom.sizes", output_dir = "h3k79/bigwigs", ncores = 6)})

# Get mouse GTF as GRanges
mouse_gtf = rtracklayer::import.gff2("gencode.vM10.annotation.gtf.gz")

# Filter for genes and export as a BED
mouse_gene_gr = mouse_gtf[mouse_gtf$type == "gene" & mouse_gtf$gene_type == "protein_coding"]
mcols(mouse_gene_gr) = NULL
rtracklayer::export.bed(mouse_gene_gr, "mouse_genes.bed")

# Get paths to all bigWigs and set names
wgbs_bws = list.files("mouse_wgbs/bigwigs", full.names = T)
names(wgbs_bws) = gsub("_UCSC.*", "", gsub(".*_WT", "WGBS", wgbs_bws))
h3k79_bws = list.files("h3k79/bigwigs", full.names = T, pattern = "UCSC")
names(h3k79_bws) = gsub("_R.*", "", gsub(".*_WT-", "", h3k79_bws))
normalized_h3k79_bws = list.files("h3k79/bigwigs", full.names = T, pattern = "input*.bw")
names(normalized_h3k79_bws) = basename(normalized_h3k79_bws)

# Get gene values for WGBS and H3K79 
system.time({gene_wgbs_values = bigwig_summarize_over_regions(bw_filepaths = wgbs_bws, gr = mouse_gene_gr, 
  statistic = "mean", column_names = names(wgbs_bws), ncores = 5)})
system.time({gene_h3k79_values = bigwig_summarize_over_regions(bw_filepaths = h3k79_bws, gr = mouse_gene_gr, 
  statistic = "mean0", column_names = names(h3k79_bws), ncores = 6)})
system.time({gene_h3k79_values_normalized = bigwig_summarize_over_regions(bw_filepaths = normalized_h3k79_bws, gr = mouse_gene_gr, 
  statistic = "mean0", column_names = names(normalized_h3k79_bws), ncores = 6)})

# Combine two tables
gene_wgbs_and_h3k79_values = cbind(gene_wgbs_values, gene_h3k79_values_normalized)

# Calculate correlations between features
gene_wgbs_and_h3k79_cors = cor(gene_wgbs_and_h3k79_values, method = "s")

gene_wgbs_and_h3k79_cor_heatmap = plotR::heatmap_without_clustering(gene_wgbs_and_h3k79_cors, number_size = 10, label_size = 10)
gene_wgbs_and_h3k79_cor_heatmap
ggsave(plot = gene_wgbs_and_h3k79_cor_heatmap, "mouse_gene_wgbs_and_h3k79_cor_heatmap.pdf", width = 9, height = 9)

# Plot 
plot(log(gene_wgbs_and_h3k79_values$h3k79me2_sample1_log2ratio_vs_input.bw), gene_wgbs_and_h3k79_values$WGBS_DMSO_d14)

#
selected_bigwigs_names = gsub("_R.*", "", gsub(".*_WT.", "", gsub("UCSC", "R2", basename(readLines("selected_bigwigs.txt")))))
cat(selected_bigwigs_names, file = "selected_bigwigs_names.txt", sep = "\n")

# Create a GRanges object with Nematostella repeats from GCF_932526225.RM.bed

# Load required packages
library(GenomicRanges)

# Load repeats for Nematostella
nvec_repeats = data.table::fread("GCF_932526225.RM.bed")
names(nvec_repeats) = c("chrom", "start", "end", "name", "alt_name", "strand")

# Create columns with family and subfamily
nvec_repeats$superfamily = gsub("/.*", "", nvec_repeats$name)
nvec_repeats$family = gsub(".*/", "", nvec_repeats$name)

# Create a GRanges object from nvec_repeats
nvec_repeats_gr = makeGRangesFromDataFrame(nvec_repeats, seqnames.field = "chrom", start.field = "start", 
  end.field = "end", strand.field = "strand", keep.extra.columns = T, starts.in.df.are.0based = T)
saveRDS(nvec_repeats_gr, "nvec_repeats_gr.rds")

# Load table for Kimura score for repeats and add repeat class and subclass
kimura_bed = data.table::fread("Nematostella_vectensis.fasta.align.kimura.bed")
names(kimura_bed) = c("chr", "start", "end", "name", "score")
kimura_bed$class = gsub("/.*", "", gsub(".*#", "", kimura_bed$name))
kimura_bed$subclass = gsub(".*#", "", gsub(".*/", "", kimura_bed$name))

# Convert kimura_bed to a GRanges and save
kimura_bed_gr = makeGRangesFromDataFrame(kimura_bed, keep.extra.columns = T)

# Get the median Kimura score for each subclass of repeats
repeat_subclass_median_kimura_score = sapply(split(kimura_bed$score, kimura_bed$subclass), median)
saveRDS(repeat_subclass_median_kimura_score, "repeat_subclass_median_kimura_score.rds")

### Create a GRanges with complete genomic annotation for Nematostella

# Load GRanges for repeats and use alt_name as transcript_id
nvec_repeats_gr = readRDS("nvec_repeats_gr.rds")

# Rename DNA superfamily to DNA Transposon and Simple_repeat superfamily to Simple Repeat
nvec_repeats_gr$superfamily[nvec_repeats_gr$superfamily == "Simple_repeat"] = "Simple Repeat"
nvec_repeats_gr$superfamily[nvec_repeats_gr$superfamily == "DNA"] = "DNA Transposon"

# Set superclass as Repetitive DNA, use superfamily for class and family for subclass
nvec_repeats_gr$superclass = "Repetitive DNA"
nvec_repeats_gr$class = nvec_repeats_gr$superfamily
nvec_repeats_gr$subclass = nvec_repeats_gr$family
nvec_repeats_gr$transcript_id = nvec_repeats_gr$alt_name

# Load GRanges with non-coding regions 
non_coding_regions_gr = readRDS("../nvec_non_coding_regions_gr.rds")

# Rename some of the types of non-coding regions in non_coding_regions_gr
non_coding_regions_gr$type = as.character(non_coding_regions_gr$type)
non_coding_regions_gr$type[non_coding_regions_gr$type == "lnc_RNA"] = "lncRNA"
non_coding_regions_gr$type[non_coding_regions_gr$type == "pseudogene"] = "Pseudogene"

# Use type for superclass, class and subclass
non_coding_regions_gr$superclass = non_coding_regions_gr$type
non_coding_regions_gr$class = non_coding_regions_gr$type
non_coding_regions_gr$subclass = non_coding_regions_gr$type

# Load GRanges with promoters, exons and introns
promoters_exons_and_introns_gr = readRDS("../nvec_pc_transcripts_promoters_exons_and_introns_gr.rds")

# Rename the regions in promoters_exons_and_introns_gr
promoters_exons_and_introns_gr$region[promoters_exons_and_introns_gr$region == "promoter"] = "Promoter"
promoters_exons_and_introns_gr$region[promoters_exons_and_introns_gr$region == "exon"] = "Exon"
promoters_exons_and_introns_gr$region[promoters_exons_and_introns_gr$region == "intron"] = "Intron"

# Set class as pc_gene body for exons and introns and pc_promoter otherwise and separate exons and introns for subclass
promoters_exons_and_introns_gr$superclass = ifelse(promoters_exons_and_introns_gr$region == "Promoter", "Promoter", "Gene Body")
promoters_exons_and_introns_gr$class =  ifelse(promoters_exons_and_introns_gr$region == "Promoter", 
  "Promoter", promoters_exons_and_introns_gr$region)
promoters_exons_and_introns_gr$subclass = promoters_exons_and_introns_gr$class

# Combine nvec_repeats_gr and promoters_exons_and_introns_gr
complete_annotation_gr = c(promoters_exons_and_introns_gr, non_coding_regions_gr, nvec_repeats_gr)

# Select necessary metadata columns and remove unused levels from class
mcols(complete_annotation_gr) = mcols(complete_annotation_gr)[c("superclass", "class", "subclass", "transcript_id")]
complete_annotation_gr$superclass = as.character(complete_annotation_gr$superclass)
complete_annotation_gr$class = as.character(complete_annotation_gr$class)

# Remove non-standard chromosomes
seqlevels(complete_annotation_gr, pruning.mode = "coarse") = seqlevels(complete_annotation_gr)[1:15]

# Sort complete_annotation_gr and save
complete_annotation_gr = sort(complete_annotation_gr, ignore.strand = T)
names(complete_annotation_gr) = NULL
saveRDS(complete_annotation_gr, "nvec_complete_annotation_gr.rds")

# For plotting convert the gene body superclass into introns and exons
complete_annotation_gr$superclass[complete_annotation_gr$superclass == "Gene Body"] = complete_annotation_gr$class[complete_annotation_gr$superclass == "Gene Body"] 

# Count the number of elements for each superclass
superclass_count = table(complete_annotation_gr$superclass)

# Convert complete_annotation_gr into a list
complete_annotation_grl = split(complete_annotation_gr, complete_annotation_gr$superclass)

# Get the length of the Nematostella genome. Is about 269 MB
nvec_geneome_length = sum(lengths(BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19::BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19)[1:15])

# Get proportion of Nematostella genome covered by each superclass of genomic region
genomic_region_class_proportions = sort(sapply(complete_annotation_grl, function(x) sum(width(reduce(x, ignore.strand = T)))))/nvec_geneome_length
genomic_region_class_proportions = data.frame(
  region = names(genomic_region_class_proportions),
  proportion = genomic_region_class_proportions, row.names = NULL
)
genomic_region_class_proportions = filter(genomic_region_class_proportions, region != "ARTEFACT")
genomic_region_class_proportions = arrange(genomic_region_class_proportions, proportion)
genomic_region_class_proportions$region = 
  factor(genomic_region_class_proportions$region, levels = genomic_region_class_proportions$region)
x_labels = as.character(genomic_region_class_proportions$region)

# Append number of repeats to each class
x_labels = paste(x_labels, paste("n =", prettyNum(superclass_count[x_labels], big.mark = ",")), sep = "\n")

# Create a barplot of the genomic proportion of different classes of genomic regions
genomic_region_class_proportions_plot = ggplot(genomic_region_class_proportions, aes(x = region, y = proportion, fill = region)) +
  geom_col(color="black", fill = colorRampPalette(RColorBrewer::brewer.pal(n = 20, name = "RdBu"))(10)) +
  scale_fill_brewer(palette="Set1")
genomic_region_class_proportions_plot = customize_ggplot_theme(genomic_region_class_proportions_plot, 
  xlab = "Genomic Region Class", ylab = "Proportion", title = "Proportion of Nematostella Genome Covered\nby Different Classes of Genomic Region",
  x_labels = x_labels) +
  coord_flip()
genomic_region_class_proportions_plot
ggsave(plot = genomic_region_class_proportions_plot, "genomic_region_class_proportions_barplot.pdf", width = 9, height = 9)

# Make a plot of the number of repeats

# Convert repeat_annotation_gr into a list
repeat_annotation_grl = split(nvec_repeats_gr, nvec_repeats_gr$class)

# Count the number of repeats for each class
repeat_count = table(nvec_repeats_gr$class)

# Get proportion of Nematostella genome covered by each class of genomic region
repeat_class_proportions = sort(sapply(repeat_annotation_grl, function(x) sum(width(reduce(x, ignore.strand = T)))))/nvec_geneome_length
repeat_class_proportions = data.frame(
  region = names(repeat_class_proportions),
  proportion = repeat_class_proportions, row.names = NULL
)
repeat_class_proportions = filter(repeat_class_proportions, region != "ARTEFACT")
repeat_class_proportions = arrange(repeat_class_proportions, proportion)
repeat_class_proportions$region = 
  factor(repeat_class_proportions$region, levels = repeat_class_proportions$region)
x_labels = as.character(repeat_class_proportions$region)

# Append number of repeats to each class
x_labels = paste(x_labels, paste("n =", prettyNum(repeat_count[x_labels], big.mark = ",")), sep = "\n")

# Create a barplot of the genomic proportion of different classes of genomic regions
repeat_class_proportions_plot = ggplot(repeat_class_proportions, aes(x = region, y = proportion, fill = region)) +
  geom_col(color="black", fill = colorRampPalette(RColorBrewer::brewer.pal(n = 20, name = "RdBu"))(10)) +
  scale_fill_brewer(palette="Set1")
repeat_class_proportions_plot = customize_ggplot_theme(repeat_class_proportions_plot, 
  xlab = "Repeat Family", ylab = "Proportion", title = "Proportion of Nematostella Genome Covered\nby Different Repeat Families",
  x_labels = x_labels) +
  coord_flip()
repeat_class_proportions_plot
ggsave(plot = repeat_class_proportions_plot, "repeat_class_proportions_barplot.pdf", width = 9, height = 9)

### Make a piechart

# Make a GRanges for the Nematostella genome
nematostella_genome_gr = GRanges(seqinfo(BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19::BSgenome.Nvectensis.NCBI.jaNemVect1.1.lambda.pUC19))[1:15]

# Define introns as only regions not overlapping exons
introns_gr = setdiff(complete_annotation_grl$Intron, complete_annotation_grl$Exon, ignore.strand = T)

# Find intronic sequences with and without repeats
introns_with_repeats = intersect(introns_gr, complete_annotation_grl$`Repetitive DNA`, ignore.strand = T)
introns_without_repeats = setdiff(introns_gr, complete_annotation_grl$`Repetitive DNA`, ignore.strand = T)

# Find intergenic regions with and without repeats
intergenic_regions = GenomicRanges::setdiff(nematostella_genome_gr, c(promoters_exons_and_introns_gr, non_coding_regions_gr), ignore.strand = T)
intergenic_regions_with_repeats = intersect(intergenic_regions, complete_annotation_grl$`Repetitive DNA`, ignore.strand = T)
intergenic_regions_without_repeats = setdiff(intergenic_regions, complete_annotation_grl$`Repetitive DNA`, ignore.strand = T)

# Create a list with the regions for the piechart
piechart_list = list(
  Exon = complete_annotation_grl$Exon, 
  "Introns (Non Repetitive)" = introns_without_repeats,
  "Intronic Repeats" = introns_with_repeats,
  "Intergenic Repeats" = intergenic_regions_with_repeats,
  "Intergenic (Non Repetitive)" = intergenic_regions_without_repeats
  )

# 
piechart_regions_class_proportions = sapply(piechart_list, function(x) sum(width(reduce(x, ignore.strand = T))))/nvec_geneome_length
piechart_regions_class_proportions = data.frame(
  region = names(piechart_regions_class_proportions),
  proportion = piechart_regions_class_proportions, row.names = NULL
)
piechart_regions_class_proportions$proportion = round(piechart_regions_class_proportions$proportion/sum(piechart_regions_class_proportions$proportion), 2)
piechart_regions_class_proportions$label = paste0(piechart_regions_class_proportions$region, " - ", piechart_regions_class_proportions$proportion, "%")
piechart_regions_class_proportions$region = factor(piechart_regions_class_proportions$region, levels = piechart_regions_class_proportions$region)

pie = ggplot(piechart_regions_class_proportions, aes(x = "", y = proportion, fill = region)) +
  geom_col(width=1, color = "white")
pie = pie + coord_polar("y", start = 0) + theme_void() + labs(fill = NULL) 
pie = pie + scale_fill_manual(labels = piechart_regions_class_proportions$label, values = rev(colorRampPalette(RColorBrewer::brewer.pal(n = 20, name = "RdBu"))(5))[c(1, 2, 4, 5, 3)])
pie = pie + theme(legend.text = element_text(size = 18), legend.key.size = unit(1.5, "cm"))

ggsave(plot = pie, "nematostella_genome_piechart.pdf", width = 16, height = 9)

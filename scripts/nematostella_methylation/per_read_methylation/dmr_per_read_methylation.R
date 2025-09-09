# Make plots of per-read methylation in Nv_Cr5_3_4mpf, NvIT16K4_sperm_Oct31 and Nv_Cr18 (an untreated cross served as control)

# Load required packages
library(dplyr)
library(GenomicRanges)
#source("../auxillary_scripts/plotting_functions.R")
source("../../auxillary_scripts/modkit_functions.R")
source("../../auxillary_scripts/alignments_functions.R")
source("../../auxillary_scripts/methylartist_functions.R")

# Load DMSO DMRs 
dmso_dmrs = readRDS("../dmso_dmrs.rds")

# Get paths to modBAMs for Nv_Cr5_3_4mpf, NvIT16K4_sperm_Oct31, Nv_Cr5 and Nv_Cr18
Nv_Cr5_modbam = "/data/SBCS-ademendoza/07-Lan/Nanopore/Nv_Cr5/Nv_Cr5_sorted.bam"
Nv_Cr5_3_4mpf_modbam = "/data/SBCS-ademendoza/07-Lan/Nanopore/Nv_Cr5_3_4mpf/Nv_Cr5_3_4mpf_sorted.bam"
NvIT16K4_sperm_Oct31_modbam = "/data/SBCS-ademendoza/07-Lan/Nanopore/NvIT16K4_sperm_Oct31/NvIT16K4_sperm_Oct31_sorted.bam"
Nv_Cr18_modbam = "/data/SBCS-ademendoza/07-Lan/Nanopore/Nv_Cr18/Nv_Cr18_sorted.bam"
Nv_Cr5_3_4mpf_haplotagged_modbam = "/data/SBCS-ademendoza/07-Lan/Nanopore/Nv_Cr5_3_4mpf/Nv_Cr5_3_4mpf_sorted_haplotagged.bam"
Nv_Cr5_3_4mpf_pedigree_haplotagged_modbam = "/data/SBCS-ademendoza/07-Lan/Nanopore/Nv_Cr5_3_4mpf/Nv_Cr5_3_4mpf_sorted_pedigree_haplotagged.bam"

# Get path to Nematostella genome FASTA
nvec_fasta = "../../nematostella_genome/Nematostella_DToL_lambda_pUC_mitochondria_originalIDs.fasta"

# Extract mod calls for DMRs from Nv_Cr5. Then find overlaps between mod calls and DMRs and indicate if mod calls are for fast or slow DMRs and save table. Took 7 minutes
system.time({Nv_Cr5_all_dmr_mod_calls = extract_mod_calls(modBAM = Nv_Cr5_modbam, genomic_regions = dmso_dmrs,
  nthreads = 12, reference_fasta = nvec_fasta, add_genomic_region_index = T)})
Nv_Cr5_all_dmr_mod_calls$dmr_status = dmso_dmrs$dmr_status[Nv_Cr5_all_dmr_mod_calls$genomic_region_index]
data.table::fwrite(Nv_Cr5_all_dmr_mod_calls, "Nv_Cr5_all_dmrs_mod_calls.tsv.gz", sep = "\t")

# Extract mod calls for DMRs from Nv_Cr5_3_4mpf. Then find overlaps between mod calls and DMRs and indicate if mod calls are for fast or slow DMRs and save table. Took 7 minutes
system.time({Nv_Cr5_3_4mpf_all_dmr_mod_calls = extract_mod_calls(modBAM = Nv_Cr5_3_4mpf_modbam, genomic_regions = dmso_dmrs,
  nthreads = 12, reference_fasta = nvec_fasta, add_genomic_region_index = T)})
Nv_Cr5_3_4mpf_all_dmr_mod_calls$dmr_status = dmso_dmrs$dmr_status[Nv_Cr5_3_4mpf_all_dmr_mod_calls$genomic_region_index]
data.table::fwrite(Nv_Cr5_3_4mpf_all_dmr_mod_calls, "Nv_Cr5_3_4mpf_all_dmrs_mod_calls.tsv.gz", sep = "\t")

# Extract mod calls for DMRs from NvIT16K4_sperm_Oct31. Then find overlaps between mod calls and DMRs and indicate if mod calls are for fast or slow DMRs and save table. Took 7 minutes
system.time({NvIT16K4_sperm_Oct31_all_dmr_mod_calls = extract_mod_calls(modBAM = NvIT16K4_sperm_Oct31_modbam, genomic_regions = dmso_dmrs, 
  nthreads = 12, reference_fasta = nvec_fasta, add_genomic_region_index = T)})
NvIT16K4_sperm_Oct31_all_dmr_mod_calls$dmr_status = dmso_dmrs$dmr_status[NvIT16K4_sperm_Oct31_all_dmr_mod_calls$genomic_region_index]
data.table::fwrite(NvIT16K4_sperm_Oct31_all_dmr_mod_calls, "NvIT16K4_sperm_Oct31_all_dmrs_mod_calls.tsv.gz", sep = "\t")

# Extract mod calls for DMRs from Nv_Cr18. Then find overlaps between mod calls and DMRs and indicate if mod calls are for fast or slow DMRs and save table. Took 5 minutes
system.time({Nv_Cr18_all_dmr_mod_calls = extract_mod_calls(modBAM = Nv_Cr18_modbam, genomic_regions = dmso_dmrs, 
  nthreads = 12, reference_fasta = nvec_fasta, add_genomic_region_index = T)})
Nv_Cr18_all_dmr_mod_calls$dmr_status = dmso_dmrs$dmr_status[Nv_Cr18_all_dmr_mod_calls$genomic_region_index]
data.table::fwrite(Nv_Cr18_all_dmr_mod_calls, "Nv_Cr18_all_dmrs_mod_calls.tsv.gz", sep = "\t")

# Read table for Nv_Cr5 and filter for quality calls for all DMRs and summarize mod calls per read. Took 6 minutes
system.time({Nv_Cr5_all_dmr_mod_calls = data.table::fread("Nv_Cr5_all_dmrs_mod_calls.tsv.gz")})
system.time({Nv_Cr5_all_dmr_mod_calls = filter_modkit_calls(Nv_Cr5_all_dmr_mod_calls)})
Nv_Cr5_all_dmr_mod_calls = dplyr::group_by(Nv_Cr5_all_dmr_mod_calls, dmr_status)
system.time({Nv_Cr5_all_dmr_mod_calls_per_read = summarize_mod_calls_per_read(Nv_Cr5_all_dmr_mod_calls)})
data.table::fwrite(Nv_Cr5_all_dmr_mod_calls_per_read, "Nv_Cr5_all_dmr_mod_calls_per_read.tsv.gz", sep = "\t")

# Read table for Nv_Cr5_3_4mpf and filter for quality calls for all DMRs and summarize mod calls per read. Took 6 minutes
system.time({Nv_Cr5_3_4mpf_all_dmr_mod_calls = data.table::fread("Nv_Cr5_3_4mpf_all_dmrs_mod_calls.tsv.gz")})
system.time({Nv_Cr5_3_4mpf_all_dmr_mod_calls = filter_modkit_calls(Nv_Cr5_3_4mpf_all_dmr_mod_calls)})
Nv_Cr5_3_4mpf_all_dmr_mod_calls = dplyr::group_by(Nv_Cr5_3_4mpf_all_dmr_mod_calls, dmr_status)
system.time({Nv_Cr5_3_4mpf_all_dmr_mod_calls_per_read = summarize_mod_calls_per_read(Nv_Cr5_3_4mpf_all_dmr_mod_calls)})
data.table::fwrite(Nv_Cr5_3_4mpf_all_dmr_mod_calls_per_read, "Nv_Cr5_3_4mpf_all_dmr_mod_calls_per_read.tsv.gz", sep = "\t")

# Read table for NvIT16K4_sperm_Oct31 and filter for quality calls for all DMRs and summarize mod calls per read. Took 6 minutes
system.time({NvIT16K4_sperm_Oct31_all_dmr_mod_calls = data.table::fread("NvIT16K4_sperm_Oct31_all_dmrs_mod_calls.tsv.gz")})
system.time({NvIT16K4_sperm_Oct31_all_dmr_mod_calls = filter_modkit_calls(NvIT16K4_sperm_Oct31_all_dmr_mod_calls)})
NvIT16K4_sperm_Oct31_all_dmr_mod_calls = dplyr::group_by(NvIT16K4_sperm_Oct31_all_dmr_mod_calls, dmr_status)
system.time({NvIT16K4_sperm_Oct31_all_dmr_mod_calls_per_read = summarize_mod_calls_per_read(NvIT16K4_sperm_Oct31_all_dmr_mod_calls)})
data.table::fwrite(NvIT16K4_sperm_Oct31_all_dmr_mod_calls_per_read, "NvIT16K4_sperm_Oct31_all_dmr_mod_calls_per_read.tsv.gz", sep = "\t")

# Read table for Nv_Cr18 and filter for quality calls for all DMRs and summarize mod calls per read. Took 6 minutes
system.time({Nv_Cr18_all_dmr_mod_calls = data.table::fread("Nv_Cr18_all_dmrs_mod_calls.tsv.gz")})
system.time({Nv_Cr18_all_dmr_mod_calls = filter_modkit_calls(Nv_Cr18_all_dmr_mod_calls)})
Nv_Cr18_all_dmr_mod_calls = dplyr::group_by(Nv_Cr18_all_dmr_mod_calls, dmr_status)
system.time({Nv_Cr18_all_dmr_mod_calls_per_read = summarize_mod_calls_per_read(Nv_Cr18_all_dmr_mod_calls)})
data.table::fwrite(Nv_Cr18_all_dmr_mod_calls_per_read, "Nv_Cr18_all_dmr_mod_calls_per_read.tsv.gz", sep = "\t")

# Load tables with per-read methylation calls for all DMRs
NvIT16K4_sperm_Oct31_all_dmr_mod_calls_per_read = data.table::fread("NvIT16K4_sperm_Oct31_all_dmr_mod_calls_per_read.tsv.gz")
Nv_Cr18_all_dmr_mod_calls_per_read = data.table::fread("Nv_Cr18_all_dmr_mod_calls_per_read.tsv.gz")
Nv_Cr5_3_4mpf_all_dmr_mod_calls_per_read = data.table::fread("Nv_Cr5_3_4mpf_all_dmr_mod_calls_per_read.tsv.gz")
Nv_Cr5_all_dmr_mod_calls_per_read = data.table::fread("Nv_Cr5_all_dmr_mod_calls_per_read.tsv.gz")

# Read in per read mod calls tables and combine into a single table
combined_per_read_mod_summary = dplyr::bind_rows(list(
  Nv_Cr18 = Nv_Cr18_all_dmr_mod_calls_per_read,
  NvIT16K4_sperm_Oct31 = NvIT16K4_sperm_Oct31_all_dmr_mod_calls_per_read,
  #Nv_Cr5_3_4mpf = Nv_Cr5_3_4mpf_all_dmr_mod_calls_per_read,
  Nv_Cr5 = Nv_Cr5_all_dmr_mod_calls_per_read
), .id = "sample")

# Round methylation to 1 decimal place
combined_per_read_mod_summary$meth_bin = round(combined_per_read_mod_summary$m, 1)

# Create a vector of sample names for the plot
sample_labels = c("Untreated x Untreated Cross", "GSK Treated Sperm", "GSK Treated x Untreated Cross")

# Specify order for samples
combined_per_read_mod_summary$sample = factor(combined_per_read_mod_summary$sample, levels = c("Nv_Cr18", "NvIT16K4_sperm_Oct31", "Nv_Cr5"))

# Create a barplot of the per-read methylation of slow DMRs
slow_dmr_per_read_barplot = ggplot(filter(combined_per_read_mod_summary, dmr_status == "Slow"), aes(x = meth_bin, y = after_stat(prop), fill = sample)) +
  geom_bar(color = "black")
slow_dmr_per_read_barplot = customize_ggplot_theme(slow_dmr_per_read_barplot, title = "Methylation of Reads Overlapping Slow DMRs", 
  xlab = "Methylation", ylab = "Proportion of Reads", 
  fill_colors = c("#A28CB1", "#D2C465", "#B1A39E"), show_legend = F,
  facet = "sample", facet_scales = "fixed", facet_nrow = 3, facet_labels = sample_labels) +
  theme(strip.background = element_rect(fill = "white"))
slow_dmr_per_read_barplot
ggsave(plot = slow_dmr_per_read_barplot, "plots/slow_dmrs_per_read_distribution_barplots_Nv_Cr5_gastrula.pdf", width = 16, height = 9)

### Make boxplots of methylation of haplotyped reads in Nv_Cr5_3_4mpf and compare to NvIT16K4_sperm_Oct31 and Nv_Cr18

# Make plot of per-read methylation for slow and fast DMRs from Nv_Cr18
Nv_Cr18_per_read_dmr_meth_boxplots = ggplot(Nv_Cr18_all_dmr_mod_calls_per_read, aes(y = m, x = dmr_status)) + 
  geom_boxplot(fill = "#A28CB1")
Nv_Cr18_per_read_dmr_meth_boxplots = customize_ggplot_theme(Nv_Cr18_per_read_dmr_meth_boxplots, 
  title = "Untreated Male X\nUntreated Female Gastrulas", plot_title_size = 18) 
Nv_Cr18_per_read_dmr_meth_boxplots

# Make plot of per-read methylation for slow and fast DMRs from NvIT16K4_sperm_Oct31
NvIT16K4_sperm_Oct31_per_read_dmr_meth_boxplots = ggplot(NvIT16K4_sperm_Oct31_all_dmr_mod_calls_per_read, aes(y = m, x = dmr_status)) + 
  geom_boxplot(fill = "#D2C465")
NvIT16K4_sperm_Oct31_per_read_dmr_meth_boxplots = customize_ggplot_theme(NvIT16K4_sperm_Oct31_per_read_dmr_meth_boxplots, 
  title = "Untreated Male X\nUntreated Female Gastrulas", plot_title_size = 18) 
NvIT16K4_sperm_Oct31_per_read_dmr_meth_boxplots

# Filter Nv_Cr5_3_4mpf_sorted_pedigree_haplotagged.bam for reads overlapping DMRs. Took 2 minutes.
sbp = ScanBamParam(what = "qname", tag = c("PS", "HP"))
system.time({dmso_dmr_read_haplotags = filter_bam_for_regions(Nv_Cr5_3_4mpf_pedigree_haplotagged_modbam, 
  genomic_regions = dmso_dmrs, param = sbp, nthreads = 20)})
saveRDS(dmso_dmr_read_haplotags, "dmso_dmr_read_haplotags.rds")
dmso_dmr_read_haplotags = readRDS("dmso_dmr_read_haplotags.rds")

# Create a vector matching reads to their haplotype
read_haplotypes = setNames(as.character(dmso_dmr_read_haplotags$tag$HP), dmso_dmr_read_haplotags$qname)

# Indicate haplotype of reads
Nv_Cr5_3_4mpf_all_dmr_mod_calls_per_read$haplotype = read_haplotypes[Nv_Cr5_3_4mpf_all_dmr_mod_calls_per_read$read_id]

# Create boxplots of methylation by haplotype and DMR status
Nv_Cr5_3_4mpf_per_read_dmr_meth_boxplots = ggplot(filter(Nv_Cr5_3_4mpf_all_dmr_mod_calls_per_read, !is.na(haplotype)), aes(fill = haplotype, y = m, x = dmr_status)) +
  geom_boxplot()
Nv_Cr5_3_4mpf_per_read_dmr_meth_boxplots = customize_ggplot_theme(Nv_Cr5_3_4mpf_per_read_dmr_meth_boxplots, 
  title = "GSK Treated Male X\nUntreated Female 4 MPF Individual", plot_title_size = 18,
  fill_title = "Haplotype", fill_labels = c("Paternal", "Maternal"), fill_colors = c("#D2C465", "#A28CB1"), show_legend = T)
Nv_Cr5_3_4mpf_per_read_dmr_meth_boxplots

# Combine boxplots for three samples together 
plotlist = list(Nv_Cr18_per_read_dmr_meth_boxplots, NvIT16K4_sperm_Oct31_per_read_dmr_meth_boxplots, Nv_Cr5_3_4mpf_per_read_dmr_meth_boxplots)
combined_boxplots = ggpubr::ggarrange(plotlist = plotlist, nrow = 1, align = "hv", common.legend = T, legend = "right")
combined_boxplots = ggpubr::annotate_figure(combined_boxplots, top = ggpubr::text_grob("Per-Read DMR Methylation Across Samples",  size = 24),
  left = ggpubr::text_grob("Proportion of Methylated CpGs per Read", color = "black", rot = 90, size = 20), 
  bottom = ggpubr::text_grob("DMR Status", color = "black", size = 20))
combined_boxplots
ggsave(plot = combined_boxplots, "plots/haplotyped_dmr_methylation_boxplots.pdf", width = 16, height = 9)

### Make methylartist plots for DMR dmr_2153 (NC_064035.1:4626543-4641394)

system.time({methylartist_region_plot(modBAM = Nv_Cr5_3_4mpf_haplotagged_modbam, region = dmso_dmrs["dmr_2153"], 
  genome_fasta = nvec_fasta, outfile = "plots/Nv_Cr5_3_4mpf_dmr_2153_methylartist_plot.svg", 
  additional_options = "--hidelegend --phased --svg")})

system.time({methylartist_locus_plot(modBAM = NvIT16K4_sperm_Oct31_modbam, region = dmso_dmrs["dmr_2153"], 
  genome_fasta = nvec_fasta, outfile = "plots/NvIT16K4_sperm_Oct31_dmr_2153_methylartist_plot.svg", additional_options = "--hidelegend --svg")})

system.time({methylartist_locus_plot(modBAM = Nv_Cr18_modbam, region = dmso_dmrs["dmr_2153"], 
  genome_fasta = nvec_fasta, outfile = "plots/Nv_Cr18_dmr_2153_methylartist_plot.svg", additional_options = "--hidelegend --svg")})

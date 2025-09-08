#EM-seq
# packages needed: BSSeeker2, fastp, bbmap, sambamba, bamUtils, samtools. cgmaptools

# create indexes files
bs_seeker2-build.py -f Genome.fasta --aligner=bowtie2 -d /indexes

## Mapping and get CG rates

# Fastq input
readR1="$1"
readR2="$2"

# Set base name for output
baseR1=$(basename "$readR1" .fq.gz) &&
baseR2=$(basename "$readR2" .fq.gz) &&

# Trim reads and QC report
fastp -i read1.fastq.gz -o read1_trimmed.fastq.gz -I read2.fastq.gz -O read2_trimmed.fastq.gz \
--html=read1_fastp_report.html &&

# Merge overlapping read pairs
bbmerge.sh in1=read1_trimmed.fastq.gz in2=read2_trimmed.fastq.gz qtrim=r \
out=read1_merged.fastq.gz outu1=read1_unmerged.fastq.gz outu2=read2_unmerged.fastq.gz \
-Xmx20g &&

echo read1_unmerged.fastq.gz
echo read2_unmerged.fastq.gz

###### Paired-end alignment
bs_seeker2-align.py \
--aligner=bowtie2 --bt2--end-to-end --bt2-p "$cores" -e 300 -X 2000 -m 4 \
-1 read1_unmerged.fastq.gz -2 read2_unmerged.fastq.gz -o read1_unmerged.PE.bam \
-d index_path/ \
-g Genome.fasta &&

###### Merged reads single-end alignment
bs_seeker2-align.py \
--aligner=bowtie2 --bt2--end-to-end --bt2-p "$cores" -e 400 -m 6 \
-i read1_merged.fastq.gz -o read2_merged.SE.bam \
-d index_path/ \
-g Genome.fasta &&


###### Sort the output bam files
sambamba sort -t "$cores" read1_unmerged.PE.bam  &&
sambamba sort -t "$cores" read1_merged.SE.bam  &&

### Clip the overlaps in the sorted BAM files using bamUtil
bam clipOverlap --in read1_unmerged.PE.sorted.bam --out read1_unmerged.PE.sorted.clipped.bam &&

#### Merge all bam files
sambamba merge -t "$cores" read1.bam read1_merged.SE.sorted.bam read1_unmerged.PE.sorted.clipped.bam &&

# Create md5 sum for output bam file
md5sum read1.bam > read1.bam.md5 &&


# Cleanup
rm fastp.json read1_trimmed.fastq.gz \
read2_trimmed.fastq.gz \
read1_merged.fastq.gz \
read1_unmerged.fastq.gz \
read2_unmerged.fastq.gz \
read1_unmerged.PE.sorted.clipped.bam \
read1_unmerged.PE.sorted.bam read1_unmerged.PE.sorted.bam.bai \
read1_merged.SE.sorted.bam read1_merged.SE.sorted.bam.bai


#if there are multiple lane reads, after mapping different lane reads separately, do
samtools merge bothL.merged.bam L1_1.bam L3_1.bam

# coverting bam file to cgmap file 

sambamba markdup -r -t 12 bothL.merged.bam bothL.dedup.bam &&
cgmaptools convert bam2cgmap \
--bam bothL.dedup.bam --genome $genome_file -o bothL

#calculate methylation levels
for i in *.CGmap.gz; 
	do name=$(echo $i | perl -pe "s/.CGmap.gz//"); 
	dul=$(zcat $i | awk '$4 == "CG" && $1 != "chrL" && $1 != "chrP" && $1 != "chrM"' | awk '{cov+=$8; mC+=$7;} END{print "mCG", cov, mC, mC/cov*100;}'); 
	echo $name $dul ; 
done > mCG_global_stats

for i in *.CGmap.gz; do 
	name=$(echo $i | perl -pe "s/.CGmap.gz//"); 
	dul=$(zcat $i | awk '$1 == "chrL"' | awk '{cov+=$8; mC+=$7;} END{print "lambda", cov, mC, mC/cov*100;}'); 
	echo $name $dul ; 
done > lambda_global_stats 

for i in *.CGmap.gz; do 
	name=$(echo $i | perl -pe "s/.CGmap.gz//"); 
	dul=$(zcat $i | awk '$1 == "chrM"' | awk '{cov+=$8; mC+=$7;} END{print "Mitochondria", cov, mC, mC/cov*100;}'); 
	echo $name $dul ; 
done > Mitochondria_global_stats

for i in *.CGmap.gz; do 
	name=$(echo $i | perl -pe "s/.CGmap.gz//"); 
	dul=$(zcat $i | awk '$1 == "chrP" && $4 == "CG"' | awk '{cov+=$8; mC+=$7;} END{print "pUC19_mCG", cov, mC, mC/cov*100;}'); 
	echo $name $dul ; 
done > pUC19_mCG_global_stats


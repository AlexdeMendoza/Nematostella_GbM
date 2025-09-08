#Mapping, calling and counting peaks for ATAC-seq

#packages required: chromap, MACS3, samtools, pigz

#use chromap to generate a index file for genome first
chromap -i -r /data/SBCS-ademendoza/Annotations/Nematostella/GCF_932526225.1_jaNemVect1.1_genomic.fna -o GCF_932526225_chromap_index

#maps ATAC-seq with chromap
# map to generate sam file
chromap --preset atac -x GCF_932526225_chromap_index -r genome.fasta -1 read1 -2 read2 --SAM -o output.sam -t $cores

# map to generate bed file (technically they are the same but this is so quick that it is fine redoing it)
chromap --preset atac -x GCF_932526225_chromap_index -r genome.fasta -1 read1 -2 read2 -o output.align.bed -t $cores

# compress bed 
pigz -p $cores output.align.bed

# make sam into a bam for bigwig generation and cleanup
samtools view -@ $cores -b output.sam -o output.bam
samtools index output.bam
bamCoverage -b output.bam -o output.bigwig -p $cores --normalizeUsing CPM
rm output.sam

# calling peaks for all samples
for i in *.align.bed.gz; do
 macs3 callpeak -t $i -n ${i%%.align.bed.gz}_model05_macs3 -f BEDPE -g $genome_size -q 0.05 -B
done

# merge reads like 
cat *narrowPeak  | cut -f 1,2,3 | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - >  all_peaks.bed
# then conver to "saf" format:
awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid="macs3Peak_"++nr;  print peakid,$1,$2,$3,"."}' all_peaks.bed  > all_peaks.saf

# count reads on peaks for each sample with featureCounts
featureCounts -p -F SAF -a all_peaks.saf --fracOverlap 0.2 -o output.allpeaks.counts output.bam -T 10 

# make counts matrix
 cut -f 1,2,3,4 output.allpeaks.counts > a
 for i in *counts; do cut -f 7 $i > $i.dul; done
 paste a *dul > all_peaks.counts

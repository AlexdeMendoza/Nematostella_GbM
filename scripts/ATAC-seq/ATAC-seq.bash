#Mapping and calling peaks for ATAC-seq

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
#Process RNA-seq data with hisat2

#package needed: hisat2, sambamba

#build index
hisat2-build genome.fasta genome_index

#map reads
hisat2 --max-intronlen 40000 --dta -N 1 -p $cores \
--summary-file output.hisat2.log -x genome_index \
-1 read1.fastq.gz -2 read2.fastq.gz -S output.hisat2.sam

sambamba view -t $cores -f bam -S -o output.hisat2.sam output.hisat2.sam

sambamba sort -t $cores -o output.hisat2.sorted.bam output.hisat2.bam

rm output.hisat2.sam output.hisat2.bam
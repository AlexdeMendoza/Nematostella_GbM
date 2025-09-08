# bash script to calculate global methylation for low passage methylome profilling by multiplexed Nanopore sequencing

#packages required: guppy; pigz; samtools; minimap2; modbam2bed


# first gather all the fast5 files and put them into a folder called fast5
# 01 basecalling by guppy 6.3.8 using GPU

nanopore_guppy guppy_basecaller -i fast5/ -s guppy_GPU_R10_PoreStella/ \
-d /data/SBCS-ademendoza/01-ademendoza/software/newGuppy/ont-guppy/data/ -c dna_r10.4.1_e8.2_400bps_modbases_5mc_cg_sup_prom.cfg \
--device auto --bam_out --recursive \
--trim_adapters --barcode_kits SQK-RBK114-24 --do_read_splitting

# 02 merge bam files for each barcode
# after running Guppy, get into the folder it has created, and then the folder named "pass"

for i in barcode{01..10}; do 
        pigz ${i}/*fastq;
        samtools cat -o ${i}.merge_unal.bam ${i}/*.bam
done

# 03 Map reads using minimap2
for i in barcode{01..10}; do
	samtools fastq -T MM,ML ${i}.merge_unal.bam | \
	minimap2 -y -ax map-ont /data/SBCS-ademendoza/Annotations/Nematostella/Nematostella_DToL_lambda_pUC_mitochondria_originalIDs.fasta - -t 32 | \
	samtools view -u - | samtools sort -@ 32 -o ${i}.merge.modmapped.bam; 
	samtools index ${i}.merge.modmapped.bam -@ 32 ;
done

# 04 call methylation

for i in barcode{01..10}; do
modbam2bed -m 5mC -t 32 \
-e -p ${i}.merge.modmapped \
--cpg \
/data/SBCS-ademendoza/Annotations/Nematostella/Nematostella_DToL_lambda_pUC_mitochondria_originalIDs.fasta \
${i}.merge.modmapped.bam > ${i}.merge.modmapped.mCG.bed ;
done;

# 05 Calculate methylation levels

for i in barcode{01..10}; do
	mCG=$(awk '$1 != "chrM" && $1 != "chrP" && $1 != "chrL" && $5>0 {can+=$12; mod+=$13} END{print mod,can+mod,100*(mod/(can+mod))}' ${i}.merge.modmapped.mCG.bed);
	mito=$(awk '$1 == "chrM" && $5>0 {can+=$12; mod+=$13} END{print mod,can+mod,100*(mod/(can+mod))}' ${i}.merge.modmapped.mCG.bed);
	lambda=$(awk '$1 == "chrL" && $5>0 {can+=$12; mod+=$13} END{print mod,can+mod,100*(mod/(can+mod))}' ${i}.merge.modmapped.mCG.bed);
	pUC19=$(awk '$1 == "chrP" && $5>0 {can+=$12; mod+=$13} END{print mod,can+mod,100*(mod/(can+mod))}' ${i}.merge.modmapped.mCG.bed);
	echo $i $mCG $mito $lambda $pUC19 > ${i}.mCG.stats;
done 

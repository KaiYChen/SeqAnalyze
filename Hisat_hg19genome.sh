#!/bin/bash
sample=(5928_5504_26679_C7UWVANXX_HP2con_ATCACG_R1 5928_5504_26680_C7UWVANXX_HP215con_CGATGT_R1 5928_5504_26681_C7UWVANXX_HP2sh3_TTAGGC_R1 5928_5504_26682_C7UWVANXX_HP215sh3_TGACCA_R1 5928_5504_26683_C7UWVANXX_HP215Flagsp3old_ACAGTG_R1 5928_5504_26684_C7UWVANXX_HP215Flagsp3new_GCCAAT_R1)
#sample=(5928_5504_26680_C7UWVANXX_HP215con_CGATGT_R1 5928_5504_26681_C7UWVANXX_HP2sh3_TTAGGC_R1 5928_5504_26682_C7UWVANXX_HP215sh3_TGACCA_R1 5928_5504_26683_C7UWVANXX_HP215Flagsp3old_ACAGTG_R1 5928_5504_26684_C7UWVANXX_HP215Flagsp3new_GCCAAT_R1)

gtf=/home/kc287/reference_geome/hg19/genes.gtf
genome_idx=/home/kc287/reference_geome/hg19/genome
for i in ${sample[@]}
do
	## hisat2
	hisat2 -p 8 -x $genome_idx -U ${i}.fastq.gz -S ./${i}.sam &> ./${i}.sam.info
	## sam to bam
	samtools view -Su ./${i}.sam | samtools sort - ${i}.sorted
	## 3. Making samtool index
	samtools index ${i}.sorted.bam
	## 4. Run HTSeq
	samtools view ${i}.sorted.bam | htseq-count - $gtf > ${i}.count.txt
done

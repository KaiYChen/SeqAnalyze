#!/bin/bash
sample=(10346359_AC_HP215Sh3 10346359_CT_HP2Cont 10346359_GA_HP215Sp3 10346359_TG_HP215Cont)

gtf=/home/kc287/reference_geome/hg19/genes.gtf
genome_idx=/home/kc287/reference_geome/hg19/genome
for i in ${sample[@]}
do
	## hisat2
	hisat2 -p 8 -x $genome_idx -U ${i}.fastq -S ./${i}.sam &> ./${i}.sam.info
	## sam to bam
	samtools view -Su ./${i}.sam | samtools sort - ${i}.sorted
	## 3. Making samtool index
	samtools index ${i}.sorted.bam
	## 4. Run HTSeq
	samtools view ${i}.sorted.bam | htseq-count - $gtf > ${i}.count.txt
done

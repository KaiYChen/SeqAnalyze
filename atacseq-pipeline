#!/bin/bash
# 2017/05/09 

files=(*.gz)
#gtf=/home/kc287/reference_geome/hg19/genes.gtf
genome_idx=~/Research/reference_geome/Homo_sapiens_UCSC_hg19/Bowtie2Index/genome
echo "genome $genome_idx"

## check if genome_idx file path is correct
#if [-f "${genome_idx}.1.bt2"];
#then
#	echo "File $genome_idx exist."
#else
#	echo "File $genome_idx does not exist" >&2
#fi

if [ ! -d ./bam ]; then
	echo "making ./bam "
	mkdir ./bam
fi
if [ ! -d ./rmdupMtUn ]; then
	echo "making ./rmdupMtUn "
	mkdir ./rmdupMtUn
fi
if [ ! -d ./macs2 ]; then
	echo "making ./macs2 "
	mkdir ./macs2
fi
if [ ! -d  ./trim ]; then
	echo "making  ./trim "
	mkdir  ./trim
fi

## pipeline
for filename in "${files[@]}"
do
	echo "processing $filename"
	basename="${filename%.fastq.gz}"
	## trim adaptors
	echo "trimming .............................. $basename"
	trim_galore $filename --output_dir ./trim
	## bowtie2 alignment
	echo "bowtie2 alignment ..................... $basename"
	bowtie2 -p 8 --very-sensitive --maxins 2000 --no-discordant -x $genome_idx -U ./trim/${basename}_trimmed.fq.gz -S ./${basename}.sam &> ./${basename}.sam.info
	## sam to sorted bam
	continue
	echo "sam to sorted bam ..................... $basename"
	samtools view -Su ./${basename}.sam | samtools sort --threads 8 -o ./bam/${basename}.sorted.bam        
	samtools index ./bam/${basename}.sorted.bam
	## remove duplication
	echo "removing duplicates/MT/Un/Random ...... $basename"
	samtools rmdup -s ./bam/${basename}.sorted.bam ./rmdupMtUn/${basename}.sorted.rmdup.bam
	samtools index ./rmdupMtUn/${basename}.sorted.rmdup.bam
	## remove MT/Un/random chromosome
	samtools idxstats ./rmdupMtUn/${basename}.sorted.rmdup.bam | cut -f 1|grep -v _random|grep -v chrM| grep chr | xargs samtools view -b ./rmdupMtUn/${basename}.sorted.rmdup.bam > ./rmdupMtUn/${basename}.sorted.rmdupMtUn.bam
	samtools index ./rmdupMtUn/${basename}.sorted.rmdupMtUn.bam
	rm ./${basename}.sam	
	## make bedGraph
	genomeCoverageBed -bg -ibam ${basename}.sorted.bam -split -g hg19.chrom.sizes >./bedgraph/${basename}sorted.bedGraph
	## make bigwig
	bedGraphToBigWig ./bedgraph/${basename}sorted.bedGraph hg19.chrom.sizes ./bedgraph/${basename}sorted.bw
	echo "MACS2 calling peaks.................... $basename"
	## call peaks
	macs2 callpeak --nomodel --broad -t ./rmdupMtUn/${basename}.sorted.rmdupMtUn.bam -f BAM -g hs --outdir ./macs2 -n ${basename}
	## annoate peaks
	echo "Homer annotating   .................... $basename"
	annotatePeaks.pl ./macs2/${basename}_peaks.broadPeak hg19 -annStats ./macs2/homer/${basename}.annStats.txt > ./macs2/homer/${basename}_peaks_annotate.txt 
done

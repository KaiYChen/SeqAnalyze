# Intall STAR
# https://github.com/alexdobin/STAR
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz
tar -xzf 2.5.2b.tar.gz
cd STAR-2.5.2b

# Alternatively, get STAR source using git
git clone https://github.com/alexdobin/STAR.git
cd STAR

# To include STAR-Fusion
git submodule update --init --recursive

# Build STAR
cd source
make STARforMacStatic

## STAR alignment
#Alternate Protocol 1: Generating Genome Indices
~/local_lib/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ~/Research/reference_geome/Homo_sapiens_UCSC_hg19/WholeGenomeFasta/genome.fa
# https://software.broadinstitute.org/gatk/documentation/article.php?id=3891
# Alignment
genomeDir=~/temp_data/reference_genome/star_hg19_genome/
runDir=~/Desktop/STAR/1pass
gtfFile=~/temp_data/reference_genome/Homo_sapiens_UCSC_hg19/genes.gtf
fastqfile=/Volumes/ShenLabData/Data/HighThroughput/Project_Epigenomic/atacSeq_HumanOrganoid/Chen_3950_170311B1/fastq/KC1_S24_L005_R1_001.fastq.gz
mkdir $runDir
cd $runDir
/opt/STAR-2.5.2b/bin/MacOSX_x86_64/STAR --genomeDir $genomeDir --readFilesIn $fastqfile \
--sjdbGTFfile $gtfFile \
--readFilesCommand gzcat \
--runThreadN 8
## new index is then created using splice junction information contained in the file SJ.out.tab from the first pass
genomeDir=~/Desktop/STAR/2pass
mkdir $genomeDir
genomeFasta=~/temp_data/reference_genome/Homo_sapiens_UCSC_hg19/WholeGenomeFasta/genome.fa
/opt/STAR-2.5.2b/bin/MacOSX_x86_64/STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $genomeFasta \
    --sjdbFileChrStartEnd ~/Desktop/STAR/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN 8
##  index is then used to produce the final alignments    
runDir=~/Desktop/STAR/2pass_final
mkdir $runDir
cd $runDir
/opt/STAR-2.5.2b/bin/MacOSX_x86_64/STAR --genomeDir $genomeDir --readFilesIn $fastqfile --runThreadN 8 --readFilesCommand gzcat

## Add read groups, sort, mark duplicates, and create index
picard AddOrReplaceReadGroups I=Aligned.out.sam O=Aligned_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample 
picard MarkDuplicates I=Aligned_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics 

## Split'N'Trim and reassign mapping qualities
referenceGenome=~/temp_data/reference_genome/Homo_sapiens_UCSC_hg19/WholeGenomeFasta/genome.fa
java -jar /opt/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T SplitNCigarReads -R $referenceGenome -I dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
## Variant Calling
java -jar /opt/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T HaplotypeCaller -R $referenceGenome  -I split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o output.split.vcf
## Variant filtering
## at least 3 SNPs that are within a window of 35 bases between them by adding -window 35 -cluster 3

java -jar /opt/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T VariantFiltration -R $referenceGenome -V output.split.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o output.split.filtered.vcf 

##https://gencore.bio.nyu.edu/variant-calling-pipeline/

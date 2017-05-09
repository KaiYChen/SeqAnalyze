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
fastqfile=./KC1_S24_L005_R1_001.fastq.gz
mkdir $runDir
cd $runDir
/opt/STAR-2.5.2b/bin/MacOSX_x86_64/STAR --genomeDir $genomeDir --readFilesIn $fastqfile \
--sjdbGTFfile $gtfFile \
--readFilesCommand gzcat \
--runThreadN 8

setwd("~/Dropbox/Course/BioinformaticStudyGroup/CountFile")
########### download libraries ##################
####### useful links ##########
### Colors in R : http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
### Graphical Parameters in R: http://www.statmethods.net/advgraphs/parameters.html 
### ggplot2 : http://ggplot2.org
### DESeq2  : https://bioconductor.org/packages/release/bioc/html/DESeq2.html

####### Install Packages & load libraries ######
# install.packages("ggplot2")
library(ggplot2)
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library(DESeq2) 
################ DESeq2 ##############
### file list
sample.fileName<-list.files(path="./htSeq/",pattern ="*.count.txt")

### sample name
sample.sampleName<-c("mock_p1","mock_p2","mock_s1","mock_s2","zika_p1","zika_p2","zika_s1","zika_s2")

### sample conditions
sample.Condition<-c(rep("mock",4),rep("zika",4))

### create a data frame including smaple name, file name, condition
sampleTable<-data.frame(sampleName=sample.sampleName, fileName=sample.fileName, condition=sample.Condition)

### obtain HTSeq input 
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory="./htseq/", design=~condition)
# DESeqDataSetFromMatrix(countData = mockRnaSeqData, colData = as.data.frame(conditions), design = ~ conditions)

### pre filtering : a minimal pre-filtering to remove rows that have only 0 or 1 read
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ]

### assign factors for comparison
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("mock","zika"))
# ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref="mock")

### run Differential expression analysis
dds<-DESeq(ddsHTSeq)
res<-results(dds)
# res<-results(dds,lfcThreshold = 0.5,altHypothesis = "greaterAbs",alpha = 0.05)#Extract results from a DESeq analysis
res<-res[order(res$padj),]
head(res)
summary(res)
mcols(res,use.names=TRUE)# More information on results columns
### output files ###############
res<-as.data.frame(res)
res$label<-rep("unss",dim(res)[1])
res$label[res$log2FoldChange > 1 & res$pvalue < 0.05]<-"ssUP"
res$label[res$log2FoldChange< -1 & res$pvalue < 0.05]<-"ssDown"
write.table(as.data.frame(res),file="HTSeq-DESeq2.out.txt",sep="\t",quote = FALSE,col.names=NA)

################ plot volcano plot ##################
ggplot(data=res, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(alpha=0.3, size= 3, shape=21,
             aes(colour=factor(label), 
                 fill = factor(label))) +
  scale_fill_manual(values=c("royalblue3","firebrick3","gray60")) + 
  scale_colour_manual(values=c("gray80", "gray80","gray80")) +
  theme(legend.position = "none",
        aspect.ratio=1,
        text = element_text(size=20),
        axis.text.y=element_text(colour="grey30"),
        axis.line = element_line(colour = "grey40")) +
  xlim(c(-10,10)) + ylim(c(0,75)) +
  xlab("log2(fold change)") + ylab("-log10(p-value)")+ ggtitle("zika vs mock")

### plotting
plotMA(res)
############### plot genes ##################
plotCounts(dds, gene="LGR5", intgroup="condition")
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", returnData=TRUE)
ggplot(d,aes(x=condition,y=count,fill=condition))+geom_boxplot()+geom_point()
############### plot heatmaps #############
## biocLite("pheatmap")
library(pheatmap)
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:100]# select genes
nt <- normTransform(dds)# defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]# create an expression matrix
pheatmap(log2.norm.counts)
# pheatmap(log2.norm.counts, cluster_rows=FALSE,show_rownames=FALSE,cluster_cols=FALSE)

##############  Principal component plot of the samples ##############
rld <- rlog(dds, blind=FALSE)# regularized log transformation
plotPCA(rld, intgroup=c("condition"))

############## sample-sample distance ##################
library(RColorBrewer)
sampleDists <- dist(t(assay(rld)))        # calculate Euclidean Distance between samples
sampleDistMatrix <- as.matrix(sampleDists)# make a distance matrix
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

############# Reporting tools
#biocLite("ReportingTools")
library(ReportingTools)
des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2',
           title = 'RNA-seq analysis of differential expression using DESeq2',
           reportDirectory = "./reports")
publish(dds,des2Report, pvalueCutoff=0.05,annotation.db="TxDb.Hsapiens.UCSC.hg19.knownGene",
        factor = colData(dds)$condition, 
        reportDir="./reports")
finish(des2Report)

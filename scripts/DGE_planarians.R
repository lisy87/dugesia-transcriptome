# Differential Expression Analysis using DESeq2
# http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start
# http://master.bioconductor.org/packages/release/bioc/html/DESeq2.html

# Installing DESeq2
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

# ERROR: Installation path not writeable, unable to update packages: mgcv, spatial
# https://stackoverflow.com/questions/41839214/installation-path-not-writable-r-unable-to-update-packages
#.libPaths()
#BiocManager::install(c("mgcv", "spatial"), lib = "/home/lisandra/R/x86_64-pc-linux-gnu-library/4.0")
#remove.packages(c("mgcv", "spatial"), lib = "/usr/lib/R/library")

library(DESeq2)
library (magrittr)

#Defining analysis name
comparation <- "name.of.the.comparation_" 

#Importing kallisto data
# https://support.bioconductor.org/p/120392/
# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

#BiocManager::install("tximportData", lib = "your PATH to library")
library(tximportData)
#BiocManager::install("tximport", lib = "your PATH to library")
library(tximport)

setwd("/media/lisandra/Data/DATA_LISY/PhD_Sexualidad/ANALISIS/TRANSCRIPTOME/NEW_ANALISIS_nov20/SEX_ASEX/NEW_ANALS_june21/ABUNDANCES")
files <- c("sample_1_kall_abundance.tsv", #importing abundance files
           "sample_2_kall_abundance.tsv",
           "sample_n_kall_abundance.tsv")

names(files) <- c("sample_1", "sample_2", "sample_n")

# Input kallisto object

txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)

setwd("Working_Directory for the comparation")

# Generate the sample condition file
# This file will determine the comparison levels.
sample_info <- read.delim(paste0(comparation,"condition_sample.txt"))
sample_info

# generate the DESeqDataSet
#DESeq.ds <- DESeqDataSetFromMatrix(countData = txi.kallisto , colData = sample_info , design = ~ condition)
DESeq.ds<- DESeqDataSetFromTximport(txi.kallisto, sample_info, ~condition)
#dds <- DESeqDataSetFromTximport(txi.kallisto, sample_info, ~condition)

# checking
colData(DESeq.ds) %>% head
assay(DESeq.ds, "counts") %>% head
rowData(DESeq.ds) %>% head

# test counts
counts(DESeq.ds) %>% str

# retain genes with more than 10 counts
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 10, ]
# test counts again
counts(DESeq.ds) %>% str

# investigate different library sizes
colSums(counts(DESeq.ds))

# calculate the size factor and add it to the data set
DESeq.ds <- estimateSizeFactors(DESeq.ds)
sizeFactors(DESeq.ds)

# check colData and corroborate that the sizeFactors have been add to it
colData(DESeq.ds)

#retrieve the _normalized_ read counts
counts.sf_normalized <- counts(DESeq.ds, normalized = TRUE )

# Log2 transformation of read counts
# transform size - factor normalized read counts to log2 scale using a pseudocount of 1
log.norm.counts <- log2(counts.sf_normalized + 1)

# Graph untransformed and log2-transformed read counts
# About boxplots: https://r-coder.com/boxplot-en-r

#png("boxplot_raw_VS_log2_readcounts")
pdf(paste0(comparation,"boxplot_raw_VS_log2_readcounts.pdf"), width = 8, height = 11)
par(mfrow=c(2,1))
par(mar=c(4,3,3,1), las = 2)
boxplot(counts.sf_normalized, notch = TRUE, main = "untransformed read counts", ylab = "read counts", cex.axis = 0.6)
boxplot(log.norm.counts, notch = TRUE, main = " log2 transformed read counts ", ylab = "log2(read counts)", cex.axis = 0.6)
dev.off()

# graph information:
# https://r-coder.com/boxplot-en-r/#La_funcion_boxplot_en_R
# http://rfunction.com/archives/1302
  
# To graph only log2-transformed read counts
# png("boxplot_log2transf_readcounts")
pdf(paste0(comparation,"boxplot_log2transf_readcounts.pdf"), width = 8, height = 11)
par(mfrow=c(1,1))
par(mar=c(6,3,3,1), las = 2)
boxplot(log.norm.counts, notch = TRUE, main = " log2 transformed read counts ", ylab = "log2(read counts)", cex.axis = 0.6)
dev.off()

# Visually exploring normalized read counts
# To get an impression of how similar read counts are between replicates, it is often insightful to simply plot the counts in a pairwise manner
par(mfrow=c(1,1))
plot(log.norm.counts[ ,1:2], cex =0.1, main = "Normalized log2 ( read counts )")

#Checking for heteroskedasticity
#BiocManager::install(c("vsn"), lib = "your PATH to libraries")
#install.packages("hexbin")
library(vsn)
library(ggplot2)
library(hexbin)

# mean-sd plot
msd_plot1 <- meanSdPlot(log.norm.counts, ranks = FALSE, plot = FALSE)
msd_plot1$gg + ggtitle ("sequencing depth normalized log2 (read counts)") + ylab ("standard deviation")

# obtain normalized log - transformed values
DESeq.rlog <- rlog(DESeq.ds, blind = TRUE)
rlog.norm.counts <- assay(DESeq.rlog)

# mean-sd plot for rlog-transformed data
library(vsn)
library(ggplot2)
msd_plot2 <- meanSdPlot(rlog.norm.counts, ranks = FALSE, plot = FALSE)
msd_plot2$gg + ggtitle("rlog-transformed read counts") + ylab("standard deviation")

# Combining graph

pdf(paste0(comparation,"pairwise_comparisons_replicate_samples.pdf"), width = 8, height = 11)
par(mfrow=c(2,1))
par(mar=c(3,2,3,1))
plot(log.norm.counts[ ,1:2], cex =0.1, main = " Normalized log2 ( read counts )")
plot(rlog.norm.counts[ ,1:2], cex =0.1, main = " Normalized r-log2 ( read counts )")
dev.off()

pdf(paste0(comparation,"mean-sd-plots.pdf"), width = 7, height = 7)
par(mfrow=c(2,1))
msd_plot1$gg + ggtitle ("sequencing depth normalized log2 (read counts)") + ylab ("standard deviation")
msd_plot2$gg + ggtitle("rlog-transformed read counts") + ylab("standard deviation")
dev.off()

# Exploring global read count patterns

#calculating the correlation between columns of a matrix
distance.m_rlog <-as.dist(1 - cor(rlog.norm.counts, method = "pearson"))
# Plotting Dendrogram
pdf(paste0(comparation,"Dendrogram.pdf"), width = 7, height = 7)
plot(hclust(distance.m_rlog), labels = colnames(rlog.norm.counts), main = "rlog-transformed read counts\ndistance : Pearson correlation ", cex.axis =0.6, cex.lab= 0.8, cex=0.7, cex.main=0.8)
dev.off()

# Principal Components Analysis (PCA)
library(DESeq2)
library (ggplot2)
P <- plotPCA(DESeq.rlog)
P <- P + theme_bw() + ggtitle("Rlog transformed counts")
P
# Add names to the points
nudge <- position_nudge(y = 1)
Q <- P + theme_bw() + ggtitle("Rlog transformed counts") + geom_text(aes(label = name), position = nudge, size = 2, color ="black")
Q

#Saving plot
pdf(paste0(comparation,"PCA.pdf"), width = 7, height = 7)
Q
dev.off()


###  DIFFERENTIAL GENE EXPRESSION ANALYSIS!!
# sexual condition as  values used as denominator for the fold change calculation
# DESeq2 uses the levels of the condition to determine the order of the comparison
str(colData(DESeq.ds)$condition)
# set sexual as the first-level-factor
colData(DESeq.ds)$condition <- relevel(colData(DESeq.ds)$condition, "sexual")

# running the DGE analysis
DESeq.ds <- DESeq(DESeq.ds)

# The results() function lets you extract the base means across samples, moderated log2 fold change, standard errors, test statistics etc. for every gene.
DGE.results <- results(DESeq.ds, independentFiltering = TRUE , alpha = 0.05)
summary(DGE.results)

#Saving results
sink(paste0(comparation,'summary_DGE.results.txt'))
cat ("*************************************\n")
cat ("summary DGE results\n")
cat ("*************************************\n")
summary(DGE.results)
cat ("*************************************\n")
cat ("DGE p < 0.05\n")
cat ("*************************************\n")
table(DGE.results$padj < 0.05)
rownames(subset(DGE.results, padj < 0.05))
sink()

upregulated <- rownames(subset(DGE.results, log2FoldChange > 0 & padj < 0.05))
downregulated <- rownames(subset(DGE.results, log2FoldChange < 0 & padj < 0.05))

write(upregulated, file = paste0(comparation,"upregulated_DGE.txt"))
write(downregulated, file = paste0(comparation,"downregulated_DGE.txt"))

# Exploratory plots following DGE analysis

# Histograms
pdf(paste0(comparation,"pvalues_freq.pdf"), width = 7, height = 7)
par(mar=c(5,4,3,1))
hist(DGE.results$pvalue, col = "grey", border = "white", xlab = " ", ylab = " ", main = "frequencies of p-values")
dev.off()

# MA plot 
pdf(paste0(comparation,"MA.pdf"), width = 7, height = 7)
par(mar=c(5,4,3,1), las = 1)
plotMA(DGE.results, alpha = 0.05 , main = " Parental vs Resistant", ylim = c(-4 ,4))
dev.off()

#Volcano plot
pdf(paste0(comparation,"volcanoplot.pdf"), width = 7, height = 7, onefile = FALSE)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(DGE.results)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
dev.off()

# Heatmaps
#install.packages("NMF")
library(NMF)

# sort the results according to the adjusted p-value
DGE.results.sorted <- DGE.results[order(DGE.results$padj) , ]

# identify genes with the desired adjusted p-value cut-off
DGEgenes <- rownames(subset(DGE.results.sorted, padj < 0.05))

#Saving DGEgenes with all data to file
DGEgenes.all <-subset(DGE.results.sorted, padj < 0.05)
write.table(DGEgenes.all, file = paste0(comparation,"DGEgenes_allinfo.csv"))

# extract the normalized read counts for DE genes into a matrix
hm.mat_DGEgenes <- log.norm.counts[DGEgenes, ]

#Saving results in a file
sink(paste0(comparation,'DGEgenes_normreadcounts.txt'))
cat ("normalized read counts for DE genes\n")
cat ("*************************************\n")
hm.mat_DGEgenes
sink()

# Plot the normalized read counts of DE genes sorted by the adjusted p-value
pdf(paste0(comparation,"heatmap_pvalue.pdf"), width = 7, height = 7, onefile = FALSE)
aheatmap(hm.mat_DGEgenes, Rowv = NA , Colv = NA)
dev.off()

# combine the heatmap with hierarchical clustering
pdf(paste0(comparation,"heatmap_clustering.pdf"), width = 7, height = 7, onefile = FALSE)
aheatmap(hm.mat_DGEgenes, Rowv = TRUE, Colv = TRUE, distfun = "euclidean", hclustfun = "average")
dev.off()

# scale the read counts per gene to emphasize the sample-type-specific differences
pdf(paste0(comparation,"heatmap_zscore.pdf"), width = 7, height = 7, onefile = FALSE)
aheatmap(hm.mat_DGEgenes, Rowv = TRUE, Colv = TRUE, distfun = "euclidean", hclustfun = "average", scale = "row")
dev.off()
# values are transformed into distances from the center of the row-specific average:
#(actual value-mean of the group) / standard deviation

# Ordered heatmap
#install.packages('pheatmap')
#https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap  
library(pheatmap)
pheatmap(hm.mat_DGEgenes, Rowv = TRUE, Colv = TRUE, distfun = "euclidean", hclustfun = "average", scale = "row", show_rownames = F, angle_col = 45)

pdf(paste0(comparation,"heatmap_zscore_order.pdf"), width = 7, height = 7, onefile = FALSE)
pheatmap(hm.mat_DGEgenes, Rowv = TRUE, Colv = TRUE, distfun = "euclidean", hclustfun = "average", scale = "row", show_rownames = F, angle_col = 45)
dev.off()

### GO TERM ENRICHMENT ANALYSIS OF DGE ###
# https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-GOSeq

# We need the results of the DESeq2 analysis, specifically the vector of deferentially expressed genes (`DGEgenes`).
# In the DGEgenes object, we have the deferentially expressed genes
# So, we need to kept this information

write(DGEgenes, file = paste0(comparation,"DGE_genes.txt"))

#From here, we can use the obtained information with a script from Trinity
# PATH to/trinityrnaseq-v2.11.0/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling --GO_assignments --lengths --background 
# For that, we need:
###############################################################################################
#
#  --factor_labeling <string>       tab delimited file with format:  factor<tab>feature_id
#   or
#  --genes_single_factor <string>   list of genes to test (can be a matrix, only the first column is used for gene IDs)
#
#  --GO_assignments <string>        extracted GO assignments with format: feature_id <tab> GO:000001,GO:00002,...
#
#  --lengths <string>               feature lengths file with format:  feature_id <tab> length
#
#  --background <string>            gene ids file that defines the full population of genes to consider in testing.
#                                   Ideally, these represent the genes that are expressed and relevant to this test
#                                   as opposed to using all genes in the genome.
#
###############################################################################################

# Obtaining the --factor_labeling <string>
#Use the upregulated_DGE.txt and the downregulated_DGE.txt files to build it
# asexual "\t" upregulated_DGE list
# sexual "\t" downregulated_DGE list

# Obtaining the --GO_assignments <string>
# use the go_annotations.txt obtained previously with 
# PATH to/Trinotate-Trinotate-v3.2.2/util/extract_GO_assignments_from_Trinotate_xls.pl

# Obtaining the --lengths <string>
# use the command line to obtain this file
# awk '/^>/ {print; next; } { seqlen = length($0); print seqlen}' Dvila_3_ref_OGanno.fasta > gene.lengths.txt
# cat gene.lengths.txt | grep ">" > names.txt
# cat gene.lengths.txt | grep -v ">" > lengths.txt
# paste names.txt lengths.txt > gene_lengths.txt
# cat gene_lengths.txt | sed 's/>//' > gene_lengths_OK.txt 

# Obtaining the --background <string>
# use the file obtained with write(DGEgenes, file = paste0(comparation,"DGE_genes.txt"))

# After, you can run the script:
# PATH to/trinityrnaseq-v2.11.0/Analysis/DifferentialExpression/run_GOseq.pl 
#                                         --factor_labeling factor_labeling.txt 
#                                         --GO_assignments go_annotations.txt 
#                                         --lengths gene_lengths_OK.txt 
#                                         --background DGE_genes.txt

# You obtain six files:
#   __runGOseq.R
#   asexual.GOseq.enriched
#   asexual.GOseq.depleted
#   sexual.GOseq.enriched
#   sexual.GOseq.depleted
#   Rplots.pdf

# Two outputs will be generated for each set of genes tested for functional enrichment, one containing the enriched categories, and another containing the depleted categories (see .enriched and .depleted files). The .enriched file will contain those categories that are found to have enriched representation among that set of genes. Similarly, the .depleted file will contain those functional categories that are depleted (under-represented) among that very same set of target genes.

# Obtain the input to REVIGO
# cat asexual.GOseq.enriched | cut -f 1-2 > "comparation"_asexual.GOseq.enriched_to_revigo.txt
# cat sexual.GOseq.enriched | cut -f 1-2 > "comparation"_sexual.GOseq.enriched_to_revigo.txt

## Using REVIGO to summarize results:
# Generating REVIGO's input
# export a list of over-represented GO terms and the corresponding p-value
# You can use the obtained file with the GO category and the pvalue in the REVIGO web: http://revigo.irb.hr/
#REVIGO generates multiple plots; the easiest way to obtain them in high quality is to download the R scripts that it offers and to run those yourself.


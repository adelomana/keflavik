rm(list = ls())

#
# -1. install libraries
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

#
# 0. load libraries
#
library(DESeq2)
library(tximport)
library(biomaRt)
library(BiocParallel)
library(crayon) 
library(ggplot2)

#
# 0. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/Users/adrian/research/keflavik/results/kallisto/kallisto.100"
results_dir = '/Users/adrian/research/keflavik/results/kallisto/deseq2'

#
# 1. generate gene to transcript mapping
#
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                        dataset="hsapiens_gene_ensembl",
                        host = 'https://www.ensembl.org',
                        verbose = TRUE)
working_attributes = c('ensembl_transcript_id', 
                       'ensembl_gene_id', 
                       'external_gene_name',
                       'gene_biotype',
                       'description')
t2g = biomaRt::getBM(attributes=working_attributes, 
                     mart=mart,
                     verbose=TRUE)
dim(t2g)
View(t2g)

#
# 2. define metadata
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
paths = file.path(dirnames, 'abundance.h5')
labels = sapply(strsplit(paths, split='/',fixed=TRUE), function(x) (x[9]))
cell_line = c(rep('M104', 6), rep('M501', 9), rep('Sk28', 9))
treatment = c(rep('empty_vector', 3), rep('setdb2_overexpression', 3), rep('miRNA1', 3), rep('miRNA2', 3), rep('miRNA_control', 3), rep('miRNA1', 3), rep('miRNA2', 3), rep('miRNA_control', 3))

metadata = data.frame(labels)
metadata$cell_line = cell_line
metadata$treatment = treatment
metadata$path = paths
View(metadata)

#
# 3. contrasts
#
threshold = 10
effect_size_threshold = log2(1.6)

#
# 3.1. contrast effect of SETDB2 overexpression in M104
#
rule = metadata$cell_line == 'M104'
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~treatment) 
dds$time = relevel(dds$treatment, ref="empty_vector")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast effect SETDB2 OE in M104:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/effect_SETDB2_OE_M104.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/effect_SETDB2_OE_M104.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('treatment')) + ggtitle('effect SETDB2 OE in M104')
ggsave(file.path(results_dir, 'effect_SETDB2_OE_M104.png'))

#
# 3.2. contrasts for M501
#

# 3.2.1. M501 miRNA1
rule = (metadata$cell_line == 'M501') & (metadata$treatment == 'miRNA1' | metadata$treatment == 'miRNA_control')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~treatment) 
dds$time = relevel(dds$treatment, ref="miRNA_control")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast effect miRNA1 in M501:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/effect_miRNA1_M501.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/effect_miRNA1_M501.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('treatment')) + ggtitle('effect miRNA1 in M501')
ggsave(file.path(results_dir, 'effect_miRNA1_M501.png'))

# 3.2.2. M501 miRNA2
rule = (metadata$cell_line == 'M501') & (metadata$treatment == 'miRNA2' | metadata$treatment == 'miRNA_control')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~treatment) 
dds$time = relevel(dds$treatment, ref="miRNA_control")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast effect miRNA2 in M501:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/effect_miRNA2_M501.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/effect_miRNA2_M501.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('treatment')) + ggtitle('effect miRNA2 in M501')
ggsave(file.path(results_dir, 'effect_miRNA2_M501.png'))

#
# 3.3. contrasts for Sk28
#

# 3.3.1. Sk28 miRNA1
rule = (metadata$cell_line == 'Sk28') & (metadata$treatment == 'miRNA1' | metadata$treatment == 'miRNA_control')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~treatment) 
dds$time = relevel(dds$treatment, ref="miRNA_control")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast effect miRNA1 in Sk28:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/effect_miRNA1_Sk28.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/effect_miRNA1_Sk28.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('treatment')) + ggtitle('effect miRNA1 in Sk28')
ggsave(file.path(results_dir, 'effect_miRNA1_Sk28.png'))

# 3.3.2. Sk28 miRNA2
rule = (metadata$cell_line == 'Sk28') & (metadata$treatment == 'miRNA2' | metadata$treatment == 'miRNA_control')
working_metadata = metadata[rule, ]
dim(working_metadata)
View(working_metadata)

txi = tximport(working_metadata$path, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~treatment) 
dds$time = relevel(dds$treatment, ref="miRNA_control")

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~1)

res = results(dds, parallel=TRUE, alpha=0.05) 
filtred_results = res[which(res$padj < 0.05 & abs(res$log2FoldChange) > effect_size_threshold), ]
sorted_filtred_results = filtred_results[order(filtred_results[["padj"]]),]
anti_results = res[which(res$padj > 0.05 | abs(res$log2FoldChange) < effect_size_threshold), ]
cat(blue(paste('contrast effect miRNA2 in Sk28:', dim(filtred_results)[1], sep=' ')), fill=TRUE)
write.table(sorted_filtred_results, 
            file=paste(results_dir, '/effect_miRNA2_Sk28.tsv', sep=''), quote=FALSE, sep='\t')
write.table(anti_results, 
            file=paste(results_dir, 
                       '/effect_miRNA2_Sk28.anti.tsv', sep=''), quote=FALSE, sep='\t')

plotPCA(rlog(dds), intgroup=c('treatment')) + ggtitle('effect miRNA2 in Sk28')
ggsave(file.path(results_dir, 'effect_miRNA2_Sk28.png'))



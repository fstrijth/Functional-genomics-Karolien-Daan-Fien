library("DESeq2")
library(ggplot2)

# Q2: clinical data vs. morphology
# 1) Compute systematically associations between clinical variables and morphological cluster counts. The purpose is to compare the magnitude of the associations of the different variables with morphology.
# 2) Discuss the association with technical variables.
# 3) For non-technical variables, redo the analysis with adjustment for the confounding technical variables, if any is reported in Q2.2. Report and discuss significant associations.


# Inspired by this tutorial:
#https://lashlock.github.io/compbio/R_presentation.html


#load data
clinical_data = read.delim(file='clinical-data.tsv',sep ='\t',header=TRUE,row.names=1)
clin_categories = clinical_data[c('COHORT', 'DTHHRDY')]
clin_categories$GTEXID = row.names(clin_categories)
rownames(clin_categories) = c()
morph_counts = t(read.delim(file='morphological-counts.tsv', sep='\t', header=TRUE, row.names=1))
row.names(morph_counts) = gsub("Mophological.cluster.", "", row.names(morph_counts))


dds = DESeqDataSetFromMatrix(countData = morph_counts, 
                             colData = clin_categories,
                             design = ~COHORT)
dds = DESeq(dds)
res = results(dds)
head(results(dds, tidy=TRUE))
summary(res)
res <- res[order(res$padj),]
head(res)


par(mfrow=c(2,3))
clusters = c("57", "33", "48", "12", "35", "5")
for (cluster in clusters) {
  plotCounts(dds, gene=cluster, intgroup='COHORT')
}

par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

library("edgeR")
library(stringr)
library("xlsx")

# Analysis based on this tutorial:
# https://bioinfo-dirty-jobs.github.io/rana2/lectures/08.rnaseq_edger/

# Read the data files
rna_counts_raw = read.delim(file='RNA-read-counts.tsv',sep='\t',header=TRUE,row.names=1)
# Each row corresponds to a gene, each column corresponds to a sample
#head(rna_counts)
# Clean the data: remove the description
rna_counts = subset(rna_counts_raw, select=-c(Description))

#rna counts uses dot in gtex id, clinical data uses line --> use dot everywhere
clin_data = read.delim(file='clinical-data.tsv',sep='\t',header=TRUE,row.names=1)
row.names(clin_data) = str_replace_all(row.names(clin_data), "-", ".")
morph_counts = read.delim(file='morphological-counts.tsv',sep='\t',header=TRUE,row.names=1)


#PCA on unfiltered data
cpm_log_rna = cpm(rna_counts, log=TRUE)
pca_unclean = prcomp(t(cpm_log_rna),scale.=TRUE)
#summary(pca_unclean)
#plot(pca_unclean)

#There seem to be two groups: one with PC1 > 0 and one with PC1 < 0
colors = rep("black", length=nrow(cpm_log_rna))
colors[pca_unclean$x[,1] > 0] = "red"
colors[pca_unclean$x[,1] < 0] = "blue"
plot(pca_unclean$x[, 1], pca_unclean$x[, 2], pch = "*", xlab = "PC1", ylab = "PC2", col=colors)
#text(pca_unclean$x[, 1], pca_unclean$x[, 2], labels = colnames(cpm_log_rna))



#Group the clinical data in terms of these two groups
groups = rep(0, length=nrow(clin_data))
for (i in 1:ncol(rna_counts)) {
  if (pca_unclean$x[,1][i] > 0){
    groups[i] = 1
  } else{
      groups[i] = 2
  }
}
clin_data$expression_group = groups

#Write to file
write.xlsx(clin_data, file = "clinical-data-with-expression-group.xlsx")


# # BELOW I TRIED TO FILTER OUT THE LOWLY EXPRESSED GENES FIRST,
# # BUT UNFILTERED LOOKS MORE INTERESTING
# # Remove genes that are unexpressed or very lowly expressed in the samples
# # Used same cutoff as tutorial: based on median log2-transformed counts
# # per gene per million mapped reads (cpm in edgeR)
# median_log2_cpm = apply(cpm_log_rna, 1, median)
# expr_cutoff = -1
# 
# # This histogram shows our cutoff
# hist(median_log2_cpm)
# abline(v = expr_cutoff, col="red", lwd=3)
# 
# sum(median_log2_cpm > expr_cutoff)
# cleaned_rna_counts = rna_counts[median_log2_cpm > expr_cutoff,]
# cpm_log_cleaned = cpm(cleaned_rna_counts, log=TRUE)
# heatmap(cor(cpm_log_cleaned))
# 
# #PCA on filtered data
# pca = prcomp(t(cpm_log_cleaned),scale.=TRUE)
# plot(pca$x[, 1], pca$x[, 2], pch = "*", xlab = "PC1", ylab = "PC2")
# text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log_cleaned))
# summary(pca)

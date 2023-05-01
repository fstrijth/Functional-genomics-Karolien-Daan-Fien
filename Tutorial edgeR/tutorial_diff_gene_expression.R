library("edgeR")
library(stringr)
library("xlsx")
library("sjmisc")

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
pca = prcomp(t(cpm_log_rna),scale.=TRUE)
summary(pca)
#plot(pca)

#Group the clinical data in terms of these two groups
groups = rep(0, length=nrow(clin_data))
for (i in 1:ncol(rna_counts)) {
  if (pca$x[,1][i] > 0){
    groups[i] = 1
  } else{
      groups[i] = 2
  }
}
clin_data$expression_group = groups

#Also add a column that is 1 if any keyword associated to thyroiditis is present
thyroiditis = rep(0, length=nrow(clin_data))
for (i in 1:nrow(clin_data)) {
  if (str_contains(clin_data$SMPTHNTS[i], c("hashimoto", "thyroiditis", "goiter", "fibrosis", "goitre", "fibrous"), ignore.case=TRUE, logic="or")) {
    thyroiditis[i] = 1
  }
}
thyroiditis
clin_data$thyroiditis = thyroiditis

#There seem to be two groups: one with PC1 > 0 and one with PC1 < 0
colors = rep("black", length=nrow(cpm_log_rna))
colors[clin_data$AGE > 50] = "red"
colors[clin_data$AGE < 50] = "blue"
plot(pca$x[, 1], pca$x[, 2], pch = "*", xlab = "PC1", ylab = "PC2", col=colors)
#text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log_rna))


#Write to file
write.xlsx(clin_data, file = "clinical-data-with-expression-group.xlsx")


#PCA on morphological data
cpm_log_morph = cpm(morph_counts, log=TRUE)
pca = prcomp(t(cpm_log_morph),scale.=TRUE)
#summary(pca)
plot(pca$x[, 1], pca$x[, 2], pch = "*", xlab = "PC1", ylab = "PC2")


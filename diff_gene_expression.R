library("edgeR")

#Completely followed this:
#https://www.reneshbedre.com/blog/edger-tutorial.html

#Read file
rna_counts = read.delim(file='RNA-read-counts.tsv',sep='\t',header=TRUE,row.names=1)
#Remove description column
rna_counts = subset(rna_counts, select=-c(Description))
#Get column names
sample_info = colnames(rna_counts)
#Create DGEList data class
dge = DGEList(counts = rna_counts, group = factor(sample_info))
#Filter out genes with low counts
keep = filterByExpr(y = dge)
dge = dge[keep, ,keep.lib.sizes=FALSE]
#Normalization and effective library sizes
dge = calcNormFactors(object=dge)
#Model fitting and estimating dispersions
#dge = estimateDisp(y = dge) FOR MULTIPLE BATCHES?

#Testing for differential gene expression
et = exactTest(object=dge)


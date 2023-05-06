library("DESeq2")
library(ggplot2)
library(corrplot)

# Q2: clinical data vs. morphology
# 1) Compute systematically associations between clinical variables and morphological
#    cluster counts. The purpose is to compare the magnitude of the associations of the 
#    different variables with morphology.
# 2) Discuss the association with technical variables.
# 3) For non-technical variables, redo the analysis with adjustment for the confounding 
#    technical variables, if any is reported in Q2.2. Report and discuss significant 
#    associations.


#load data
clinical_data = read.delim(file='clinical-data.tsv',sep ='\t',header=TRUE,row.names=1)
#DTHHRDY is a categorical variable, set its type to character to that it can be viewed as such by deseq2
clinical_data$DTHHRDY <- as.character(clinical_data$DTHHRDY)
morph_counts = t(read.delim(file='morphological-counts.tsv', sep='\t', header=TRUE, row.names=1))
#Take out "Morphological.cluster" from column names for brevity's sake
row.names(morph_counts) = gsub("Mophological.cluster.", "", row.names(morph_counts))


#Automatically generates interesting plots in function of var
#All plots are saved under morph_plots/ and file names end in the variable name
#In case of a categorical variable, there are some extra plots that can be
#generated lile plotCounts and volcano plot. These are not interesting
#for the other variables because they are typically made in a 
#"cancer/no cancer" type situation, this can be an interesting style
#of analysis for e.g. postmortem/donor, but not for weight or height.

morph_analysis <- function(var, formula, categorical) {
  
  #construct DESEQDataset object
  #we don't normalize the counts since the tutorial from bioconductor
  #says deseq2 does not expect that: 
  #https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts
  
  #dds holds the read counts. The design formula expresses the variables
  #which will be used in modeling, it is used to estimate the dispersion
  #and to estimate the log2 fold changes of the model
  
  #We will not be pre-filtering low count clusters since we are only working
  #with 64 clusters instead of hundreds or thousands of genes
  #https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering
  dds = DESeqDataSetFromMatrix(countData = morph_counts, 
                               colData = clinical_data,
                               design = formula)
  #Run DESeq function
  dds = DESeq(dds)
  res = results(dds)
  #summary of the differential morphology expression
  summary(res)
  
  #Plotting adjusted p-values for each morphological cluster
  #We will go with the default cutoff value 0.1 for padj
  #https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#p-values-and-adjusted-p-values
  #The more dots below the line, the more clusters are significantly linked
  #to the chosen clinical/technical variable(s)
  png(paste("morph_plots/adjusted_P_values/adjusted_P_values_",var,".png", sep=""))
  plot(rownames(morph_counts), res$padj)
  abline(h=0.1, col="red")
  dev.off()
  
  #Normalized counts of cluster
  #Shows the log2 fold changes attributable to a given variable over the mean
  #of normalized counts for all the samples in the data set
  #https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#ma-plot
  png(paste("morph_plots/normalized_counts/normalized_counts_",var,".png", sep=""))
  plotMA(res)
  #TODO: bigger dots
  dev.off()
  
  #sort results by p-value
  res <- na.omit(res[order(res$padj),])
  
  #Save list of most significant clusters
  #cutoff: padj has to be smaller than 0.1
  sign_clusters = res[res$padj < 0.1,]
  write.csv(sign_clusters, paste("Morph_plots/cluster_lists/lowest_padj_",var,".csv",sep=""))
  
  #Plot the counts of the clusters with the 10 lowest adjusted p-values, one plot per cluster
  #https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#plot-counts
  if (categorical) {
    clusters = rownames(head(res,10))
    for (cluster in clusters) {
      png(paste("morph_plots/counts_by_",var,"/",cluster,".png", sep=""))
      plotCounts(dds, gene=cluster, intgroup=var, main=paste("Cluster ",cluster,sep=""))
      dev.off()
    }
  }
  
  #Volcano plot
  #Blue = padj < 0.1
  #Red = padj < 0.1 and fold change > 2
  #https://en.wikipedia.org/wiki/Volcano_plot_(statistics)
  if (categorical) {
    png(paste("morph_plots/volcano/volcano_",var,".png", sep=""))
    par(mfrow=c(1,1))
    with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
    with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
    with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
    dev.off()
  }
  
  #PCA
  #TODO???
}

#Q2.1: clinical variables
morph_analysis("AGE", ~ AGE, FALSE)
morph_analysis("HGHT", ~ HGHT, FALSE)
morph_analysis("WGHT", ~ WGHT, FALSE)

#Q2.2: Technical variables
morph_analysis("COHORT", ~ COHORT, TRUE)
morph_analysis("DTHHRDY", ~ DTHHRDY, TRUE)
morph_analysis("TRISCHD", ~TRISCHD, FALSE)

#Q2.3: Accounting for confounding variables
#TODO

#Plotting correlation
#TODO: only for relevant morphological clusters
png("test.png")
corrplot(cor(t(morph_counts)))
dev.off()



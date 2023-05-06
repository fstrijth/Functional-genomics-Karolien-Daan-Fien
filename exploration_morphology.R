library("DESeq2")
library(ggplot2)

# Q2: clinical data vs. morphology
# 1) Compute systematically associations between clinical variables and morphological
#    cluster counts. The purpose is to compare the magnitude of the associations of the 
#    different variables with morphology.
# 2) Discuss the association with technical variables.
# 3) For non-technical variables, redo the analysis with adjustment for the confounding 
#    technical variables, if any is reported in Q2.2. Report and discuss significant 
#    associations.


# Inspired by this tutorial:
#https://lashlock.github.io/compbio/R_presentation.html


#load data
clinical_data = read.delim(file='clinical-data.tsv',sep ='\t',header=TRUE,row.names=1)
#DTHHRDY is a categorical variable, set its type to character to that it can be viewed as such by deseq2
clinical_data$DTHHRDY <- as.character(clinical_data$DTHHRDY)
morph_counts = t(read.delim(file='morphological-counts.tsv', sep='\t', header=TRUE, row.names=1))
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
  dds = DESeqDataSetFromMatrix(countData = morph_counts, 
                               colData = clinical_data,
                               design = formula)
  #Run DESeq function
  dds = DESeq(dds)
  res = results(dds)
  #summary of the differential morphology expression
  summary(res)
  
  #Plotting adjusted p-values for each morphological cluster
  png(paste("morph_plots/adjusted_P_values/adjusted_P_values_",var,".png", sep=""))
  plot(rownames(morph_counts), res$padj)
  abline(h=0.1, col="red") #cutoff?
  dev.off()
  
  #Normalized counts of cluster
  png(paste("morph_plots/normalized_counts/normalized_counts_",var,".png", sep=""))
  plotMA(res)
  #TODO: bigger dots
  dev.off()
  
  #sort results by p-value
  res <- res[order(res$padj),]
  
  #TODO: save list of most significant clusters
  
  #Plot the counts of the clusters with the lowest adjusted p-values
  if (categorical) {
    clusters = rownames(head(res,10))
    for (cluster in clusters) {
      png(paste("morph_plots/counts_by_",var,"/",cluster,".png", sep=""))
      plotCounts(dds, gene=cluster, intgroup=var, main=paste("Cluster ",cluster,sep=""))
      dev.off()
    }
  }
  
  #Volcano plot
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



morph_analysis("COHORT", ~ COHORT, TRUE)
morph_analysis("DTHHRDY", ~ DTHHRDY, TRUE)
morph_analysis("AGE", ~ AGE, FALSE)
morph_analysis("HGHT", ~ HGHT, FALSE)
morph_analysis("WGHT", ~ WGHT, FALSE)
morph_analysis("TRISCHD", ~TRISCHD, FALSE)



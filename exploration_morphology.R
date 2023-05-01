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
#In case of a categorical variable, deseq2's plotCounts can be called, if you want that
#, set pltcnt to TRUE
morph_analysis <- function(var, formula, pltcnt) {
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
  png(paste("morph_plots/adjusted_P_values_",var,".png", sep=""))
  plot(rownames(morph_counts), res$padj)
  abline(h=0.1, col="red") #cutoff?
  dev.off()
  
  #Normalized counts of cluster
  png(paste("morph_plots/normalized_counts_",var,".png", sep=""))
  plotMA(res)
  dev.off()
  
  #sort results by p-value
  res <- res[order(res$padj),]

  
  #Plot the counts of the clusters with the lowest adjusted p-values
  if (pltcnt) {
    png(paste("morph_plots/counts_by_",var,".png", sep=""))
    par(mfrow=c(3,4))
    clusters = rownames(head(res,10))
    for (cluster in clusters) {
      plotCounts(dds, gene=cluster, intgroup=var)
    }
    dev.off()
  }
  
  #Volcano plot
  png(paste("morph_plots/volcano_",var,".png", sep=""))
  par(mfrow=c(1,1))
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
  with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  dev.off()
  
  #PCA
  #TODO???
}

morph_analysis("COHORT", ~ COHORT, TRUE)
morph_analysis("DTHHRDY", ~ DTHHRDY, TRUE)
morph_analysis("AGE", ~ AGE, FALSE)
morph_analysis("HGHT", ~ HGHT, FALSE)
morph_analysis("WGHT", ~ WGHT, FALSE)



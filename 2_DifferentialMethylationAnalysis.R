#!/usr/bin/env Rscript
# Author: Marta Espinosa Camarena
# email: mespinosa@imim.es
# Date: 23th Oct 2023
# Differential Methylation Analysis
# https://www.biostars.org/p/465230/

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(dplyr)
library(tidyr)
library(rtracklayer)
library(minfiData)
library(limma)
# library(RGChannelSet)
library(missMethyl)
library(DMRcate)
library(Gviz)
library(ggplot2)
library(wacolors)
library(edgeR)
library(gridExtra)

########################## reading cancer types ################################
# 450k paired patients bvalue files
inpath <- "/users/genomics/marta/TCGA_methylation/analysis/0_bvalues"
ctypes = list.files(inpath, pattern="^TCGA", full.names=F)
# dir.create("/users/genomics/marta/TCGA_methylation/plots/UCSC_Groups")

dfList <- list()  ## create empty list
plots_UCSC_Groups <- list()
deff.meth.annot.interesting.List <- list()
DeltaBvaluesList <- list()
################################################################################


for(f in 1:length(ctypes)) {
  print(ctypes[f])
  bval = readRDS(file.path(inpath,ctypes[f],"bval.RDS"))
  
  ##### compute differences of beta values (average tumor - average normal) ####
  DeltaBvalues = bval 
  DeltaBvalues$probeID = rownames(DeltaBvalues)
  DeltaBvalues = DeltaBvalues %>% pivot_longer(cols=-c(probeID), names_to = "sample", values_to = "bval") %>% 
    mutate(normal_tumor = case_when(grepl("tumor", sample) ~ "tumor",
                                    grepl("normal", sample) ~ "normal")) 
  
  normals = DeltaBvalues %>% subset(normal_tumor == "normal") %>% group_by(probeID) %>% summarise(average_normal = mean(bval))
  tumors = DeltaBvalues %>% subset(normal_tumor == "tumor") %>% group_by(probeID) %>% summarise(average_tumor = mean(bval))
  
  DeltaBvalues = merge(normals, tumors, by="probeID")
  DeltaBvalues$difference = DeltaBvalues$average_tumor - DeltaBvalues$average_normal
  DeltaBvalues = DeltaBvalues %>% mutate(hypo_hyper = case_when(difference > 0.2 ~ "hyper",
                                                                difference < -0.2 ~ "hypo"))
  DeltaBvaluesList[[ctypes[f]]] = DeltaBvalues
  ##############################################################################
  
  mval = readRDS(file.path("/users/genomics/marta/TCGA_methylation/analysis/1_Mvalues",ctypes[f],"mval.RDS"))
  
  samples = data.frame("sample" = colnames(mval))
  samples = samples %>% mutate(sample_type = case_when(grepl("normal", sample) ~ "normal",
                                                       grepl("tumor", sample) ~ "tumor"))
  
  #_____________ DMC analysis________________#
  design <- model.matrix(~ sample_type, data = samples)
  # fit the linear model 
  fit <- lmFit(mval, design)
  fit2 <- eBayes(fit)
  
  # extracting significantly methylated probes
  deff.meth = topTable(fit2, coef=ncol(design), sort.by="p",number = nrow(mval), adjust.method = "BY", p.value = 0.01)
  dfList[[ctypes[f]]] <- deff.meth
  
  #### Visualization
  #### plot the top 10 most significantly differentially methylated CpGs
  # par(mfrow=c(2,5))
  # sapply(rownames(deff.meth)[1:10], function(cpg){
  #   plotCpg(bval, cpg=cpg, pheno=samples$sample_type, ylab = "Beta values")
  # })
  
  ## making a volcano plot
  ## making dataset
  dat <- data.frame(foldchange = fit[["coefficients"]][,2], logPvalue =  -log10(fit2[["p.value"]][,2]))
  dat$threshold <- as.factor(abs(dat$foldchange) < 0.4)
  
  # dir.create("/users/genomics/marta/TCGA_methylation/plots/DMCAnalysis", recursive = T)
  ## Visualization
  # cols <- c("TRUE" = "grey", "FALSE" = "blue")
  # ggplot(data=dat, aes(x=foldchange, y = logPvalue, color=threshold)) +
  #   geom_point(size=1.2) +
  #   scale_colour_manual(values = cols) +
  #   geom_vline(xintercept = 0.4, colour="#990000", linetype="dashed") + 
  #   geom_vline(xintercept = - 0.4, colour="#990000", linetype="dashed") +
  #   theme(legend.position="none") +
  #   xlab("Fold Change") +
  #   ylab("-log10 p value") +
  #   ggtitle(ctypes[f]) +
  #   theme_bw() +
  #   theme(legend.position = "none")
  # ggsave(paste0("/users/genomics/marta/TCGA_methylation/plots/DMCAnalysis/",ctypes[f],"_volcanoPlot.pdf"), width=5.88, height=4.84)
  # ggsave(paste0("/users/genomics/marta/TCGA_methylation/plots/DMCAnalysis/",ctypes[f],"_volcanoPlot.png"), width=5.88, height=4.84)
}

# export
saveRDS(dfList, file = "/users/genomics/marta/TCGA_methylation/analysis/DE_CpGs.RDS")
saveRDS(deff.meth.annot.interesting.List, file = "/users/genomics/marta/TCGA_methylation/analysis/DE_CpGs_annotated.RDS")
saveRDS(DeltaBvaluesList, file = "/users/genomics/marta/TCGA_methylation/analysis/DeltaBetaValues.RDS")


## get intersection of genes with differentially methylated CpGs
# Step 1: Extract the unique gene names from each dataframe
gene_name_lists <- lapply(deff.meth.annot.interesting.List, function(df) unique(df$UCSC_RefGene_Name))

# Initialize a variable to store the maximum count of dataframes containing a gene name
max_count <- 0

# Initialize a variable to store the gene names that meet the condition
max_gene_names <- character(0)

# Loop through unique gene names and count how many dataframes contain each name
for (gene_name in unique(unlist(gene_name_lists))) {
  count <- sum(sapply(gene_name_lists, function(gene_list) gene_name %in% gene_list))
  if (count >= max_count) {
    # If the count is greater than or equal to the current maximum, update max_count and max_gene_names
    max_count <- count
    max_gene_names <- c(max_gene_names, gene_name)
  }
}

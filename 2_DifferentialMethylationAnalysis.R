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
# dir.create("/users/genomics/marta/TCGA_methylation/plots/Hyper_Hypo")

dfList <- list()  ## create empty list
plots_UCSC_Groups <- list()
deff.meth.annot.interesting.List <- list()
DeltaBvaluesList <- list()

################################################################################


for(f in 1:length(ctypes)) {
  print(ctypes[f])
  bval = readRDS(file.path(inpath,ctypes[f],"bval.RDS"))
  
  ##### compute differences of beta values tumor vs normal ####
  DeltaBvalues = bval 
  DeltaBvalues$probeID = rownames(DeltaBvalues)
  DeltaBvalues = DeltaBvalues %>% pivot_longer(cols=-c(probeID), names_to = "sample", values_to = "bval") %>% 
    mutate(normal_tumor = case_when(grepl("tumor", sample) ~ "tumor",
                                    grepl("normal", sample) ~ "normal")) 
  
  ### option A) Average Tumors - Average Normals
  normals = DeltaBvalues %>% subset(normal_tumor == "normal") %>% group_by(probeID) %>% summarise(average_normal = mean(bval))
  tumors = DeltaBvalues %>% subset(normal_tumor == "tumor") %>% group_by(probeID) %>% summarise(average_tumor = mean(bval))
  
  DeltaBvalues_A = merge(normals, tumors, by="probeID")
  DeltaBvalues_A$diff_bval = DeltaBvalues_A$average_tumor - DeltaBvalues_A$average_normal
  DeltaBvalues_A = DeltaBvalues_A %>% mutate(hypo_hyper = case_when(diff_bval > 0.2 ~ "Hypermethylated",
                                                                diff_bval < -0.2 ~ "Hypomethylated",
                                                                TRUE ~ "Intermediate"))
  table(DeltaBvalues_A$hypo_hyper)
  
  ### option B) Average(Tumors - Normals)
  DeltaBvalues_B = DeltaBvalues
  DeltaBvalues_B$patient = gsub("_.*","",DeltaBvalues_B$sample)
  
  DeltaBvalues_B = DeltaBvalues_B %>% select(-sample) %>% group_by(probeID,patient) %>% 
    summarize(diff_bval = mean(bval[normal_tumor == "tumor"]) - mean(bval[normal_tumor == "normal"])) %>% distinct()
  DeltaBvalues_B_final = DeltaBvalues_B %>% group_by(probeID) %>% summarise(mean_diff_bval = mean(diff_bval))
  DeltaBvalues_B_final = DeltaBvalues_B_final %>% mutate(hypo_hyper = case_when(mean_diff_bval > 0.2 ~ "Hypermethylated",
                                                                                mean_diff_bval < -0.2 ~ "Hypomethylated",
                                                                    TRUE ~ "Intermediate"))
  table(DeltaBvalues_B_final$hypo_hyper)
  
  
  ### option C) Tumors - Normals in 90% of the patients
  DeltaBvalues_C = DeltaBvalues_B
  num_patients = nrow(DeltaBvalues_C %>% ungroup() %>% select(patient) %>% distinct)
  
  DeltaBvalues_C_final = DeltaBvalues_C %>% group_by(probeID) %>%
    mutate(hypo_hyper = case_when(
      sum(diff_bval > 0.2) > num_patients*0.9 ~ "Hypermethylated",
      sum(diff_bval < -0.2) > num_patients*0.9 ~ "Hypomethylated",
      TRUE ~ "Intermediate"))
  
  DeltaBvalues_C_final = DeltaBvalues_C_final %>% group_by(probeID, hypo_hyper) %>% summarise(mean_diff_bval = mean(diff_bval)) %>% select(probeID, mean_diff_bval, hypo_hyper) %>% distinct()
  
  table(DeltaBvalues_C_final$hypo_hyper)
  DeltaBvaluesList[[ctypes[f]]] = DeltaBvalues_C_final
  
  ggplot(to_plot, aes(x=factor(hypo_hyper, level=c('Hypomethylated', 'Intermediate', 'Hypermethylated')), y=mean_diff_bval, fill=hypo_hyper)) +
    geom_boxplot()+
    labs(y="Mean(Tumor BValues - Normal BValues)",
         x="Level of methylation in >90% of the patients",
         title=paste0("Methylation level | ", ctypes[f])) +
    scale_color_manual("Methylation state", values=c("#2c7fb8","#feb24c","#bdbdbd")) +
    scale_fill_manual("Methylation state", values=c("#2c7fb8","#feb24c","#bdbdbd")) +
    theme_bw()  
  ggsave(paste0("/users/genomics/marta/TCGA_methylation/plots/Hyper_Hypo/boxplot_HyperHypoInt_",ctypes[f],".pdf"), width=6.18, height=4.52)
  ggsave(paste0("/users/genomics/marta/TCGA_methylation/plots/Hyper_Hypo/boxplot_HyperHypoInt_",ctypes[f],".png"), width=6.18, height=4.52)
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
saveRDS(DeltaBvaluesList, file = "/users/genomics/marta/TCGA_methylation/analysis/DeltaBetaValues.RDS")


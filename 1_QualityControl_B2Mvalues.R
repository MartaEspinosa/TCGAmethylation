#!/usr/bin/env Rscript
# Author: Marta Espinosa Camarena
# email: mespinosa@imim.es
# Date: 23th Oct 2023
# Quality control of beta-values and conversion to M values
# https://www.biostars.org/p/465230/

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(limma)
library(missMethyl)
library(DMRcate)
library(Gviz)
library(ggplot2)
library(RColorBrewer)
library(edgeR)
library(dplyr)

########################## reading cancer types ################################
# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# 450k paired patients bvalue files
inpath <- "/users/genomics/marta/TCGA_methylation/analysis/0_bvalues"
files = list.files(inpath, recursive=T, pattern = "bvalues_450Kpaired.csv", full.names=T)
ctypes = list.files(inpath, pattern="^TCGA", full.names=F)
################################################################################

#___________inspecting methylation data_______________#

summary = data.frame("initial_CpGs" = numeric(),
                "cancer_type" = factor(),
                "no_NAs" = numeric(),
                "noSNPs_MAF" = numeric(),
                "noMultiMapping" = numeric(),
                stringsAsFactors = F)

for(f in 1:length(files)) {
  print(ctypes[f])
  # read beta-values
  met = read.csv(files[f])
  rownames(met) = met$CpG
  met$CpG = NULL
  
  raw_num_cases = nrow(met)
  
  # chose those has no NA values in rows
  met = met[complete.cases(met), ]
  met %>% head
  noNA_num_cases = nrow(met)
  
  ## remove probes that match chromosomes X and Y 
  keep <- !(row.names(met) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
  table(keep)
  met <- met[keep, ]
  rm(keep) # remove no further needed probes
  
  ## remove SNPs overlapped probe
  table (is.na(ann450k$Probe_rs))
  # probes without snp
  no.snp.probe <- ann450k$Name[is.na(ann450k$Probe_rs)]
  snp.probe <- ann450k[!is.na(ann450k$Probe_rs), ]
  #snps with maf <= 0.05
  snp5.probe <- snp.probe$Name[snp.probe$Probe_maf <= 0.05]
  # filter met
  met <- met[row.names(met) %in% c(no.snp.probe, snp5.probe), ]
  rm(no.snp.probe, snp.probe, snp5.probe) #remove no-further needed dataset
  noSNPs_MAF_num_cases = nrow(met)
  
  ## Removing probes that have been demonstrated to map to multiple places in the genome.
  # list adapted from https://www.tandfonline.com/doi/full/10.4161/epi.23470
  crs.reac <- read.csv("/users/genomics/marta/TCGA_methylation/cross_reactive_probe.chen2013.csv")
  crs.reac <- crs.reac$TargetID[-1]
  
  # filtre met
  met <- met[ -which(row.names(met) %in% crs.reac), ]
  filtered_num_cases = nrow(met)
  bval <- met
  
  tmp = data.frame("initial_CpGs" = raw_num_cases,
                   "cancer_type" = ctypes[f],
                   "no_NAs" = noNA_num_cases,
                   "noSNPs_MAF" = noSNPs_MAF_num_cases,
                   "noMultiMapping" = filtered_num_cases,
                   stringsAsFactors = F)
  summary = rbind(summary, tmp)
  
  ## converting beta values to m_values
  ## m = log2(beta/1-beta)
  mval <- t(apply(met, 1, function(x) log2(x/(1-x))))
  
  #______________saving/loading_____________________#
  # save data sets
  dir.create(file.path("/users/genomics/marta/TCGA_methylation/analysis/1_Mvalues",ctypes[f]), recursive = T)
  
  saveRDS(mval, file = file.path("/users/genomics/marta/TCGA_methylation/analysis/1_Mvalues",ctypes[f],"mval.RDS"))
  write.csv(mval, file = file.path("/users/genomics/marta/TCGA_methylation/analysis/1_Mvalues",ctypes[f],"Mvalues_filtered.csv"), quote=F)
  
  saveRDS(bval, file = file.path(inpath,ctypes[f],"bval.RDS"))
  write.csv(bval, file = file.path(inpath,ctypes[f],"bvalues_filtered.csv"), quote=F)

}
summary %>% head
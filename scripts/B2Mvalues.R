#!/usr/bin/env Rscript
# Author: Marta Espinosa Camarena
# email: mespinosa@imim.es
# Date: 19th Oct 2023
# This scripts passes from B-values to M-values

# Run like this:
# Rscript B2Mvalues.R $inpath $outpath

######################## libraries to be loaded #################################
library(dplyr)
library(tidyr)
library(data.table)
library(lumi)
options(warn=-1)

################################################################################

############################## reading inputs ##################################
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("Please provide the path until the folder with the b values and the output folder", call.=FALSE)
}
################################################################################

# Read arguments
inpath <- args[1]
outpath <- args[2]

# inpath <- "/users/genomics/marta/TCGA_methylation/analysis/0_bvalues/TCGA-BLCA"
# outpath <- "/users/genomics/marta/TCGA_methylation/analysis/1_mvalues/TCGA-BLCA"
dir.create(outpath, recursive = T)

################################################
### Read input
################################################
bvalue_files = list.files(inpath, full.names=T)
patientIDs = list.files(inpath, full.names=F)
patientIDs = gsub("_betavalues.txt","",patientIDs)

for(i in 1:length(bvalue_files)) {
  # name of outfile
  outfile = paste0(outpath,"/",patientIDs[i],"_mvalues.txt")
  
  # read beta-values
  bvalue = as.data.frame(fread(bvalue_files[i]))
  row.names(bvalue) = bvalue$CpG
  bvalue$CpG = NULL
  bvalue = as.matrix(bvalue)
  
  # convert to m-values
  mvalue = as.data.frame(beta2m(bvalue))
  colnames(mvalue) = "Mvalues"
  mvalue$CpG <- rownames(mvalue)
  # reorder
  mvalue <- mvalue[, c("CpG", "Mvalues")]
  
  # save output Mvalues
  write.table(mvalue, outfile, row.names = F, quote=F, sep="\t")
}

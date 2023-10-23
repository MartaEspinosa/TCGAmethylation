#!/usr/bin/env Rscript
# Author: Marta Espinosa Camarena
# email: mespinosa@imim.es
# Date: 19th Oct 2023
# This scripts passes from individual files with B or M values to merged ones. If a different number of CpGs are found, NAs are used to fill blank spaces

# Run like this:
# Rscript onFile_CancerType.R $inpath

######################## libraries to be loaded #################################

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(data.table))
options(warn=-1)

################################################################################

############################## reading inputs ##################################
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("Please provide the path until the folder with the values", call.=FALSE)
}
################################################################################

# Read arguments
inpath <- args[1]
# inpath <- "/users/genomics/marta/TCGA_methylation/analysis/0_bvalues/TCGA-BLCA"

if (grepl("mvalues", inpath)) {
  values <- "mvalues"
} else if (grepl("bvalues", inpath)) {
  values <- "bvalues"
}

################################################
### Read input
################################################
files = list.files(inpath, pattern="^TCGA", full.names=T)
patientIDs = list.files(inpath, pattern="^TCGA", full.names=F)
if (values == "mvalues") {
  patientIDs = gsub("_mvalues.txt","",patientIDs)
} else if (values == "bvalues") {
  patientIDs = gsub("_betavalues.txt","",patientIDs)
}

df = data.frame("CpG" = character(),
                stringsAsFactors = F)

for(i in 1:length(files)) {
  # read file
  input = as.data.frame(fread(files[i]))
  names(input)[2] = patientIDs[i]
  
  # get one file where each column has the Bvalues per patient for all the CpGs. Keep NAs 
  df = merge(df, input, by="CpG", all.x = T, all.y = T)
  
}
df %>% head

write.csv(df, paste0(inpath,"/",values,".csv"), row.names =F, quote=F)

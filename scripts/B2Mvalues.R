#!/usr/bin/env Rscript
# Author: Marta Espinosa Camarena
# email: mespinosa@imim.es
# Date: 19th Oct 2023
# This scripts passes from B-values to M-values

# Run like this:
# Rscript B2Mvalues.R 

######################## libraries to be loaded #################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(data.table)
library(lumi)
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

################################################
### Read input
################################################
bvalue_files = list.file(inpath, all.files=T, full.names=T)
patientIDs = list.file(bvalue_files, all.files=T, full.names=F)
patientIDs = gsub("_betavalues.txt","",patientIDs)
patientIDs

for(i in 1:len(bvalue_files)) {
  # read beta-values
  bvalue = fread(bvalue_files[i])
  # convert to m-values
  mvalue = beta2m(bvalue)
  print(mvalue)
}

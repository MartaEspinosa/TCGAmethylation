#!/usr/bin/env Rscript
# Author: Marta Espinosa Camarena
# email: mespinosa@imim.es
# Date: 23th Oct 2023
# Quality control of beta-values

# Run like this:
# Rscript 
library(tidyr)
library(dplyr)
library(ggplot2)

########################## reading cancer types ################################
inpath <- "/users/genomics/marta/TCGA_methylation/analysis/0_bvalues"
ctypes_path = list.files(inpath, full.names=T)
ctypes = list.files(inpath)

patients_450 = read.csv("/users/genomics/marta/TCGA_methylation/results/patients_450K.csv")
patients_450 = patients_450$case_ID
patients_450 = unique(patients_450)
################################################################################

################################ betavalues ####################################
# df = data.frame("Cpg" = character(),
#                 "sample" = character(),
#                 "betavalue" = numeric(),
#                 "cancer_type" = factor(),
#                 "normal_tumor" = factor(),
#                 stringsAsFactors = F)
# for(c in 1:length(ctypes)) {
#   input = as.data.frame(fread(file.path(ctypes_path[c],"bvalues.csv"),sep=","))
#   
#   long = input %>% pivot_longer(cols=-c("CpG"), names_to = "sample", values_to = "betavalue") %>%
#     mutate("cancer_type" = ctypes[c]) %>%
#     mutate("normal_tumor" = case_when(grepl("tumor", sample) ~ "tumor",
#                                       grepl("normal", sample) ~ "normal"))
#   df = rbind(df, long)
# }
# df %>% head
# saveRDS(df, file=file.path(inpath,"merged_bvalues_longdf.rds"))

readRDS(file=file.path(inpath,"merged_bvalues_longdf.rds"))
ggplot(df, aes(x=betavalue, color=normal_tumor)) +
  geom_density()+
  facet_wrap(~ cancer_type, scales="free")
################################################################################

######################### sample quality control ###############################
########## input merged bvalues #######
merged_bvalues = list.files(inpath, recursive=T, pattern = "bvalues.csv", full.names=T)

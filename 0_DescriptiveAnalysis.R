#!/usr/bin/env Rscript
# Author: Marta Espinosa Camarena
# email: mespinosa@imim.es
# Date: 24th Oct 2023
# Descriptive analysis of paired methylation data from TCGA

library(dplyr)
library(ggplot2)
library(VennDiagram)

############################## input data ######################################
clinical = "/datasets/marta/TCGA/MET_matched_cancers_biospecimen/clinical_merged_modified_bothRNAMET.tsv"
################################################################################

####### reading data ####### 
metadata = read.csv(clinical)
metadata %>% head
metadata = metadata %>% 
  select(File.Name, Project.ID, caseID, Sample.Type, age_at_index, ethnicity, gender, race, vital_status, ajcc_pathologic_stage, primary_diagnosis, prior_malignancy, prior_treatment, treatment_type)
metadata$age_at_index = as.numeric(metadata$age_at_index)

####### Gender ####### 
table(metadata$gender, metadata$Project.ID)
gender = metadata %>% select(Project.ID, caseID, gender) %>% distinct()
ggplot(gender, aes(x=Project.ID, fill=gender)) +
  geom_bar(position = "fill") +
  scale_fill_manual("", values=c("#f1a340","#998ec3"), labels=c("Female","Male")) +
  labs(title="Gender balance",
       y="") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        legend.position = "top")

####### VitalStatus ####### 
table(metadata$vital_status, metadata$Project.ID)
vital_status = metadata %>% select(Project.ID, caseID, vital_status) %>% distinct()
ggplot(vital_status, aes(x=Project.ID, fill=vital_status)) +
  geom_bar(position = "fill") +
  scale_fill_manual("", values=c("#1a9850","#d73027","#999999"), labels=c("Alive","Dead","Not Reported")) +
  labs(title="Vital status",
       y="") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        legend.position = "top")

####### Age ####### 
age = metadata %>% select(Project.ID, caseID, age_at_index) %>% distinct()
age_summary = age %>% group_by(Project.ID) %>% 
  mutate(mean = mean(age_at_index),
         median = median(age_at_index),
         max = max(age_at_index),
         min = min(age_at_index)) %>% 
  select(-c(caseID, age_at_index)) %>%
  distinct()
age_summary

ggplot(age_summary, aes(x=Project.ID, y=mean)) +
  geom_point(size=5, color="#4d9221") +
  geom_segment(aes(x=Project.ID, xend=Project.ID, y=0, yend=mean), color="#4d9221") +
  labs(title="Age distribution",
       y="Age (mean)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90)) 

ggplot(age, aes(x=Project.ID, y=age_at_index, fill=Project.ID, color=Project.ID)) +
  geom_violin() +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange", color="black", size=0.25) +
  labs(title="Age distribution",
       y="") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90)) 

####### Ethnicity | Race ####### 
table(metadata$ethnicity, metadata$Project.ID)
table(metadata$race, metadata$Project.ID)

####### Prior Malignancy | Prior Treatment ####### 
table(metadata$prior_malignancy, metadata$Project.ID)
table(metadata$prior_treatment, metadata$Project.ID)

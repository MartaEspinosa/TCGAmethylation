#!/usr/bin/env Rscript
# Author: Marta Espinosa Camarena
# email: mespinosa@imim.es
# Date: 6th Nov 2023
# Annotation of DMC to Functional Genomic Regions. Result from script 2_DifferentialMethylationAnalysis

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(dplyr)
library(tidyr)
library(rtracklayer)
library(minfiData)
library(limma)
library(missMethyl)
library(DMRcate)
library(Gviz)
library(ggplot2)
library(wacolors)
library(edgeR)
library(gridExtra)

################################ importing #####################################
biomart = read.csv("/genomics/users/marta/genomes/transcript_gene.txt")
biomart = biomart %>% select(Gene.name, Gene.type)
names(biomart) = c("gene_name","gene_type")

dfList = readRDS("/users/genomics/marta/TCGA_methylation/analysis/DE_CpGs.RDS")
DeltaBvaluesList = readRDS("/users/genomics/marta/TCGA_methylation/analysis/DeltaBetaValues.RDS")

ctypes <- names(dfList)
DMC_manual_genes_interestingList = list()
DeltaBvalues_annot_biotypeList = list()
# dir.create("/users/genomics/marta/TCGA_methylation/plots/UCSC_Groups")
################################################################################


############################## BASIC ANNOTATION ################################
chain=import.chain("/genomics/users/marta/genomes/hg19ToHg38.over.chain") # file downloaded from UCSC

rg <- get(data("RGsetEx")) 
rg %>% head # RGsetEx
man <- as.data.frame(getAnnotation(rg))
colnames(man)

## liftover to hg38
annRanges= makeGRangesFromDataFrame(man[,1:4],
                                    keep.extra.columns=T,
                                    ignore.strand=FALSE,
                                    seqinfo=NULL,
                                    seqnames.field=c("seqnames", "seqname",
                                                     "chromosome", "chrom",
                                                     "chr", "chromosome_name",
                                                     "seqid"),
                                    start.field="pos",
                                    end.field="pos",
                                    strand.field="strand",
                                    starts.in.df.are.0based=FALSE)
hg38Locs=liftOver(annRanges,chain)

hg38LocsDF=data.frame(hg38Locs)
rownames(hg38LocsDF)=hg38LocsDF$group_name
pos38=start(unlist(hg38Locs))
man=data.frame(man,"pos.hg19"=man$pos)
man$pos=rep(NA,dim(man)[1])
man[hg38LocsDF$Name,"pos"]=hg38LocsDF[,"start"]
man %>% head

## split GeneNames and GeneGroups
as.data.frame(table(man$UCSC_RefGene_Name==""))
gene.regionv <- unlist(strsplit(man$UCSC_RefGene_Group, ";"))
unique(gene.regionv)

man_wrt_CGI = man %>% select(chr, Name, Relation_to_Island, Islands_Name)
names(man_wrt_CGI)[2] = "probeID"

man_wrt_genes = man %>% select(chr, pos, strand, Name, UCSC_RefGene_Name, UCSC_RefGene_Group)
names(man_wrt_genes)[5] = "gene_name"
names(man_wrt_genes)[4] = "probeID"

man_wrt_genes %>% head
man_wrt_genes = man_wrt_genes %>%
  separate_rows(gene_name, UCSC_RefGene_Group, sep=";")

man_wrt_genes_biomart = merge(man_wrt_genes, biomart, by="gene_name", all.x=T)
man_wrt_genes_biomart = man_wrt_genes_biomart %>% distinct()
man_wrt_genes_biomart %>% head ## too few lncRNAs

################################################################################


for(f in 1:length(ctypes)) {
  print(ctypes[f])
  #_____________ DMC analysis________________#
  ## input data
  deff.meth = dfList[[ctypes[f]]]
  deff.meth$probeID = rownames(deff.meth)
  
  DeltaBvalues = DeltaBvaluesList[[ctypes[f]]]
  
  DMC = merge(deff.meth, DeltaBvalues, by="probeID")
  if(nrow(DMC) > 0) {
    # merge with the manual (CGI-related)
    DMC_manual_CGI = merge(DMC, man_wrt_CGI, by.x="probeID", all.x=T)
    DMC_manual_CGI = DMC_manual_CGI %>% distinct()
    
    ggplot(DMC_manual_CGI %>% subset(hypo_hyper != "Intermediate") %>% subset(Relation_to_Island != "OpenSea"), aes(x=Relation_to_Island, fill=hypo_hyper)) +
      geom_bar() +
      # scale_y_continuous(trans="log2") +
      scale_fill_manual("Methylation state", values=c("#2c7fb8","#feb24c")) +
      labs(x="",
           y="Number of CpGs",
           title=paste0("Distribution of DMCs in different regions related to CGIs | ",ctypes[f])) +
      theme_bw() 
    ggsave(paste0("/users/genomics/marta/TCGA_methylation/plots/Hyper_Hypo/distr_DMCs_CGIregions_",ctypes[f],".png"), width=6.18, height=2.84)
    ggsave(paste0("/users/genomics/marta/TCGA_methylation/plots/Hyper_Hypo/distr_DMCs_CGIregions_",ctypes[f],".pdf"), width=6.18, height=2.84)
    
    # merge with the manual (gene-related)
    DMC_manual_genes = merge(DMC, man_wrt_genes_biomart, by.x="probeID", all.x=T)
    DMC_manual_genes = DMC_manual_genes %>% distinct()  
    
    # ggplot(DMC_manual_genes %>% subset(hypo_hyper != "Intermediate"), aes(x=UCSC_RefGene_Group, fill=hypo_hyper)) +
    #   geom_bar() +
    #   scale_fill_manual("Methylation state", values=c("#2c7fb8","#feb24c")) +   
    #   theme_bw()
    
    # subset those that are in gene types we are interested in
    DMC_manual_genes_interesting = DMC_manual_genes %>% subset(gene_type == "lncRNA" | gene_type == "processed_pseudogene" | gene_type == "protein_coding")
    DMC_manual_genes_interesting = DMC_manual_genes_interesting %>% mutate(biotype = case_when(gene_type == "processed_pseudogene" ~ "lncRNA",TRUE ~ gene_type))
    DMC_manual_genes_interesting = DMC_manual_genes_interesting %>%
      mutate(promoter = case_when(UCSC_RefGene_Group == "5'UTR" ~ "TRUE",
                                  UCSC_RefGene_Group == "TSS1500" ~ "TRUE",
                                  UCSC_RefGene_Group == "TSS200" ~ "TRUE",
                                  UCSC_RefGene_Group == "Body" ~ "FALSE",
                                  UCSC_RefGene_Group == "3'UTR" ~ "FALSE",
                                  UCSC_RefGene_Group == "1stExon" ~ "FALSE"))
    DMC_manual_genes_interestingList[[ctypes[f]]] = DMC_manual_genes_interesting
      
    ### plots UCSC Groups
    p1 <- ggplot(DMC_manual_genes_interesting, aes(x=biotype, fill=UCSC_RefGene_Group)) +
      geom_bar(position="fill") +
      labs(x="Biotype",
           y="",
           title=paste0("CpG position | ",ctypes[f])) +
      scale_fill_wa_d(wacolors$sound_sunset, reverse=TRUE) +
      theme_bw()
    ggsave(paste0("/users/genomics/marta/TCGA_methylation/plots/UCSC_Groups/UCSC_Groups_",ctypes[f],".png"), plot=p1, width=5.02, height=5.02)
    ggsave(paste0("/users/genomics/marta/TCGA_methylation/plots/UCSC_Groups/UCSC_Groups_",ctypes[f],".pdf"), plot=p1, width=5.02, height=5.02)
    
    p2 <- ggplot(DMC_manual_genes_interesting %>% subset(hypo_hyper != "Intermediate"), aes(x=UCSC_RefGene_Group, fill=hypo_hyper)) +
      geom_bar() +
      # scale_y_continuous(trans="log2") +
      labs(y="Number of CpGs",
           x="",
           title=paste0("Distribution of DMCs across gene regions | ",ctypes[f]))+
      scale_fill_manual("Methylation state", values=c("#2c7fb8","#feb24c")) +   
      theme_bw() +
      theme(legend.position = "top") +
      facet_wrap(~ biotype, scales="free_y", ncol=1)
    ggsave(paste0("/users/genomics/marta/TCGA_methylation/plots/Hyper_Hypo/distr_DMCs_GeneRegions_",ctypes[f],".png"), plot=p2, width=5.02, height=5.02)
    ggsave(paste0("/users/genomics/marta/TCGA_methylation/plots/Hyper_Hypo/distr_DMCs_GeneRegions_",ctypes[f],".pdf"), plot=p2, width=5.02, height=5.02)
      
    # how many of these are annotated in a lncRNA?
    DMC_manual_genes_lncRNA = DMC_manual_genes_interesting %>% subset(biotype == "lncRNA")
    length(unique(DMC_manual_genes_lncRNA$gene_name))
    # how many of these are annotated in a protein coding?
    DMC_manual_genes_PCGs = DMC_manual_genes_interesting %>% subset(gene_type == "protein_coding")
    length(unique(DMC_manual_genes_PCGs$gene_name))
  
    print(table(DMC_manual_genes_interesting$biotype, DMC_manual_genes_interesting$hypo_hyper))
    
    p3 <- ggplot(DMC_manual_genes_interesting %>% subset(hypo_hyper != "Intermediate"), aes(x=biotype, fill=hypo_hyper)) +
      geom_bar(position="fill") +
      labs(y="",
           x="",
           title=paste0("Hypo-Hyper | ",ctypes[f],"\n(Tumor Samples wrt Normal Adjacent Samples) ")) +
      scale_fill_manual("Difference Bvalues", values=c("#2c7fb8","#feb24c")) +
      theme_bw() +
      theme(legend.position = "top")
    ggsave(paste0("/users/genomics/marta/TCGA_methylation/plots/Hyper_Hypo/HyperHypo_biotype_",ctypes[f],".png"), plot=p3, width=4.45, height=5.02)
    ggsave(paste0("/users/genomics/marta/TCGA_methylation/plots/Hyper_Hypo/HyperHypo_biotype_",ctypes[f],".pdf"), plot=p3, width=4.45, height=5.02)
  }
}


## get intersection of genes with differentially methylated CpGs
# Step 1: Extract the unique gene names from each dataframe
gene_name_lists <- lapply(DMC_manual_genes_interestingList, function(df) unique(df$gene_name))

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

biomart %>% subset(gene_name %in% max_gene_names) %>% distinct()


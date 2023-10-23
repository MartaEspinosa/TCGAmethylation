# TCGAmethylation

## 0.preprocessing

**STEP 1**

- Create a cancer-type file with
  - Project ID
  - Case ID
  - Normal ID
  - Tumor ID
  - Gender
  - Age


**STEP 2**
Generate 1 file per cancer type where rows are CpGs and columns are beta-values corresponding to each sample.

## 1. Quality Control

The `EASIER` package is used.

- Remove NAs
- Exclude CpGs that meet exclusion conditions

A) Control and non-CpG probes + probes targeting CpGs in sexual chromosomes
  Control probes (“control_probes”): technical control probes that do not correspond
  to CpGs, such as bisulfite conversion I, bisulfite conversion II, extension, hybridization
  and negative. Classified as “rs” in the filtering variable named “probeType”;

  - Non-cpg probes (“noncpg_probes”): probes targeting non-cpg sites classified as
    “ch” in the filtering variable named “probeType”;

  - Sex chromosomes (“Sex”): probes targeting sex chromosomes. Select this to avoid
    misleading results due to differences in sex-chromosome dosage on the human methylome.
    Filtering variable “Sex”;

B) Probes with hybridizing problems
  - Poor mapping probes (“MASK_mapping”): Probes that have poor quality mapping
  to the target genomic location as indicated in the array’s manifest file based on genome
  build GRCh37 and GRCh38 (for example due to the presence of INDELs (Insertion–
  deletion mutations present in the genome);

  - Cross-hybridising probes (“MASK_sub30”): The sequence of the last 30bp at the
  3’ end of the probe is non-unique (problematic because the beta value of such probes
  is more likely to represent a combination of multiple sites and not the level of initially
  targeted CpG sites); (Zhou, Laird, and Shen 2017) recommend 30bp, but in this script
  there is the possibility to adapt this to probes with non-unique 25bp, or 35bp, or 40bp,
  or 45bp 3’-subsequences (“MASK_sub25”, “MASK_sub40”, “MASK_sub45”).

C) Probes affected by the presence of SNPs
  - “MASK_extBase”: Probes with a SNP altering the CpG dinucleotide sequence context
  and hence the ability of target cytosines to be methylated (regardless of the MAF);
  - “MASK_typeINextBaseSwitch”: Probes with a SNP in the extension base that
  causes a colour channel switch from the official annotation (regardless of the MAF);
  - “MASK_snp5_GMAF1p”: probes with SNPs at the last 5bp of the 3’ end of the
  probe, with an average minor allele frequency (MAF) >1%, by ethnic group;
  - “MASK_snp5_common”: probes with SNPs at the last 5bp of the 3’ end of the
  probe, with any average minor allele frequency (MAF) (can be <1%), by ethnic group;

D) Probes giving inconsistent results between arrays
  - “Unrel_450_EPIC_blood”: These are probes that are known to yield different results
  for the 450K and EPIC array in BLOOD, suggesting that results are unreliable for at
  least one of these arrays. CpGs based on (Solomon et al. 2018).

  - “Unrel_450_EPIC_pla_restrict” or “Unrel_450_EPIC_pla”: These are probes
  that are known to yield different results for the 450K and EPIC array in PLACENTA,
  suggesting that results are unreliable for at least one of these arrays. CpGs based
  on(Fernandez-Jimenez et al. 2019) .

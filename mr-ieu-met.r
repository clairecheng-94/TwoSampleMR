#LOAD PACKAGES 
library(TwoSampleMR)
library(data.table)

#Read table 
#Make sure names dictate which ones are exposure and outcome data 
#outcome
df <- fread('/home/cfc85413/PUFAS/ieu.out.tsv', header=TRUE, sep="\t")

#Omega 3 data table
dv<- fread('/home/cfc85413/PUFAS/CConvertP2.out.tsv', header=TRUE, sep="\t") 

#Reading Exposure data 
exposireO3<-read_exposure_data(filename ='/home/cfc85413/PUFAS/CConvertP2.out.tsv',
  sep = "\t" ,snp_col = "SNP",  
  beta_col ="BETA", 
  se_col = "SE", 
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "AF", 
  pval_col = "P",
  samplesize_col = 114999) 
  
outcome_dat <- read_outcome_data(snps = exposireO3$SNP,filename = '/home/cfc85413/PUFAS/ieu.out.tsv', sep = "\t" ,snp_col = "SNP",   
  beta_col ="BETA",  
  se_col = "SE", 
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "AF", 
  pval_col = "P",
  samplesize_col = 82315)
  
 #Harmonising the data
 res<-harmonise_data(exposireO3, outcome_dat)

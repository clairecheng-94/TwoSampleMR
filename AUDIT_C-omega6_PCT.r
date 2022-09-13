#Switch into the interact node
interact -c 4 --mem 50G

#go to scratch directory
cd /scratch/cfc85413/PUFAS

#Load R
ml R/4.1.3-foss-2020b
R

#Load the libraries
library(TwoSampleMR)
library(data.table)
library (dplyr) #make sure to load this one especially 

#choose either one of the options: fread or read_exposure_data
#read the outcome data
outcome_dat <- fread('')

#read the exposure data
exposure_dat<- fread('/scratch/cfc85413/PUFAS/UKB_Omega_6_pct.a1effect.munge.rmInDels.uniq.tsv.gz')

#read the exposure data
#filename can be the variant that has already been set 
exposure03<-read_exposure_data(filename ='/scratch/cfc85413/PUFAS/UKB_Omega_6_pct.a1effect.munge.rmInDels.uniq.tsv.gz', clump=FALSE, sep = "\t" , snp_col = "SNP", beta_col ="BETA",  se_col = "SE",effect_allele_col = "A1",  other_allele_col = "A2", eaf_col = "FRQ",pval_col = "P", chr_col="CHR", pos_col="BP")

#Filter out the p value by lower than 5e-8 (GWAS significant threshold) 
sigificant_exposure <- filter(exposure03,pval.exposure <= 5e-8) 


#Clump data: read the exposure data 
clump_dat <- clump_data(sigificant_exposure,  clump_kb = 10000, clump_r2 = 0.001,  clump_p1 = 5e-8,  clump_p2 = 5e-8, pop= "EUR")
  

#Using clumped data in the read_outcome_data function, this file does not have allele frequencey, still proceed forth 
outcome_dat <- read_outcome_data(snps = clump_dat$SNP,filename = '/scratch/cfc85413/PUFAS/AUDIT_C_30336701.a1effect.munge.rmInDels.uniq.tsv.gz',  sep="\t", snp_col= "SNP",  beta_col ="BETA", se_col = "SE", effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P", chr_col="CHR",pos_col= "BP")


#harmonize data 
res<-harmonise_data(clump_dat, outcome_dat) 

#Yitang's removal of genetic instruments, this filters out the new column of 'MR_Keep', you want to keep the TRUE data 
res_true<-filter(res, (mr_keep.exposure + mr_keep + mr_keep.outcome) > 0)
             
#Rename the exposure ID and the outcome ID to the Omega 6 and AUDI_C    
names(res_TRUE1)[names(res_TRUE1) == 'id.exposure']<- "Exposure"
names(res_TRUE1)[names(res_TRUE1)== 'id.outcome'] <- "Outcome"

#Do Mr





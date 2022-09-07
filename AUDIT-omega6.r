#Switch into the interact node
interact -c 4 --mem50G

#go to scratch directory
cd /scratch/cfc85413/PUFAS

#Load R
ml R/4.1.3-foss-2020b
R

#Load the libraries
library(TwoSampleMR)
library(data.table)

#read the outcome data
outcome_dat <- fread('/scratch/cfc85413/PUFAS/AUDIT_C_30336701.a1effect.munge.rmInDels.uniq.tsv.gz')

#read the exposure data
exposure_dat<- fread('/scratch/cfc85413/PUFAS/UKB_Omega_6_pct.a1effect.munge.rmInDels.uniq.tsv.gz')

#read the exposure data
exposure03<-read_exposure_data(filename ='/scratch/cfc85413/PUFAS/UKB_Omega_6_pct.a1effect.munge.rmInDels.uniq.tsv.gz',
+ sep = "\t" ,snp_col = "SNP",  
+   beta_col ="BETA", 
+   se_col = "SE", 
+ effect_allele_col = "A1",
+ other_allele_col = "A2",
+ eaf_col = "FRQ",
+ pval_col = "P",
+ samplesize_col = 1052486)

#Clump data
 clump_data<-fread('UKB_Omega_6_pct.a1effect.munge.rmInDels.uniq.tsv.gz', header=TRUE, sep="\t")

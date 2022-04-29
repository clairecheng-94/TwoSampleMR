#LOAD PACKAGES 
library(TwoSampleMR)
library(data.table)

#Read table 
#Make sure names dictate which ones are exposure and outcome data 
#outcome
outcome_dat <- fread('/home/cfc85413/PUFAS/ieu.out.tsv', header=TRUE, sep="\t")

#Omega 3 data table
exposure_dat<- fread('/home/cfc85413/PUFAS/CConvertP2.out.tsv', header=TRUE, sep="\t") 



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
  
#Clump Data
clump_data<-fread('/home/cfc85413/PUFAS/Omega_3_extract_exp_data_clump.txt', header=TRUE, sep="\t")


#Used clump data in the read_outcome_data 
outcome_dat <- read_outcome_data(snps = clump_data$SNP,filename = '/home/cfc85413/PUFAS/ieu.out.tsv', sep = "\t" ,snp_col = "SNP",   
  beta_col ="BETA",  
  se_col = "SE", 
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "AF", 
  pval_col = "P",
  samplesize_col = 82315)
  
 #Harmonising the data
 res<-harmonise_data(exposireO3, outcome_dat)

#sensitivity analysis 
 mr_heterogeneity(res)

mr_pleiotropy_test(res)
res_loo<-mr_leaveoneout(res)
head(res_loo)

funnel<-mr_singlesnp(
 res,
 parameters = default_parameters(),
 single_method = "mr_wald_ratio",
 all_method = c("mr_ivw", "mr_egger_regression")
 )

mr_funnel_plot(funnel)




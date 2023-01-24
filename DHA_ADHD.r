## Packages that are loaded 
 # version for TwoSample MR is 0.5.6
 # version for R package  is R/4.2.1-foss-2020b
library(TwoSampleMR) 
library (dplyr) 
library(data.table)

## Read Exposure data 
 # path to this file /scratch/cfc85413/PUFAS/PSYCH/UKB_DHA.a1effect.munge.rmInDels.uniq.tsv.gz
exposure_dat<-read_exposure_data(filename ='/scratch/cfc85413/PUFAS/PSYCH/UKB_DHA.a1effect.munge.rmInDels.uniq.tsv.gz', clump=FALSE, sep = "\t" , snp_col = "SNP", beta_col ="BETA",  se_col = "SE",effect_allele_col = "A1",  other_allele_col = "A2", eaf_col = "FRQ",pval_col = "P", chr_col="CHR", pos_col="BP")

## Filter 
sigificant_exposure <- exposure_dat[exposure_dat$pval.exposure<=5e-8,]

## Clumnping
clump_dat <- clump_data(sigificant_exposure,  clump_kb = 10000, clump_r2 = 0.001,  clump_p1 = 5e-8,  clump_p2 = 5e-8, pop= "EUR")

## Read Outcome data
 # Full path: /scratch/cfc85413/PUFAS/PSYCH/ADHD_30478444.a1effect.munge.rmInDels.uniq.tsv.gz
outcome_dat <- read_outcome_data(snps = clump_dat$SNP,filename = '/scratch/cfc85413/PUFAS/PSYCH/ADHD_30478444.a1effect.munge.rmInDels.uniq.tsv.gz',  sep="\t", snp_col= "SNP",  beta_col ="BETA", se_col = "SE", effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P", chr_col="CHR",pos_col= "BP",eaf_col = "FRQ")

## Harmonise the data 
harmonise_res<-harmonise_data(clump_dat, outcome_dat)

## MR = KEEP function 
  # to exclude all false results for Mr KEEP 
harmonise_res<-harmonise_res[harmonise_res$mr_keep==TRUE,]

## Senstivity analysis 
res_heterogenity<- mr_heterogeneity(harmonise_res)
res_ple<-mr_pleiotropy_test(harmonise_res)
res_leave<-mr_leaveoneout(harmonise_res, parameters = default_parameters(), method = mr_ivw)

## Change names
harmonise_res$id.exposure <- "DHA"
harmonise_res$exposure <- "DHA"
harmonise_res$outcome <- "ADHD"
harmonise_res$id.outcome <- "ADHD"

## Perform the MR function
mr_res<- mr(harmonise_res,parameters = default_parameters(), method_list = subset(mr_method_list(), use_by_default)$obj)

## IVW function 
mr_ivw_mre_claire <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))

  ivw.res <- summary(lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
  b <- ivw.res$coef["b_exp","Estimate"]
  # se <- ivw.res$coef["b_exp","Std. Error"]
  se <- ivw.res$coef["b_exp","Std. Error"]/min(ivw.res$sigma,1)
  pval <- 2 * pnorm(abs(b/se), lower.tail=FALSE)
  Q_df <- length(b_exp) - 1
  Q <- ivw.res$sigma^2 * Q_df
  Q_pval <- pchisq(Q, Q_df, lower.tail=FALSE)
  # from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
  # Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
  return(list(b = b, se = se, pval = pval, nsnp=length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}


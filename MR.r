#Load the libraries≈
library(TwoSampleMR)
library(data.table)
library(dplyr)

#read the exposure data
#filename can be the variant that has already been set
Args <- commandArgs()
exposure_data <-paste('/scratch/cfc85413/PUFAS/PSYCH/', Args[6], sep ="")
output_data <- paste('/scratch/cfc85413/PUFAS/PSYCH/', Args[7], sep = "")
exposure_name <- Args[8]
outcome_name <- Args[9]
output_file <-paste('/scratch/cfc85413/PUFAS/data/table/reverse_output/', Args[10], sep ="")
output_table <-paste('/scratch/cfc85413/PUFAS/data/table/reverse_table/', Args[11], sep ="")


exposure_dat <- read_exposure_data(filename=exposure_data,clump=FALSE,sep="\t",snp_col="SNP",beta_col="BETA",se_col="SE",effect_allele_col="A1",other_allele_col="A2", eaf_col="FRQ",pval_col="P", $

#Filter out the p value by lower than 5e-8 (GWAS significant threshold)
#sigificant_exposure <- filter(exposure_dat,pval.exposure <= 5e-8)
sigificant_exposure <- exposure_dat[exposure_dat$pval.exposure<=5e-8,]


#Clump data: read the exposure data
clump_dat <- clump_data(sigificant_exposure,  clump_kb = 10000, clump_r2 = 0.001,  clump_p1 = 5e-8,  clump_p2 = 5e-8, pop= "EUR")


#Using clumped data in the read_outcome_data function, this file does not have allele frequencey, still proceed forth
outcome_dat <- read_outcome_data(snps = clump_dat$SNP,filename =output_data,  sep="\t", snp_col= "SNP",  beta_col ="BETA", se_col = "SE",
effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P", chr_col="CHR",pos_col= "BP")

#harmonize data
harmonise_res<-harmonise_data(clump_dat, outcome_dat)

#Yitang's removal of genetic instruments, this filters out the new column of 'MR_Keep', you want to keep the TRUE data
harmonise_res<-harmonise_res[harmonise_res$mr_keep==TRUE,]

#Sensitively analysis
#res_heterogenity<- mr_heterogeneity(harmonise_res, parameters = default_parameters(), method_list = subset(mr_method_list(), heterogeneity_test & use_by_default)$obj)
res_heterogenity<- mr_heterogeneity(harmonise_res)

res_ple<-mr_pleiotropy_test(harmonise_res)
res_leave<-mr_leaveoneout(harmonise_res, parameters = default_parameters(), method = mr_ivw)

#Rename the exposure ID and the outcome ID
harmonise_res$id.exposure <- exposure_name
harmonise_res$id.outcome <- outcome_name

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

IVW_new_result <- mr_ivw_mre_claire(b_exp=harmonise_res$beta.exposure,b_out=harmonise_res$beta.outcome, se_exp=harmonise_res$se.exposure, se_out=harmonise_res$se.outcome)

#Do Mr
mr_res<- mr(harmonise_res,parameters = default_parameters(), method_list = subset(mr_method_list(), use_by_default)$obj)
write.table(mr, file = output_file, col.names=FALSE,row.names =F, quote =FALSE, sep ='\t')

#making the one row table to add to the header
b_MRE =mr_res[mr_res$method=='MR Egger',]$b
se_MRE =mr_res[mr_res$method=='MR Egger',]$se
pval_MRE =mr_res[mr_res$method=='MR Egger',]$pval
b_WMed =mr_res[mr_res$method=='Weighted median',]$b
se_WMed =mr_res[mr_res$method=='Weighted median',]$se
pval_WMed =mr_res[mr_res$method=='Weighted median',]$pval
b_IVW =IVW_new_result$b
se_IVW =IVW_new_result$se
pval_IVW =IVW_new_result$pval
b_SM =mr_res[mr_res$method=='Simple mode',]$b
se_SM =mr_res[mr_res$method=='Simple mode',]$se
pval_SM =mr_res[mr_res$method=='Simple mode',]$pval
b_WMODE =mr_res[mr_res$method=='Weighted mode',]$b
se_WMODE =mr_res[mr_res$method=='Weighted mode',]$se
pval_WMODE= mr_res[mr_res$method=='Weighted mode',]$pval
Q_val_het_Egger=res_heterogenity[res_heterogenity$method=='MR Egger',]$Q
Q_val_het_IVW=IVW_new_result$Q
Q_df_het_Egger=res_heterogenity[res_heterogenity$method=='MR Egger',]$Q_df
Q_df_het_IVW=IVW_new_result$Q_df
Q_pval_het_Egger=res_heterogenity[res_heterogenity$method=='MR Egger',]$Q_pval
Q_pval_het_IVW=IVW_new_result$Q_pval
Egger_intercept= res_ple$egger_intercept
pval_intercept= res_ple$pval

new_row <- c(exposure_name,outcome_name,b_MRE,se_MRE,pval_MRE,b_WMed,se_WMed,pval_WMed,b_IVW,se_IVW,pval_IVW,b_SM,se_SM,pval_SM,b_WMODE,se_WMODE,pval_WMODE,Q_val_het_Egger, Q_val_het_IVW,
 Q_df_het_Egger, Q_df_het_IVW,Q_pval_het_Egger, Q_pval_het_IVW, Egger_intercept,pval_intercept)

header<-c("exposure_name","outcome_name","b_MRE","se_MRE","pval_MRE","b_WMed","se_WMed","pval_WMed","b_IVW","se_IVW","pval_IVW","b_SM","se_SM","pval_SM","b_WMODE",
"se_WMODE","pval_WMODE", "Q_val_het_Egger", "Q_val_het_IVW", "Q_df_het_Egger","Q_df_het_IVW", "Q_pval_het_Egger", "Q_pval_het_IVW", "Egger_intercept", "pval_intercept")

total<-as.data.frame(rbind(header,new_row))

write.table(total[2,], file=output_table, col.names=FALSE,row.names =F, quote =FALSE, sep ='\t')


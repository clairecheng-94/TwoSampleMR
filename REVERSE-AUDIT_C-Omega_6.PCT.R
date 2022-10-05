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

#read the exposure data
exposure_dat<-read_exposure_data(filename ='/scratch/cfc85413/PUFAS/AUDIT_C_30336701.a1effect.munge.rmInDels.uniq.tsv.gz', clump=FALSE, sep = "\t" , snp_col = "SNP", beta_col ="BETA",  se_col = "SE",effect_allele_col = "A1",  other_allele_col = "A2",pval_col = "P", chr_col="CHR", pos_col="BP")

#Filter out the p value by lower than 5e-8 (GWAS significant threshold) 
sigificant_exposure <- filter(exposure_dat,pval.exposure <= 5e-8) 


#Clump data: read the exposure data 
clump_dat <- clump_data(sigificant_exposure,  clump_kb = 10000, clump_r2 = 0.001,  clump_p1 = 5e-8,  clump_p2 = 5e-8, pop= "EUR")

#Using clumped data in the read_outcome_data function, this file does not have allele frequencey, still proceed forth 
outcome_dat <- read_outcome_data(snps = clump_dat$SNP,filename = '/scratch/cfc85413/PUFAS/UKB_Omega_6_pct.a1effect.munge.rmInDels.uniq.tsv.gz', sep="\t", snp_col= "SNP",  beta_col ="BETA", se_col = "SE", effect_allele_col = "A1",other_allele_col = "A2",eaf_col = "FRQ",pval_col = "P", chr_col="CHR",pos_col= "BP")

#harmonize data 
res<-harmonise_data(clump_dat, outcome_dat) 


#Yitang's removal of genetic instruments, this filters out the new column of 'MR_Keep', you want to keep the TRUE data 
##not correct res_true<-filter(res, (mr_keep.exposure + mr_keep + mr_keep.outcome) > 0)
res<-res[res$mr_keep==TRUE,]

#Rename the exposure ID and the outcome ID to the Omega 6 and AUDI_C within rows 
res$id.exposure <-"AUDIT_C" 
res$id.outcome <- "Omega-6.pct"

res$exposure <- "AUDIT_C"
res$outcome <- "Omega-6.pct"

#Sensitively analysis 
res_heterogenity<- mr_heterogeneity(res, parameters = default_parameters(), method_list = subset(mr_method_list(), heterogeneity_test & use_by_default)$obj) 
res_ple<-mr_pleiotropy_test(res)
res_leave<-mr_leaveoneout(res, parameters = default_parameters(), method = mr_ivw)

#plots 
##funnel plot 
res_singlesnap<-mr_singlesnp(res, parameters = default_parameters(), single_method = "mr_wald_ratio", all_method = c("mr_ivw", "mr_egger_regression"))

#funnel plot
pdf("REVERSE_AUDIT_C-OMEGA6_PCT.funnelplot.pdf")
mr_funnel_plot(res_singlesnap)
dev.off()

##forest plot 
pdf("REVERSE_AUDIT_C-OMEGA6_PCT.forestplot.pdf")
forest_plot<-mr_forest_plot(res_singlesnap, parameters = default_parameters(),single_method = "mr_wald_ratio",all_method = c("mr_ivw", "mr_egger_regression"))
dev.off()

##leave one out 
pdf("REVERSE_AUDIT_C-OMEGA6_PCT.leaveoneoutplot.pdf")
res_leaveone<-mr_leaveoneout(res,parameters = default_parameters(), method = mr_ivw)
res_leaveone_plot<-mr_leaveoneout_plot(res_leaveone)
dev.off()

#Do Mr
mr_res<- mr(res,parameters = default_parameters(), method_list = subset(mr_method_list(), use_by_default)$obj) 

#scatter plot 
pdf("REVERSE_AUDIT_C-OMEGA6_PCT.SCATTERPLOT.pdf")
z <- exposure_dat[ ,("beta.exposure")]
y <- outcome_dat[ ,("beta.outcome")]
C <-plot(z, y, main = "AUDIT_C vs Omega6_PCT", xlab = "Omega6 beta values", ylab = "AUDIT_C beta values", pch=19, frame= FALSE)
abline(lm (z~y, data= C), col="blue")
pdf("AUDIT_C-OMEGA6_PCT.SCATTERPLOT.pdf")
dev.off()





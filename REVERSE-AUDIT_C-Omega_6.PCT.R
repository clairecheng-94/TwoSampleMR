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
harmonized_res<-harmonise_data(clump_dat, outcome_dat) 


#Yitang's removal of genetic instruments, this filters out the new column of 'MR_Keep', you want to keep the TRUE data 
##not correct res_true<-filter(res, (mr_keep.exposure + mr_keep + mr_keep.outcome) > 0)
harmonized_res<-harmonized_res[harmonized_res$mr_keep==TRUE,]

#Rename the exposure ID and the outcome ID to the Omega 6 and AUDI_C within rows 
harmonized_res$id.exposure <-"AUDIT_C" 
harmonized_res$id.outcome <- "Omega-6.pct"

harmonized_res$exposure <- "AUDIT_C"
harmonized_res$outcome <- "Omega-6.pct"

#Sensitively analysis 
res_heterogenity<- mr_heterogeneity(harmonized_res, parameters = default_parameters(), method_list = subset(mr_method_list(), heterogeneity_test & use_by_default)$obj) 
res_ple<-mr_pleiotropy_test(harmonized_res)
res_leave<-mr_leaveoneout(harmonized_res, parameters = default_parameters(), method = mr_ivw)

#plots 
##funnel plot 
res_singlesnap<-mr_singlesnp(harmonized_res, parameters = default_parameters(), single_method = "mr_wald_ratio", all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

#funnel plot
pdf("REVERSE_AUDIT_C-OMEGA6_PCT.funnelplot.pdf")
mr_funnel_plot(res_singlesnap)
dev.off()

##forest plot 
pdf("REVERSE_AUDIT_C-OMEGA6_PCT.forestplot.pdf")
mr_forest_plot(res_singlesnap, exponentiate = FALSE)
dev.off()

##leave one out 
pdf("REVERSE_AUDIT_C-OMEGA6_PCT.leaveoneoutplot.pdf")
mr_leaveoneout_plot(res_leave)
dev.off()


#Do Mr
mr_res<- mr(harmonized_res,parameters = default_parameters(), method_list = subset(mr_method_list(), use_by_default)$obj) 

#scatter plot 
pdf("REVERSE_AUDIT_C-OMEGA6_PCT.SCATTERPLOT.pdf")
mr_scatter_plot(mr_res, harmonized_res)
dev.off()

#funnel plot
pdf("REVERSE_AUDIT_C-OMEGA6_PCT.funnelplot.pdf")
mr_funnel_plot(res_singlesnap)
dev.off()

##forest plot 
pdf("REVERSE_AUDIT_C-OMEGA6_PCT.forestplot.pdf")
forest_plot<-mr_forest_plot(res_singlesnap,  exponentiate = FALSE)
dev.off()

##leave one out 
pdf("REVERSE_AUDIT_C-OMEGA6_PCT.leaveoneoutplot.pdf")
res_leaveone_plot<-mr_leaveoneout_plot(res_leave)
dev.off()

#Do Mr
mr_res<- mr(harmonized_res,parameters = default_parameters(), method_list = subset(mr_method_list(), use_by_default)$obj) 







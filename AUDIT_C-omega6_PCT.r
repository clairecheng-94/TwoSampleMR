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
#filename can be the variant that has already been set 
exposure_dat<-read_exposure_data(filename ='/scratch/cfc85413/PUFAS/UKB_Omega_6_pct.a1effect.munge.rmInDels.uniq.tsv.gz', clump=FALSE, sep = "\t" , snp_col = "SNP", beta_col ="BETA",  se_col = "SE",effect_allele_col = "A1",  other_allele_col = "A2", eaf_col = "FRQ",pval_col = "P", chr_col="CHR", pos_col="BP")

#Filter out the p value by lower than 5e-8 (GWAS significant threshold) 
sigificant_exposure <- filter(exposure_dat,pval.exposure <= 5e-8) 


#Clump data: read the exposure data 
clump_dat <- clump_data(sigificant_exposure,  clump_kb = 10000, clump_r2 = 0.001,  clump_p1 = 5e-8,  clump_p2 = 5e-8, pop= "EUR")
  

#Using clumped data in the read_outcome_data function, this file does not have allele frequencey, still proceed forth 
outcome_dat <- read_outcome_data(snps = clump_dat$SNP,filename = '/scratch/cfc85413/PUFAS/AUDIT_C_30336701.a1effect.munge.rmInDels.uniq.tsv.gz',  sep="\t", snp_col= "SNP",  beta_col ="BETA", se_col = "SE", effect_allele_col = "A1",other_allele_col = "A2",pval_col = "P", chr_col="CHR",pos_col= "BP")


#harmonize data 
harmonise_res<-harmonise_data(clump_dat, outcome_dat) 

#Sensitively analysis 
res_heterogenity<- mr_heterogeneity(harmonise_res, parameters = default_parameters(), method_list = subset(mr_method_list(), heterogeneity_test & use_by_default)$obj) 
write.table(res_heterogenity, file="AUDIT_C_ON_OMEGA-6_PCT.heterogenity.table", sep="\t", quote=FALSE, col.names=T, row.names=F)
res_ple<-mr_pleiotropy_test(harmonise_res)
res_leave<-mr_leaveoneout(harmonise_res, parameters = default_parameters(), method = mr_ivw)


#Yitang's removal of genetic instruments, this filters out the new column of 'MR_Keep', you want to keep the TRUE data 
##not correct res_true<-filter(res, (mr_keep.exposure + mr_keep + mr_keep.outcome) > 0)
harmonise_res<-harmonise_res[harmonise_res$mr_keep==TRUE,]
             
#Rename the exposure ID and the outcome ID to the Omega 6 and AUDI_C within rows 
harmonise_res$id.exposure <- "Omega-6.pct"
harmonise_res$id.outcome <- "AUDIT_C"

#change names for exposure and outcome
harmonise_res$exposure <- "Omega-6.pct"
harmonise_res$outcome <- "AUDIT_C"
write.table(harmonise_res, file = "AUDIT_C_ON_OMEGA-6_PCT.harmonized.table", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
#plots 

res_singlesnap<-mr_singlesnp(harmonise_res, parameters = default_parameters(), single_method = "mr_wald_ratio", all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
write.table(mr_res, file = "AUDIT_C_ON_OMEGA-6_PCT.MR.table", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

#funnel plot
pdf("AUDIT_C-OMEGA6_PCT.funnelplot.pdf")
mr_funnel_plot(res_singlesnap)
dev.off()

##forest plot 
pdf("AUDIT_C-OMEGA6_PCT.forestplot.pdf")
mr_forest_plot(res_singlesnap, exponentiate = FALSE)
dev.off()

##leave one out 
pdf("AUDIT_C-OMEGA6_PCT.leaveoneoutplot.pdf")
mr_leaveoneout_plot(res_leave)
dev.off()


#Do Mr
mr_res<- mr(harmonise_res,parameters = default_parameters(), method_list = subset(mr_method_list(), use_by_default)$obj) 

#scatter plot 
pdf("AUDIT_C-OMEGA6_PCT.SCATTERPLOT.pdf")
mr_scatter_plot(mr_res, harmonise_res)
dev.off()

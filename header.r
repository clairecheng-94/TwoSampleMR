exposure_name="Omega-3"
outcome_name="AUDIT_P"

b_Egger=mr_res[mr_res$method=='MR Egger',]$b
se_Egger=mr_res[mr_res$method=='MR Egger',]$se
pval_Egger=mr_res[mr_res$method=='MR Egger',]$pval
b_WMed=mr_res[mr_res$method=='Weighted median',]$b
se_WMed=mr_res[mr_res$method=='Weighted median',]$se
pval_WMed=mr_res[mr_res$method=='Weighted median',]$pval
b_IVW=mr_res[mr_res$method=='Inverse variance weighted',]$b
se_IVW=mr_res[mr_res$method=='Inverse variance weighted',]$se
pval_IVW=mr_res[mr_res$method=='Inverse variance weighted',]$pval
b_SM=mr_res[mr_res$method=='Simple mode',]$b
se_SM=mr_res[mr_res$method=='Simple mode',]$se
pval_SM=mr_res[mr_res$method=='Simple mode',]$pval
b_WMODE=mr_res[mr_res$method=='Weighted mode',]$b
se_WMODE=mr_res[mr_res$method=='Weighted mode',]$se
pval_WMODE=mr_res[mr_res$method=='Weighted mode',]$pval

new_row <- c(exposure_name,outcome_name,b_Egger,se_Egger,pval_Egger,b_WMed,se_WMed,pval_WMed,b_IVW,se_IVW,pval_IVW,b_SM,se_SM,pval_SM,b_WMODE,se_WMODE,pval_WMODE)
header<-c("exposure_name","outcome_name","b_Egger","se_Egger","pval_Egger","b_WMed","se_WMed","pval_WMed","b_IVW","se_IVW","pval_IVW","b_SM","se_SM","pval_SM","b_WMODE","se_WMODE","pval_WMODE")

 total<-as.data.frame(rbind(header,new_row))

#change file name each time 
 write.table(total[2,], file='/scratch/cfc85413/PUFAS/data/table/Omega-3_', col.names=FALSE,row.names =F, quote =FALSE, sep ='\t')

#for loop to combine files
for file in `ls *.txt`;do cat $file >> total.csv;done

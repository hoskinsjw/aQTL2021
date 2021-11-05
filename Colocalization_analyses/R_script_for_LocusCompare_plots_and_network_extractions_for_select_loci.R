### This script is for extra analyses on the most interesting select loci (1p36.1, 7q32 and 12p13.1)

#install.packages("devtools")
#library(devtools)
#install_github("jrs95/hyprcoloc", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = F)
#browseVignettes("hyprcoloc") The install kept failing when trying to build the vignettes, so I disabled that.
#devtools::install_github("boxiangliu/locuscomparer")
library(hyprcoloc)
library(locuscomparer)

setwd("YOUR WORKING DIRECTORY")

# Read in the LD matrices, and the GWAS and QTL data. The file locations are relative to your working directory, so adjust accordingly.
bmi=read.table("./Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt",sep = "\t",header = T)
t2d=read.table("./Mahajan.NatGenet2018b.T2Dbmiadj.European.with.rsIDs.txt",sep = "\t",header = T)
hdl=read.table("./jointGwasMc_HDL.txt",sep = "\t",header = T)
triG=read.table("./jointGwasMc_TG.txt",sep = "\t",header = T)
cis_eqtl=read.table("./Eurobats_adipose_select_loci_cis-eQTLs_from_INT_logTPM.txt",sep = "\t",header = T)
cis_aqtl=read.table("./Eurobats_adipose_select_loci_cis-aQTLs_from_unnormalized_activities.txt",sep = "\t",header = T)
trans_bmi_eqtl=read.table("./Eurobats_adipose_select_loci_trans-eQTLs_for_BMI_MRs.txt",sep = "\t",header = T)
trans_bmi_aqtl=read.table("./Eurobats_adipose_select_loci_trans-aQTLs_for_BMI_MRs.txt",sep = "\t",header = T)
trans_t2d_eqtl=read.table("./Eurobats_adipose_select_loci_trans-eQTLs_for_HOMA-IR_MRs.txt",sep = "\t",header = T)
trans_t2d_aqtl=read.table("./Eurobats_adipose_select_loci_trans-aQTLs_for_HOMA-IR_MRs.txt",sep = "\t",header = T)
trans_hdl_eqtl=read.table("./Eurobats_adipose_select_loci_trans-eQTLs_for_HDL_MRs.txt",sep = "\t",header = T)
trans_hdl_aqtl=read.table("./Eurobats_adipose_select_loci_trans-aQTLs_for_HDL_MRs.txt",sep = "\t",header = T)
trans_triG_eqtl=read.table("./Eurobats_adipose_select_loci_trans-eQTLs_for_TriG_MRs.txt",sep = "\t",header = T)
trans_triG_aqtl=read.table("./Eurobats_adipose_select_loci_trans-aQTLs_for_TriG_MRs.txt",sep = "\t",header = T)
ld_files=c("Eurobats_chr1p36.1_LD_matrix.txt","Eurobats_chr7q32_LD_matrix.txt","Eurobats_chr12p13.33_LD_matrix.txt","Eurobats_chr12p13.1_LD_matrix.txt")

ld=list()
index=1
for(i in ld_files){
  ld[[index]]=read.table(paste("./",i,sep=""),sep = "\t",header = F)
  rownames(ld[[index]])=ld[[index]][,3]
  ld[[index]]=ld[[index]][,-c(1:5)]
  colnames(ld[[index]])=rownames(ld[[index]])
  index=index+1
}

# The HDL GWAS data has coordinates for hg18 and hg 19, but I need to have CHR and POS columns (based on hg19) instead.
colnames(hdl)=c("CHR","POS","SNP","A1","A2","BETA","SE","N","P","Freq.A1.1000G.EUR")
hdl$CHR=gsub("chr","",hdl$CHR)
hdl$CHR=as.numeric(gsub(":.*","",hdl$CHR))
# This introduced NAs, but only for 3 SNPs without rsIDs (labeled only as ".")
hdl=hdl[!is.na(hdl$CHR),]
hdl$POS=as.numeric(gsub("chr.*:","",hdl$POS))

# The TriG GWAS data has coordinates for hg18 and hg 19, but I need to have CHR and POS columns (based on hg19) instead.
colnames(triG)=c("CHR","POS","SNP","A1","A2","BETA","SE","N","P","Freq.A1.1000G.EUR")
triG$CHR=gsub("chr","",triG$CHR)
triG$CHR=as.numeric(gsub(":.*","",triG$CHR))
# This introduced NAs, but only for 3 SNPs without rsIDs (labeled only as ".")
triG=triG[!is.na(triG$CHR),]
triG$POS=as.numeric(gsub("chr.*:","",triG$POS))

# Filter GWAS, QTL and LD data to the same SNPs
filt_bmi=bmi[na.omit(match(c(rownames(ld[[1]]),rownames(ld[[2]]),rownames(ld[[3]]),rownames(ld[[4]])),bmi$SNP)),]
filt_bmi=filt_bmi[na.omit(match(unique(cis_eqtl$snps),filt_bmi$SNP)),]
filt_bmi=filt_bmi[na.omit(match(t2d$rsID,filt_bmi$SNP)),]
filt_bmi=filt_bmi[na.omit(match(hdl$SNP,filt_bmi$SNP)),]
filt_bmi=filt_bmi[na.omit(match(triG$SNP,filt_bmi$SNP)),]
filt_t2d=t2d[na.omit(match(filt_bmi$SNP,t2d$rsID)),]
filt_hdl=hdl[na.omit(match(filt_bmi$SNP,hdl$SNP)),]
filt_triG=triG[na.omit(match(filt_bmi$SNP,triG$SNP)),]
filt_cis_eqtl=cis_eqtl[cis_eqtl$snps %in% filt_bmi$SNP,]
filt_cis_aqtl=cis_aqtl[cis_aqtl$snps %in% filt_bmi$SNP,]
filt_trans_bmi_eqtl=trans_bmi_eqtl[trans_bmi_eqtl$snps %in% filt_bmi$SNP,]
filt_trans_bmi_aqtl=trans_bmi_aqtl[trans_bmi_aqtl$snps %in% filt_bmi$SNP,]
filt_trans_t2d_eqtl=trans_t2d_eqtl[trans_t2d_eqtl$snps %in% filt_t2d$rsID,]
filt_trans_t2d_aqtl=trans_t2d_aqtl[trans_t2d_aqtl$snps %in% filt_t2d$rsID,]
filt_trans_hdl_eqtl=trans_hdl_eqtl[trans_hdl_eqtl$snps %in% filt_hdl$SNP,]
filt_trans_hdl_aqtl=trans_hdl_aqtl[trans_hdl_aqtl$snps %in% filt_hdl$SNP,]
filt_trans_triG_eqtl=trans_triG_eqtl[trans_triG_eqtl$snps %in% filt_triG$SNP,]
filt_trans_triG_aqtl=trans_triG_aqtl[trans_triG_aqtl$snps %in% filt_triG$SNP,]
filt_ld=list()
filt_ld[[1]]=ld[[1]][filt_bmi$SNP[filt_bmi$CHR==1],filt_bmi$SNP[filt_bmi$CHR==1]]
filt_ld[[2]]=ld[[2]][filt_bmi$SNP[filt_bmi$CHR==7],filt_bmi$SNP[filt_bmi$CHR==7]]
filt_ld[[3]]=ld[[3]][filt_bmi$SNP[filt_bmi$CHR==12 & filt_bmi$POS<1400000],filt_bmi$SNP[filt_bmi$CHR==12 & filt_bmi$POS<1400000]]
filt_ld[[4]]=ld[[4]][filt_bmi$SNP[filt_bmi$CHR==12 & filt_bmi$POS>1400000],filt_bmi$SNP[filt_bmi$CHR==12 & filt_bmi$POS>1400000]]

# Let's free up some memory by dropping the huge trans-QTL data.frames
rm(trans_bmi_eqtl)
rm(trans_bmi_aqtl)
rm(trans_t2d_eqtl)
rm(trans_t2d_aqtl)
rm(trans_hdl_eqtl)
rm(trans_hdl_aqtl)
rm(trans_triG_eqtl)
rm(trans_triG_aqtl)

# Add chromosome and position to the QTLs for sorting
filt_cis_eqtl$chr=filt_bmi[match(filt_cis_eqtl$snps,filt_bmi$SNP),1]
filt_cis_aqtl$chr=filt_bmi[match(filt_cis_aqtl$snps,filt_bmi$SNP),1]
filt_trans_bmi_eqtl$chr=filt_bmi[match(filt_trans_bmi_eqtl$snps,filt_bmi$SNP),1]
filt_trans_bmi_aqtl$chr=filt_bmi[match(filt_trans_bmi_aqtl$snps,filt_bmi$SNP),1]
filt_trans_t2d_eqtl$chr=filt_bmi[match(filt_trans_t2d_eqtl$snps,filt_bmi$SNP),1]
filt_trans_t2d_aqtl$chr=filt_bmi[match(filt_trans_t2d_aqtl$snps,filt_bmi$SNP),1]
filt_trans_hdl_eqtl$chr=filt_bmi[match(filt_trans_hdl_eqtl$snps,filt_bmi$SNP),1]
filt_trans_hdl_aqtl$chr=filt_bmi[match(filt_trans_hdl_aqtl$snps,filt_bmi$SNP),1]
filt_trans_triG_eqtl$chr=filt_bmi[match(filt_trans_triG_eqtl$snps,filt_bmi$SNP),1]
filt_trans_triG_aqtl$chr=filt_bmi[match(filt_trans_triG_aqtl$snps,filt_bmi$SNP),1]
filt_cis_eqtl$position=filt_bmi[match(filt_cis_eqtl$snps,filt_bmi$SNP),2]
filt_cis_aqtl$position=filt_bmi[match(filt_cis_aqtl$snps,filt_bmi$SNP),2]
filt_trans_bmi_eqtl$position=filt_bmi[match(filt_trans_bmi_eqtl$snps,filt_bmi$SNP),2]
filt_trans_bmi_aqtl$position=filt_bmi[match(filt_trans_bmi_aqtl$snps,filt_bmi$SNP),2]
filt_trans_t2d_eqtl$position=filt_bmi[match(filt_trans_t2d_eqtl$snps,filt_bmi$SNP),2]
filt_trans_t2d_aqtl$position=filt_bmi[match(filt_trans_t2d_aqtl$snps,filt_bmi$SNP),2]
filt_trans_hdl_eqtl$position=filt_bmi[match(filt_trans_hdl_eqtl$snps,filt_bmi$SNP),2]
filt_trans_hdl_aqtl$position=filt_bmi[match(filt_trans_hdl_aqtl$snps,filt_bmi$SNP),2]
filt_trans_triG_eqtl$position=filt_bmi[match(filt_trans_triG_eqtl$snps,filt_bmi$SNP),2]
filt_trans_triG_aqtl$position=filt_bmi[match(filt_trans_triG_aqtl$snps,filt_bmi$SNP),2]

# Sort by chr and position
filt_bmi=filt_bmi[order(filt_bmi$CHR,filt_bmi$POS),]
filt_t2d=filt_t2d[order(filt_t2d$Chr,filt_t2d$Pos),]
filt_hdl=filt_hdl[order(filt_hdl$CHR,filt_hdl$POS),]
filt_triG=filt_triG[order(filt_triG$CHR,filt_triG$POS),]
filt_cis_eqtl=filt_cis_eqtl[order(filt_cis_eqtl$chr,filt_cis_eqtl$position),]
filt_cis_aqtl=filt_cis_aqtl[order(filt_cis_aqtl$chr,filt_cis_aqtl$position),]
filt_trans_bmi_eqtl=filt_trans_bmi_eqtl[order(filt_trans_bmi_eqtl$chr,filt_trans_bmi_eqtl$position),]
filt_trans_bmi_aqtl=filt_trans_bmi_aqtl[order(filt_trans_bmi_aqtl$chr,filt_trans_bmi_aqtl$position),]
filt_trans_t2d_eqtl=filt_trans_t2d_eqtl[order(filt_trans_t2d_eqtl$chr,filt_trans_t2d_eqtl$position),]
filt_trans_t2d_aqtl=filt_trans_t2d_aqtl[order(filt_trans_t2d_aqtl$chr,filt_trans_t2d_aqtl$position),]
filt_trans_hdl_eqtl=filt_trans_hdl_eqtl[order(filt_trans_hdl_eqtl$chr,filt_trans_hdl_eqtl$position),]
filt_trans_hdl_aqtl=filt_trans_hdl_aqtl[order(filt_trans_hdl_aqtl$chr,filt_trans_hdl_aqtl$position),]
filt_trans_triG_eqtl=filt_trans_triG_eqtl[order(filt_trans_triG_eqtl$chr,filt_trans_triG_eqtl$position),]
filt_trans_triG_eqtl=filt_trans_triG_eqtl[order(filt_trans_triG_eqtl$chr,filt_trans_triG_eqtl$position),]

### LocusCompare plots

## 1p36.1
# First, grab the necessary P-values for the SNPs used in the HyPrColoc analyses for the traits of interest
bmi1=filt_bmi[filt_bmi$CHR==1,c(3,9)]
loc1_eEPHB2=filt_cis_eqtl[filt_cis_eqtl$gene=="EPHB2",c(1,4)]
loc1_aEPHB2=filt_cis_aqtl[filt_cis_aqtl$gene=="EPHB2",c(1,4)]
loc1_eZNF436=filt_cis_eqtl[filt_cis_eqtl$gene=="ZNF436",c(1,4)]
loc1_aZNF436=filt_cis_aqtl[filt_cis_aqtl$gene=="ZNF436",c(1,4)]
loc1_eTCEA3=filt_cis_eqtl[filt_cis_eqtl$gene=="TCEA3",c(1,4)]
loc1_aTCEA3=filt_cis_aqtl[filt_cis_aqtl$gene=="TCEA3",c(1,4)]
loc1_eLASP1=filt_trans_bmi_eqtl[filt_trans_bmi_eqtl$gene=="LASP1" & filt_trans_bmi_eqtl$chr==1,c(1,4)]
loc1_aLASP1=filt_trans_bmi_aqtl[filt_trans_bmi_aqtl$gene=="LASP1" & filt_trans_bmi_aqtl$chr==1,c(1,4)]
loc1_eRASSF4=filt_trans_bmi_eqtl[filt_trans_bmi_eqtl$gene=="RASSF4" & filt_trans_bmi_eqtl$chr==1,c(1,4)]
loc1_aRASSF4=filt_trans_bmi_aqtl[filt_trans_bmi_aqtl$gene=="RASSF4" & filt_trans_bmi_aqtl$chr==1,c(1,4)]
loc1_aGNA14=filt_trans_bmi_aqtl[filt_trans_bmi_aqtl$gene=="GNA14" & filt_trans_bmi_aqtl$chr==1,c(1,4)]
loc1_aDOK5=filt_trans_bmi_aqtl[filt_trans_bmi_aqtl$gene=="DOK5" & filt_trans_bmi_aqtl$chr==1,c(1,4)]
colnames(bmi1)=c("rsid","pval")
colnames(loc1_eEPHB2)=c("rsid","pval")
colnames(loc1_aEPHB2)=c("rsid","pval")
colnames(loc1_eZNF436)=c("rsid","pval")
colnames(loc1_aZNF436)=c("rsid","pval")
colnames(loc1_eTCEA3)=c("rsid","pval")
colnames(loc1_aTCEA3)=c("rsid","pval")
colnames(loc1_eLASP1)=c("rsid","pval")
colnames(loc1_aLASP1)=c("rsid","pval")
colnames(loc1_eRASSF4)=c("rsid","pval")
colnames(loc1_aRASSF4)=c("rsid","pval")
colnames(loc1_aGNA14)=c("rsid","pval")
colnames(loc1_aDOK5)=c("rsid","pval")
rownames(bmi1)=bmi1$rsid
rownames(loc1_eEPHB2)=loc1_eEPHB2$rsid
rownames(loc1_aEPHB2)=loc1_aEPHB2$rsid
rownames(loc1_eZNF436)=loc1_eZNF436$rsid
rownames(loc1_aZNF436)=loc1_aZNF436$rsid
rownames(loc1_eTCEA3)=loc1_eTCEA3$rsid
rownames(loc1_aTCEA3)=loc1_aTCEA3$rsid
rownames(loc1_eLASP1)=loc1_eLASP1$rsid
rownames(loc1_aLASP1)=loc1_aLASP1$rsid
rownames(loc1_eRASSF4)=loc1_eRASSF4$rsid
rownames(loc1_aRASSF4)=loc1_aRASSF4$rsid
rownames(loc1_aGNA14)=loc1_aGNA14$rsid
rownames(loc1_aDOK5)=loc1_aDOK5$rsid

# Check out some relevant LocusCompare plots before picking the which to write to file
locuscompare(in_fn1=bmi1,in_fn2=loc1_eEPHB2,title1 = "BMI GWAS", title2 = "EPHB2 cis-eQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aEPHB2,title1 = "BMI GWAS", title2 = "EPHB2 cis-aQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aEPHB2,title1 = "BMI GWAS", title2 = "EPHB2 cis-aQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aEPHB2,title1 = "BMI GWAS", title2 = "EPHB2 cis-aQTL",snp = "rs12408468") # Potential 3rd BMI signal?
locuscompare(in_fn1=loc1_eEPHB2,in_fn2=loc1_aEPHB2,title1 = "EPHB2 cis-eQTL", title2 = "EPHB2 cis-aQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_eZNF436,title1 = "BMI GWAS", title2 = "ZNF436 cis-eQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aZNF436,title1 = "BMI GWAS", title2 = "ZNF436 cis-aQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=loc1_eZNF436,in_fn2=loc1_aZNF436,title1 = "ZNF436 cis-eQTL", title2 = "ZNF436 cis-aQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aZNF436,title1 = "BMI GWAS", title2 = "ZNF436 cis-aQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_eTCEA3,title1 = "BMI GWAS", title2 = "TCEA3 cis-eQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aTCEA3,title1 = "BMI GWAS", title2 = "TCEA3 cis-aQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_eLASP1,title1 = "BMI GWAS", title2 = "LASP1 trans-eQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aLASP1,title1 = "BMI GWAS", title2 = "LASP1 trans-aQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=loc1_aEPHB2,in_fn2=loc1_aLASP1,title1 = "EPHB2 cis-aQTL", title2 = "LASP1 trans-aQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=loc1_aEPHB2,in_fn2=loc1_aLASP1,title1 = "EPHB2 cis-aQTL", title2 = "LASP1 trans-aQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_eRASSF4,title1 = "BMI GWAS", title2 = "RASSF4 trans-eQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aRASSF4,title1 = "BMI GWAS", title2 = "RASSF4 trans-aQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=loc1_aEPHB2,in_fn2=loc1_aRASSF4,title1 = "EPHB2 cis-aQTL", title2 = "RASSF4 trans-aQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=loc1_aEPHB2,in_fn2=loc1_aRASSF4,title1 = "EPHB2 cis-aQTL", title2 = "RASSF4 trans-aQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aLASP1,title1 = "BMI GWAS", title2 = "LASP1 trans-aQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aRASSF4,title1 = "BMI GWAS", title2 = "RASSF4 trans-aQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aGNA14,title1 = "BMI GWAS", title2 = "GNA14 trans-aQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aGNA14,title1 = "BMI GWAS", title2 = "GNA14 trans-aQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=loc1_aEPHB2,in_fn2=loc1_aGNA14,title1 = "EPHB2 cis-aQTL", title2 = "GNA14 trans-aQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=loc1_aEPHB2,in_fn2=loc1_aGNA14,title1 = "EPHB2 cis-aQTL", title2 = "GNA14 trans-aQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aDOK5,title1 = "BMI GWAS", title2 = "DOK5 trans-aQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=bmi1,in_fn2=loc1_aDOK5,title1 = "BMI GWAS", title2 = "DOK5 trans-aQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=loc1_aEPHB2,in_fn2=loc1_aDOK5,title1 = "EPHB2 cis-aQTL", title2 = "DOK5 trans-aQTL",snp = "rs6692586") # Top GWAS SNP
locuscompare(in_fn1=loc1_aEPHB2,in_fn2=loc1_aDOK5,title1 = "EPHB2 cis-aQTL", title2 = "DOK5 trans-aQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
locuscompare(in_fn1=loc1_eEPHB2,in_fn2=loc1_aEPHB2,title1 = "EPHB2 cis-eQTL", title2 = "EPHB2 cis-aQTL",snp = "rs4654828") # Top multi-QTL cluster SNP
# This locus is complex and consequently difficult to interpret. The EPHB2 cis-aQTL and various BMI MR trans-aQTL signals suggest 2 functional signals
# represented by rs6692586 (the top BMI GWAS SNP) and rs4654828 (the top trans-aQTL signal for many BMI MRs). The EPHB2 cis-aQTL has these two SNPs
# at roughly equal strength while rs6692586 is clearly stronger for BMI and rs4654828 is clearly stronger for the trans-aQTLs. Perhaps the best
# hypothetical explanations for these observations is that rs6692586 operates in cis thru effects on EPHB2 expression and activity, while rs4654828
# has an alternative proximal effect that distally affects the activities of many correlated BMI MRs, including EPHB2, which shows up as a bump in
# in the EPHB2 aQTL signal. The proximal effect of rs4654828 might be on ZNF436 activity, but this probably cannot be mediated via expression levels
# since there is an extremely strong cis-eQTL for ZNF436 at this locus that does not overlap the BMI signal or the cis-aQTL signal. I looked into the
# the position of rs4654828, but it is quite far away from ZNF436 in a LACTBL1 intron. LACTBL1 is apparently not expressed in our adipose tissue, so
# it is hard to imagine how it could be mediating the effect on BMI within adipose. It is best expressed in testis, which does have a sig eQTL
# between rs4654828-LACTBL1, but this also doesn't seem relevant to BMI. So if rs4654828 does affect ZNF436 activity in adipose, it is mediated
# some other way. Interestingly, rs4654828 in GTEx does show sig eQTLs with ZNF436 in other tissue types (Skin, Aorta, Tibial Artery, Esophagus,
# Tibial Nerve and Thyroid). It's hard to imagine how ZNF436 expression effects in other tissues could be relevant to ZNF436 activity in adipose.
# There are also rs4654828-TCEA3 eQTLs in Skin and Skeletal Muscle, and a TCEA3 splicing QTL in skin. The TCEA3 eQTL/aQTL LocusCompare plots do
# not look like TCEA3 is relevant to either BMI signal in adipose. Regardless, this sort of scenario might manifest epistatic effects on BMI and 
# EPHB2 between these two SNPs. This is easy enough to test for EPHB2 activity, but I can't test it for BMI.

# Let's write some to PDFs
pdf("rs6692586-EPHB2_eQTL_and_BMI_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi1,in_fn2=loc1_eEPHB2,title1 = "BMI GWAS", title2 = "EPHB2 cis-eQTL",snp = "rs6692586")
dev.off()
pdf("rs6692586-EPHB2_aQTL_and_BMI_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi1,in_fn2=loc1_aEPHB2,title1 = "BMI GWAS", title2 = "EPHB2 cis-aQTL",snp = "rs6692586")
dev.off()
pdf("rs6692586-EPHB2_eQTL_and_aQTL_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=loc1_eEPHB2,in_fn2=loc1_aEPHB2,title1 = "EPHB2 cis-eQTL", title2 = "EPHB2 cis-aQTL",snp = "rs6692586")
dev.off()
pdf("rs4654828-EPHB2_aQTL_and_BMI_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi1,in_fn2=loc1_aEPHB2,title1 = "BMI GWAS", title2 = "EPHB2 cis-aQTL",snp = "rs4654828")
dev.off()
pdf("rs6692586-ZNF436_aQTL_and_BMI_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi1,in_fn2=loc1_aZNF436,title1 = "BMI GWAS", title2 = "ZNF436 cis-aQTL",snp = "rs6692586")
dev.off()
pdf("rs4654828-ZNF436_aQTL_and_BMI_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi1,in_fn2=loc1_aZNF436,title1 = "BMI GWAS", title2 = "ZNF436 cis-aQTL",snp = "rs4654828")
dev.off()
pdf("rs6692586-DOK5_aQTL_and_BMI_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi1,in_fn2=loc1_aDOK5,title1 = "BMI GWAS", title2 = "DOK5 cis-aQTL",snp = "rs6692586")
dev.off()
pdf("rs4654828-DOK5_aQTL_and_BMI_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi1,in_fn2=loc1_aDOK5,title1 = "BMI GWAS", title2 = "DOK5 cis-aQTL",snp = "rs4654828")
dev.off()
pdf("rs6692586-RASSF4_aQTL_and_BMI_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi1,in_fn2=loc1_aRASSF4,title1 = "BMI GWAS", title2 = "RASSF4 cis-aQTL",snp = "rs6692586")
dev.off()
pdf("rs4654828-RASSF4_aQTL_and_BMI_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi1,in_fn2=loc1_aRASSF4,title1 = "BMI GWAS", title2 = "RASSF4 cis-aQTL",snp = "rs4654828")
dev.off()
pdf("rs6692586-GNA14_aQTL_and_BMI_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi1,in_fn2=loc1_aGNA14,title1 = "BMI GWAS", title2 = "GNA14 cis-aQTL",snp = "rs6692586")
dev.off()
pdf("rs4654828-GNA14_aQTL_and_BMI_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi1,in_fn2=loc1_aGNA14,title1 = "BMI GWAS", title2 = "GNA14 cis-aQTL",snp = "rs4654828")
dev.off()
pdf("rs6692586-LASP1_aQTL_and_BMI_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi1,in_fn2=loc1_aLASP1,title1 = "BMI GWAS", title2 = "LASP1 cis-aQTL",snp = "rs6692586")
dev.off()
pdf("rs4654828-LASP1_aQTL_and_BMI_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi1,in_fn2=loc1_aLASP1,title1 = "BMI GWAS", title2 = "LASP1 cis-aQTL",snp = "rs4654828")
dev.off()
pdf("rs6692586-LASP1_aQTL_and_EPHB2_aQTL_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=loc1_aEPHB2,in_fn2=loc1_aLASP1,title1 = "EPHB2 cis-aQTL", title2 = "LASP1 cis-aQTL",snp = "rs6692586")
dev.off()
pdf("rs6692586-GNA14_aQTL_and_EPHB2_aQTL_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=loc1_aEPHB2,in_fn2=loc1_aGNA14,title1 = "EPHB2 cis-aQTL", title2 = "GNA14 cis-aQTL",snp = "rs6692586")
dev.off()
pdf("rs6692586-RASSF4_aQTL_and_EPHB2_aQTL_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=loc1_aEPHB2,in_fn2=loc1_aRASSF4,title1 = "EPHB2 cis-aQTL", title2 = "RASSF4 cis-aQTL",snp = "rs6692586")
dev.off()
pdf("rs6692586-DOK5_aQTL_and_EPHB2_aQTL_1p36_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=loc1_aEPHB2,in_fn2=loc1_aDOK5,title1 = "EPHB2 cis-aQTL", title2 = "DOK5 cis-aQTL",snp = "rs6692586")
dev.off()

## 7q32
# First, grab the necessary P-values for the SNPs used in the HyPrColoc analyses for the traits of interest
bmi2=filt_bmi[filt_bmi$CHR==7,c(3,9)]
t2d2=filt_t2d[filt_t2d$Chr==7,c(1,9)]
hdl2=filt_hdl[filt_hdl$CHR==7,c(3,9)]
triG2=filt_triG[filt_triG$CHR==7,c(3,9)]
loc2_eLINC=filt_cis_eqtl[filt_cis_eqtl$gene=="LINC-PINT",c(1,4)]
loc2_aLINC=filt_cis_aqtl[filt_cis_aqtl$gene=="LINC-PINT",c(1,4)]
loc2_eKLF14=filt_cis_eqtl[filt_cis_eqtl$gene=="KLF14",c(1,4)]
loc2_aKLF14=filt_cis_aqtl[filt_cis_aqtl$gene=="KLF14",c(1,4)]
loc2_eAC=filt_cis_eqtl[filt_cis_eqtl$gene=="AC016831.7",c(1,4)]
loc2_eTBX4=filt_trans_bmi_eqtl[filt_trans_bmi_eqtl$gene=="TBX4" & filt_trans_bmi_eqtl$chr==7,c(1,4)]
loc2_aTBX4=filt_trans_bmi_aqtl[filt_trans_bmi_aqtl$gene=="TBX4" & filt_trans_bmi_aqtl$chr==7,c(1,4)]
loc2_eGNB1=filt_trans_t2d_eqtl[filt_trans_t2d_eqtl$gene=="GNB1" & filt_trans_t2d_eqtl$chr==7,c(1,4)]
loc2_aGNB1=filt_trans_t2d_aqtl[filt_trans_t2d_aqtl$gene=="GNB1" & filt_trans_t2d_aqtl$chr==7,c(1,4)]
loc2_eESR2=filt_trans_hdl_eqtl[filt_trans_hdl_eqtl$gene=="ESR2" & filt_trans_hdl_eqtl$chr==7,c(1,4)]
loc2_aESR2=filt_trans_hdl_aqtl[filt_trans_hdl_aqtl$gene=="ESR2" & filt_trans_hdl_aqtl$chr==7,c(1,4)]
loc2_eNR2F1=filt_trans_hdl_eqtl[filt_trans_hdl_eqtl$gene=="NR2F1" & filt_trans_hdl_eqtl$chr==7,c(1,4)]
loc2_aNR2F1=filt_trans_hdl_aqtl[filt_trans_hdl_aqtl$gene=="NR2F1" & filt_trans_hdl_aqtl$chr==7,c(1,4)]
loc2_eAGT=filt_trans_triG_eqtl[filt_trans_triG_eqtl$gene=="AGT" & filt_trans_triG_eqtl$chr==7,c(1,4)]
loc2_aAGT=filt_trans_triG_aqtl[filt_trans_triG_aqtl$gene=="AGT" & filt_trans_triG_aqtl$chr==7,c(1,4)]
loc2_eRABIF=filt_trans_triG_eqtl[filt_trans_triG_eqtl$gene=="RABIF" & filt_trans_triG_eqtl$chr==7,c(1,4)]
loc2_aRABIF=filt_trans_triG_aqtl[filt_trans_triG_aqtl$gene=="RABIF" & filt_trans_triG_aqtl$chr==7,c(1,4)]
colnames(bmi2)=c("rsid","pval")
colnames(t2d2)=c("rsid","pval")
colnames(hdl2)=c("rsid","pval")
colnames(triG2)=c("rsid","pval")
colnames(loc2_eLINC)=c("rsid","pval")
colnames(loc2_aLINC)=c("rsid","pval")
colnames(loc2_eKLF14)=c("rsid","pval")
colnames(loc2_aKLF14)=c("rsid","pval")
colnames(loc2_eAC)=c("rsid","pval")
colnames(loc2_eTBX4)=c("rsid","pval")
colnames(loc2_aTBX4)=c("rsid","pval")
colnames(loc2_eGNB1)=c("rsid","pval")
colnames(loc2_aGNB1)=c("rsid","pval")
colnames(loc2_eESR2)=c("rsid","pval")
colnames(loc2_aESR2)=c("rsid","pval")
colnames(loc2_eNR2F1)=c("rsid","pval")
colnames(loc2_aNR2F1)=c("rsid","pval")
colnames(loc2_eAGT)=c("rsid","pval")
colnames(loc2_aAGT)=c("rsid","pval")
colnames(loc2_eRABIF)=c("rsid","pval")
colnames(loc2_aRABIF)=c("rsid","pval")
rownames(bmi2)=bmi2$rsid
rownames(t2d2)=t2d2$rsid
rownames(hdl2)=hdl2$rsid
rownames(triG2)=triG2$rsid
rownames(loc2_eLINC)=loc2_eLINC$rsid
rownames(loc2_aLINC)=loc2_aLINC$rsid
rownames(loc2_eKLF14)=loc2_eKLF14$rsid
rownames(loc2_aKLF14)=loc2_aKLF14$rsid
rownames(loc2_eAC)=loc2_eAC$rsid
rownames(loc2_eTBX4)=loc2_eTBX4$rsid
rownames(loc2_aTBX4)=loc2_aTBX4$rsid
rownames(loc2_eGNB1)=loc2_eGNB1$rsid
rownames(loc2_aGNB1)=loc2_aGNB1$rsid
rownames(loc2_eESR2)=loc2_eESR2$rsid
rownames(loc2_aESR2)=loc2_aESR2$rsid
rownames(loc2_eNR2F1)=loc2_eNR2F1$rsid
rownames(loc2_aNR2F1)=loc2_aNR2F1$rsid
rownames(loc2_eAGT)=loc2_eAGT$rsid
rownames(loc2_aAGT)=loc2_aAGT$rsid
rownames(loc2_eRABIF)=loc2_eRABIF$rsid
rownames(loc2_aRABIF)=loc2_aRABIF$rsid

# Check out some relevant LocusCompare plots before picking the which to write to file
locuscompare(in_fn1=bmi2,in_fn2=t2d2,title1 = "BMI GWAS", title2 = "T2D GWAS",snp = "rs972283") # Top BMI GWAS SNP
locuscompare(in_fn1=bmi2,in_fn2=t2d2,title1 = "BMI GWAS", title2 = "T2D GWAS",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
locuscompare(in_fn1=bmi2,in_fn2=hdl2,title1 = "BMI GWAS", title2 = "HDL GWAS",snp = "rs972283") # Top BMI GWAS SNP
locuscompare(in_fn1=bmi2,in_fn2=hdl2,title1 = "BMI GWAS", title2 = "HDL GWAS",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=bmi2,in_fn2=triG2,title1 = "BMI GWAS", title2 = "Triglycerides GWAS",snp = "rs972283") # Top BMI GWAS SNP
locuscompare(in_fn1=bmi2,in_fn2=triG2,title1 = "BMI GWAS", title2 = "Triglycerides GWAS",snp = "rs287621") # Top TriG GWAS SNP
locuscompare(in_fn1=t2d2,in_fn2=hdl2,title1 = "T2D GWAS", title2 = "HDL GWAS",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
locuscompare(in_fn1=t2d2,in_fn2=hdl2,title1 = "T2D GWAS", title2 = "HDL GWAS",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=t2d2,in_fn2=triG2,title1 = "T2D GWAS", title2 = "Triglycerides GWAS",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
locuscompare(in_fn1=t2d2,in_fn2=triG2,title1 = "T2D GWAS", title2 = "Triglycerides GWAS",snp = "rs287621") # Top TriG GWAS SNP
locuscompare(in_fn1=hdl2,in_fn2=triG2,title1 = "HDL GWAS", title2 = "Triglycerides GWAS",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=hdl2,in_fn2=triG2,title1 = "HDL GWAS", title2 = "Triglycerides GWAS",snp = "rs287621") # Top TriG GWAS SNP
locuscompare(in_fn1=bmi2,in_fn2=loc2_eLINC,title1 = "BMI GWAS", title2 = "LINC-PINT cis-eQTL",snp = "rs972283") # Top BMI GWAS SNP
locuscompare(in_fn1=loc2_eLINC,in_fn2=loc2_aLINC,title1 = "LINC-PINT cis-eQTL", title2 = "LINC-PINT cis-aQTL",snp = "rs972283") # Top BMI GWAS SNP
locuscompare(in_fn1=t2d2,in_fn2=loc2_eLINC,title1 = "T2D GWAS", title2 = "LINC-PINT cis-eQTL",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
locuscompare(in_fn1=hdl2,in_fn2=loc2_eLINC,title1 = "HDL GWAS", title2 = "LINC-PINT cis-eQTL",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=triG2,in_fn2=loc2_eLINC,title1 = "Triglycerides GWAS", title2 = "LINC-PINT cis-eQTL",snp = "rs287621") # Top TriG GWAS SNP
locuscompare(in_fn1=bmi2,in_fn2=loc2_eKLF14,title1 = "BMI GWAS", title2 = "KLF14 cis-eQTL",snp = "rs972283") # Top BMI GWAS SNP
locuscompare(in_fn1=t2d2,in_fn2=loc2_eKLF14,title1 = "T2D GWAS", title2 = "KLF14 cis-eQTL",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
locuscompare(in_fn1=hdl2,in_fn2=loc2_eKLF14,title1 = "HDL GWAS", title2 = "KLF14 cis-eQTL",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=triG2,in_fn2=loc2_eKLF14,title1 = "Triglycerides GWAS", title2 = "KLF14 cis-eQTL",snp = "rs287621") # Top TriG GWAS SNP
locuscompare(in_fn1=bmi2,in_fn2=loc2_aKLF14,title1 = "BMI GWAS", title2 = "KLF14 cis-aQTL",snp = "rs972283") # Top BMI GWAS SNP
locuscompare(in_fn1=t2d2,in_fn2=loc2_aKLF14,title1 = "T2D GWAS", title2 = "KLF14 cis-aQTL",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
locuscompare(in_fn1=hdl2,in_fn2=loc2_aKLF14,title1 = "HDL GWAS", title2 = "KLF14 cis-aQTL",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=triG2,in_fn2=loc2_aKLF14,title1 = "Triglycerides GWAS", title2 = "KLF14 cis-aQTL",snp = "rs287621") # Top TriG GWAS SNP
locuscompare(in_fn1=loc2_eKLF14,in_fn2=loc2_aKLF14,title1 = "KLF14 cis-eQTL", title2 = "KLF14 cis-aQTL",snp = "rs972283") # Top BMI GWAS SNP
locuscompare(in_fn1=loc2_eKLF14,in_fn2=loc2_aKLF14,title1 = "KLF14 cis-eQTL", title2 = "KLF14 cis-aQTL",snp = "rs738134") # Top BMI GWAS SNP
locuscompare(in_fn1=loc2_eKLF14,in_fn2=loc2_aKLF14,title1 = "KLF14 cis-eQTL", title2 = "KLF14 cis-aQTL",snp = "rs287621") # Top BMI GWAS SNP
locuscompare(in_fn1=bmi2,in_fn2=loc2_eAC,title1 = "BMI GWAS", title2 = "AC016831.7 cis-eQTL",snp = "rs972283") # Top BMI GWAS SNP
locuscompare(in_fn1=t2d2,in_fn2=loc2_eAC,title1 = "T2D GWAS", title2 = "AC016831.7 cis-eQTL",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
locuscompare(in_fn1=hdl2,in_fn2=loc2_eAC,title1 = "HDL GWAS", title2 = "AC016831.7 cis-eQTL",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=triG2,in_fn2=loc2_eAC,title1 = "Triglycerides GWAS", title2 = "AC016831.7 cis-eQTL",snp = "rs287621") # Top TriG GWAS SNP
locuscompare(in_fn1=bmi2,in_fn2=loc2_eTBX4,title1 = "BMI GWAS", title2 = "TBX4 trans-eQTL",snp = "rs972283") # Top BMI GWAS SNP
locuscompare(in_fn1=t2d2,in_fn2=loc2_eTBX4,title1 = "T2D GWAS", title2 = "TBX4 trans-eQTL",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
locuscompare(in_fn1=hdl2,in_fn2=loc2_eTBX4,title1 = "HDL GWAS", title2 = "TBX4 trans-eQTL",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=triG2,in_fn2=loc2_eTBX4,title1 = "Triglycerides GWAS", title2 = "TBX4 trans-eQTL",snp = "rs287621") # Top HDL GWAS SNP
locuscompare(in_fn1=bmi2,in_fn2=loc2_aTBX4,title1 = "BMI GWAS", title2 = "TBX4 trans-aQTL",snp = "rs972283") # Top BMI GWAS SNP
locuscompare(in_fn1=t2d2,in_fn2=loc2_aTBX4,title1 = "T2D GWAS", title2 = "TBX4 trans-aQTL",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
locuscompare(in_fn1=hdl2,in_fn2=loc2_aTBX4,title1 = "HDL GWAS", title2 = "TBX4 trans-aQTL",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=triG2,in_fn2=loc2_aTBX4,title1 = "Triglycerides GWAS", title2 = "TBX4 trans-aQTL",snp = "rs287621") # Top HDL GWAS SNP
locuscompare(in_fn1=t2d2,in_fn2=loc2_eGNB1,title1 = "T2D GWAS", title2 = "GNB1 trans-eQTL",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
locuscompare(in_fn1=t2d2,in_fn2=loc2_aGNB1,title1 = "T2D GWAS", title2 = "GNB1 trans-aQTL",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
locuscompare(in_fn1=hdl2,in_fn2=loc2_eESR2,title1 = "HDL GWAS", title2 = "ESR2 trans-eQTL",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=hdl2,in_fn2=loc2_aESR2,title1 = "HDL GWAS", title2 = "ESR2 trans-aQTL",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=hdl2,in_fn2=loc2_eNR2F1,title1 = "HDL GWAS", title2 = "NR2F1 trans-eQTL",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=hdl2,in_fn2=loc2_aNR2F1,title1 = "HDL GWAS", title2 = "NR2F1 trans-aQTL",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=loc2_eNR2F1,in_fn2=loc2_aNR2F1,title1 = "NR2F1 trans-eQTL", title2 = "NR2F1 trans-aQTL",snp = "rs11765979") # Top HDL GWAS SNP
locuscompare(in_fn1=triG2,in_fn2=loc2_eRABIF,title1 = "Triglycerides GWAS", title2 = "RABIF trans-eQTL",snp = "rs287621") # Top TriG GWAS SNP
locuscompare(in_fn1=triG2,in_fn2=loc2_aRABIF,title1 = "Triglycerides GWAS", title2 = "RABIF trans-aQTL",snp = "rs287621") # Top TriG GWAS SNP
locuscompare(in_fn1=triG2,in_fn2=loc2_eAGT,title1 = "Triglycerides GWAS", title2 = "AGT trans-eQTL",snp = "rs287621") # Top TriG GWAS SNP
locuscompare(in_fn1=triG2,in_fn2=loc2_aAGT,title1 = "Triglycerides GWAS", title2 = "AGT trans-aQTL",snp = "rs287621") # Top TriG GWAS SNP
locuscompare(in_fn1=loc2_eAGT,in_fn2=loc2_aAGT,title1 = "AGT trans-eQTL", title2 = "AGT trans-aQTL",snp = "rs287621") # Top TriG GWAS SNP
locuscompare(in_fn1=loc2_eAGT,in_fn2=loc2_aAGT,title1 = "AGT trans-eQTL", title2 = "AGT trans-aQTL",snp = "rs972283") # Top TriG GWAS SNP
locuscompare(in_fn1=loc2_eAGT,in_fn2=loc2_aAGT,title1 = "AGT trans-eQTL", title2 = "AGT trans-aQTL",snp = "rs738134") # Top TriG GWAS SNP

# I never ran HyPrColoc on just cis-e_KLF14 and cis-a_KLF14 alone. First I need to calculate the SEs and format the data.
filt_cis_eqtl$SE=filt_cis_eqtl$beta/filt_cis_eqtl$statistic
filt_cis_aqtl$SE=filt_cis_aqtl$beta/filt_cis_aqtl$statistic
all(filt_cis_eqtl$snps[filt_cis_eqtl$gene=="KLF14"]==filt_cis_aqtl$snps[filt_cis_aqtl$gene=="KLF14"]) # TRUE
betas2=cbind("cis-e_KLF14"=filt_cis_eqtl[filt_cis_eqtl$gene=="KLF14","beta"],"cis-a_KLF14"=filt_cis_aqtl[filt_cis_aqtl$gene=="KLF14","beta"])
ses2=cbind("cis-e_KLF14"=filt_cis_eqtl[filt_cis_eqtl$gene=="KLF14","SE"],"cis-a_KLF14"=filt_cis_aqtl[filt_cis_aqtl$gene=="KLF14","SE"])
rownames(betas2)=filt_cis_aqtl$snps[filt_cis_aqtl$gene=="KLF14"]
rownames(ses2)=filt_cis_aqtl$snps[filt_cis_aqtl$gene=="KLF14"]
all(rownames(betas2)==rownames(filt_ld[[2]])) # TRUE
all(rownames(ses2)==rownames(filt_ld[[2]])) # TRUE
eKLF14_aKLF14=hyprcoloc(as.matrix(betas2),as.matrix(ses2),
                        trait.names=colnames(betas2),snp.id=rownames(betas2),ld.matrix = filt_ld[[2]],
                        trait.subset = c("cis-e_KLF14","cis-a_KLF14"),snpscores = T)
# KLF14 eQTL and aQTL colocalize with a PP=0.9094 that is best explained by rs4731702.

# Now let's do the same sort of analysis for NR2F1 and AGT.
filt_trans_hdl_eqtl$SE=filt_trans_hdl_eqtl$beta/filt_trans_hdl_eqtl$statistic
filt_trans_hdl_aqtl$SE=filt_trans_hdl_aqtl$beta/filt_trans_hdl_aqtl$statistic
temp_e=filt_trans_hdl_eqtl[filt_trans_hdl_eqtl$chr==7 & filt_trans_hdl_eqtl$gene=="NR2F1",]
temp_a=filt_trans_hdl_aqtl[filt_trans_hdl_aqtl$chr==7 & filt_trans_hdl_aqtl$gene=="NR2F1",]
all(rownames(betas2)==temp_e$snps) # TRUE
all(rownames(betas2)==temp_a$snps) # TRUE
betas2=cbind(betas2,"trans-e_NR2F1"=temp_e$beta,"trans-a_NR2F1"=temp_a$beta)
ses2=cbind(ses2,"trans-e_NR2F1"=temp_e$SE,"trans-a_NR2F1"=temp_a$SE)

filt_trans_triG_eqtl$SE=filt_trans_triG_eqtl$beta/filt_trans_triG_eqtl$statistic
filt_trans_triG_aqtl$SE=filt_trans_triG_aqtl$beta/filt_trans_triG_aqtl$statistic
temp_e=filt_trans_triG_eqtl[filt_trans_triG_eqtl$chr==7 & filt_trans_triG_eqtl$gene=="AGT",]
temp_a=filt_trans_triG_aqtl[filt_trans_triG_aqtl$chr==7 & filt_trans_triG_aqtl$gene=="AGT",]
all(rownames(betas2)==temp_e$snps) # TRUE
all(rownames(betas2)==temp_a$snps) # FALSE
temp_a=temp_a[match(rownames(betas2),temp_a$snps),]
all(rownames(betas2)==temp_a$snps) # TRUE
betas2=cbind(betas2,"trans-e_AGT"=temp_e$beta,"trans-a_AGT"=temp_a$beta)
ses2=cbind(ses2,"trans-e_AGT"=temp_e$SE,"trans-a_AGT"=temp_a$SE)

eNR2F1_aNR2F1=hyprcoloc(as.matrix(betas2),as.matrix(ses2),
                        trait.names=colnames(betas2),snp.id=rownames(betas2),ld.matrix = filt_ld[[2]],
                        trait.subset = c("trans-e_NR2F1","trans-a_NR2F1"),snpscores = T)
# NR2F1 eQTL and aQTL colocalize with a PP=0.8700 that is best explained by rs738134.

eAGT_aAGT=hyprcoloc(as.matrix(betas2),as.matrix(ses2),
                        trait.names=colnames(betas2),snp.id=rownames(betas2),ld.matrix = filt_ld[[2]],
                        trait.subset = c("trans-e_AGT","trans-a_AGT"),snpscores = T)
# AGT eQTL and aQTL colocalize with a PP=0.6451 that is best explained by rs11765979.

# Let's write some to PDFs
pdf("BMI_T2D_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=bmi2,in_fn2=t2d2,title1 = "BMI GWAS", title2 = "T2D GWAS",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("BMI_T2D_7q32_LocusCompare_rs287621.pdf",width = 10)
locuscompare(in_fn1=bmi2,in_fn2=t2d2,title1 = "BMI GWAS", title2 = "T2D GWAS",snp = "rs287621") # Top TriG GWAS SNP
dev.off()
pdf("BMI_T2D_7q32_LocusCompare_rs738134.pdf",width = 10)
locuscompare(in_fn1=bmi2,in_fn2=t2d2,title1 = "BMI GWAS", title2 = "T2D GWAS",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
dev.off()
pdf("BMI_HDL_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=bmi2,in_fn2=hdl2,title1 = "BMI GWAS", title2 = "HDL GWAS",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("BMI_HDL_7q32_LocusCompare_rs287621.pdf",width = 10)
locuscompare(in_fn1=bmi2,in_fn2=hdl2,title1 = "BMI GWAS", title2 = "HDL GWAS",snp = "rs287621") # Top TriG GWAS SNP
dev.off()
pdf("BMI_HDL_7q32_LocusCompare_rs11765979.pdf",width = 10)
locuscompare(in_fn1=bmi2,in_fn2=hdl2,title1 = "BMI GWAS", title2 = "HDL GWAS",snp = "rs11765979") # Top HDL GWAS SNP
dev.off()
pdf("BMI_TriG_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=bmi2,in_fn2=triG2,title1 = "BMI GWAS", title2 = "Triglycerides GWAS",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("BMI_TriG_7q32_LocusCompare_rs287621.pdf",width = 10)
locuscompare(in_fn1=bmi2,in_fn2=triG2,title1 = "BMI GWAS", title2 = "Triglycerides GWAS",snp = "rs287621") # Top TriG GWAS SNP
dev.off()
pdf("T2D_HDL_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=t2d2,in_fn2=hdl2,title1 = "T2D GWAS", title2 = "HDL GWAS",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("T2D_HDL_7q32_LocusCompare_rs738134.pdf",width = 10)
locuscompare(in_fn1=t2d2,in_fn2=hdl2,title1 = "T2D GWAS", title2 = "HDL GWAS",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
dev.off()
pdf("T2D_HDL_7q32_LocusCompare_rs11765979.pdf",width = 10)
locuscompare(in_fn1=t2d2,in_fn2=hdl2,title1 = "T2D GWAS", title2 = "HDL GWAS",snp = "rs11765979") # Top HDL GWAS SNP
dev.off()
pdf("T2D_TriG_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=t2d2,in_fn2=triG2,title1 = "T2D GWAS", title2 = "Triglycerides GWAS",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("T2D_TriG_7q32_LocusCompare_rs738134.pdf",width = 10)
locuscompare(in_fn1=t2d2,in_fn2=triG2,title1 = "T2D GWAS", title2 = "Triglycerides GWAS",snp = "rs738134") # Near top T2D GWAS SNP (top SNP did not overlap other GWAS)
dev.off()
pdf("T2D_TriG_7q32_LocusCompare_rs287621.pdf",width = 10)
locuscompare(in_fn1=t2d2,in_fn2=triG2,title1 = "T2D GWAS", title2 = "Triglycerides GWAS",snp = "rs287621") # Top TriG GWAS SNP
dev.off()
pdf("HDL_TriG_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=hdl2,in_fn2=triG2,title1 = "HDL GWAS", title2 = "Triglycerides GWAS",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("HDL_TriG_7q32_LocusCompare_rs11765979.pdf",width = 10)
locuscompare(in_fn1=hdl2,in_fn2=triG2,title1 = "HDL GWAS", title2 = "Triglycerides GWAS",snp = "rs11765979") # Top HDL GWAS SNP
dev.off()
pdf("HDL_TriG_7q32_LocusCompare_rs287621.pdf",width = 10)
locuscompare(in_fn1=hdl2,in_fn2=triG2,title1 = "HDL GWAS", title2 = "Triglycerides GWAS",snp = "rs287621") # Top TriG GWAS SNP
dev.off()
pdf("./BMI/BMI_e_LINC-PINT_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=bmi2,in_fn2=loc2_eLINC,title1 = "BMI GWAS", title2 = "LINC-PINT cis-eQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./T2D/T2D_e_LINC-PINT_7q32_LocusCompare_rs738134.pdf",width = 10)
locuscompare(in_fn1=t2d2,in_fn2=loc2_eLINC,title1 = "T2D GWAS", title2 = "LINC-PINT cis-eQTL",snp = "rs738134") # Near top T2D GWAS SNP
dev.off()
pdf("./HDL/HDL_e_LINC-PINT_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=hdl2,in_fn2=loc2_eLINC,title1 = "HDL GWAS", title2 = "LINC-PINT cis-eQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./Triglycerides/TriG_e_LINC-PINT_7q32_LocusCompare_rs287621.pdf",width = 10)
locuscompare(in_fn1=triG2,in_fn2=loc2_eLINC,title1 = "Triglycerides GWAS", title2 = "LINC-PINT cis-eQTL",snp = "rs287621") # Top TriG GWAS SNP
dev.off()
pdf("./LINC-PINT_eQTL_aQTL_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=loc2_eLINC,in_fn2=loc2_aLINC,title1 = "LINC-PINT cis-eQTL", title2 = "LINC-PINT cis-aQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./BMI/BMI_e_AC016831.7_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=bmi2,in_fn2=loc2_eAC,title1 = "BMI GWAS", title2 = "AC016831.7 cis-eQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./T2D/T2D_e_AC016831.7_7q32_LocusCompare_rs738134.pdf",width = 10)
locuscompare(in_fn1=t2d2,in_fn2=loc2_eAC,title1 = "T2D GWAS", title2 = "AC016831.7 cis-eQTL",snp = "rs738134") # Near top T2D GWAS SNP
dev.off()
pdf("./HDL/HDL_e_AC016831.7_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=hdl2,in_fn2=loc2_eAC,title1 = "HDL GWAS", title2 = "AC016831.7 cis-eQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./Triglycerides/TriG_e_AC016831.7_7q32_LocusCompare_rs287621.pdf",width = 10)
locuscompare(in_fn1=triG2,in_fn2=loc2_eAC,title1 = "Triglycerides GWAS", title2 = "AC016831.7 cis-eQTL",snp = "rs287621") # Top TriG GWAS SNP
dev.off()
pdf("./BMI/BMI_e_KLF14_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=bmi2,in_fn2=loc2_eKLF14,title1 = "BMI GWAS", title2 = "KLF14 cis-eQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./T2D/T2D_e_KLF14_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=t2d2,in_fn2=loc2_eKLF14,title1 = "T2D GWAS", title2 = "KLF14 cis-eQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./T2D/T2D_e_KLF14_7q32_LocusCompare_rs738134.pdf",width = 10)
locuscompare(in_fn1=t2d2,in_fn2=loc2_eKLF14,title1 = "T2D GWAS", title2 = "KLF14 cis-eQTL",snp = "rs738134") # Near top T2D GWAS SNP
dev.off()
pdf("./HDL/HDL_e_KLF14_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=hdl2,in_fn2=loc2_eKLF14,title1 = "HDL GWAS", title2 = "KLF14 cis-eQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./Triglycerides/TriG_e_KLF14_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=triG2,in_fn2=loc2_eKLF14,title1 = "Triglycerides GWAS", title2 = "KLF14 cis-eQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./Triglycerides/TriG_e_KLF14_7q32_LocusCompare_rs287621.pdf",width = 10)
locuscompare(in_fn1=triG2,in_fn2=loc2_eKLF14,title1 = "Triglycerides GWAS", title2 = "KLF14 cis-eQTL",snp = "rs287621") # Top TriG GWAS SNP
dev.off()
pdf("./BMI/BMI_a_KLF14_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=bmi2,in_fn2=loc2_aKLF14,title1 = "BMI GWAS", title2 = "KLF14 cis-aQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./T2D/T2D_a_KLF14_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=t2d2,in_fn2=loc2_aKLF14,title1 = "T2D GWAS", title2 = "KLF14 cis-aQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./T2D/T2D_a_KLF14_7q32_LocusCompare_rs738134.pdf",width = 10)
locuscompare(in_fn1=t2d2,in_fn2=loc2_aKLF14,title1 = "T2D GWAS", title2 = "KLF14 cis-aQTL",snp = "rs738134") # Near top T2D GWAS SNP
dev.off()
pdf("./HDL/HDL_a_KLF14_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=hdl2,in_fn2=loc2_aKLF14,title1 = "HDL GWAS", title2 = "KLF14 cis-aQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./Triglycerides/TriG_a_KLF14_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=triG2,in_fn2=loc2_aKLF14,title1 = "Triglycerides GWAS", title2 = "KLF14 cis-aQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./Triglycerides/TriG_a_KLF14_7q32_LocusCompare_rs287621.pdf",width = 10)
locuscompare(in_fn1=triG2,in_fn2=loc2_aKLF14,title1 = "Triglycerides GWAS", title2 = "KLF14 cis-aQTL",snp = "rs287621") # Top TriG GWAS SNP
dev.off()
pdf("./KLF14_eQTL_aQTL_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=loc2_eKLF14,in_fn2=loc2_aKLF14,title1 = "KLF14 cis-eQTL", title2 = "KLF14 cis-aQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./BMI/BMI_e_TBX4_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=bmi2,in_fn2=loc2_eTBX4,title1 = "BMI GWAS", title2 = "TBX4 trans-eQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./T2D/T2D_e_TBX4_7q32_LocusCompare_rs738134.pdf",width = 10)
locuscompare(in_fn1=t2d2,in_fn2=loc2_eTBX4,title1 = "T2D GWAS", title2 = "TBX4 trans-eQTL",snp = "rs738134") # Near top T2D GWAS SNP
dev.off()
pdf("./HDL/HDL_e_TBX4_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=hdl2,in_fn2=loc2_eTBX4,title1 = "HDL GWAS", title2 = "TBX4 trans-eQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./Triglycerides/TriG_e_TBX4_7q32_LocusCompare_rs287621.pdf",width = 10)
locuscompare(in_fn1=triG2,in_fn2=loc2_eTBX4,title1 = "Triglycerides GWAS", title2 = "TBX4 trans-eQTL",snp = "rs287621") # Top TriG GWAS SNP
dev.off()
pdf("./TBX4_eQTL_aQTL_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=loc2_eTBX4,in_fn2=loc2_aTBX4,title1 = "TBX4 cis-eQTL", title2 = "TBX4 cis-aQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./HDL/HDL_e_NR2F1_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=hdl2,in_fn2=loc2_eNR2F1,title1 = "HDL GWAS", title2 = "NR2F1 trans-eQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./HDL/HDL_a_NR2F1_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=hdl2,in_fn2=loc2_aNR2F1,title1 = "HDL GWAS", title2 = "NR2F1 trans-aQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./NR2F1_eQTL_aQTL_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=loc2_eNR2F1,in_fn2=loc2_aNR2F1,title1 = "NR2F1 trans-eQTL", title2 = "NR2F1 trans-aQTL",snp = "rs972283") # Top BMI GWAS SNP
dev.off()
pdf("./Triglycerides/TriG_e_AGT_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=triG2,in_fn2=loc2_eAGT,title1 = "Triglycerides GWAS", title2 = "AGT trans-eQTL",snp = "rs287621") # Top BMI GWAS SNP
dev.off()
pdf("./Triglycerides/TriG_a_AGT_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=triG2,in_fn2=loc2_aAGT,title1 = "Triglycerides GWAS", title2 = "AGT trans-aQTL",snp = "rs287621") # Top BMI GWAS SNP
dev.off()
pdf("./AGT_eQTL_aQTL_7q32_LocusCompare_rs972283.pdf",width = 10)
locuscompare(in_fn1=loc2_eAGT,in_fn2=loc2_aAGT,title1 = "AGT trans-eQTL", title2 = "AGT trans-aQTL",snp = "rs287621") # Top BMI GWAS SNP
dev.off()

## 12p13.1
# First, grab the necessary P-values for the SNPs used in the HyPrColoc analyses for the traits of interest
bmi4=filt_bmi[match(rownames(filt_ld[[4]]),filt_bmi$SNP),c(3,9)]
loc4_eANG=filt_trans_bmi_eqtl[filt_trans_bmi_eqtl$gene=="ANG",c(1,4)]
loc4_aANG=filt_trans_bmi_aqtl[filt_trans_bmi_aqtl$gene=="ANG",c(1,4)]
loc4_eID2=filt_trans_bmi_eqtl[filt_trans_bmi_eqtl$gene=="ID2",c(1,4)]
loc4_aID2=filt_trans_bmi_aqtl[filt_trans_bmi_aqtl$gene=="ID2",c(1,4)]
loc4_ePTPRJ=filt_trans_bmi_eqtl[filt_trans_bmi_eqtl$gene=="PTPRJ",c(1,4)]
loc4_aPTPRJ=filt_trans_bmi_aqtl[filt_trans_bmi_aqtl$gene=="PTPRJ",c(1,4)]
loc4_eTENM4=filt_trans_bmi_eqtl[filt_trans_bmi_eqtl$gene=="TENM4",c(1,4)]
loc4_aTENM4=filt_trans_bmi_aqtl[filt_trans_bmi_aqtl$gene=="TENM4",c(1,4)]
loc4_eEPHB2=filt_trans_bmi_eqtl[filt_trans_bmi_eqtl$gene=="EPHB2",c(1,4)]
loc4_aEPHB2=filt_trans_bmi_aqtl[filt_trans_bmi_aqtl$gene=="EPHB2",c(1,4)]
colnames(bmi4)=c("rsid","pval")
colnames(loc4_eANG)=c("rsid","pval")
colnames(loc4_aANG)=c("rsid","pval")
colnames(loc4_eID2)=c("rsid","pval")
colnames(loc4_aID2)=c("rsid","pval")
colnames(loc4_ePTPRJ)=c("rsid","pval")
colnames(loc4_aPTPRJ)=c("rsid","pval")
colnames(loc4_eTENM4)=c("rsid","pval")
colnames(loc4_aTENM4)=c("rsid","pval")
colnames(loc4_eEPHB2)=c("rsid","pval")
colnames(loc4_aEPHB2)=c("rsid","pval")
rownames(bmi4)=bmi4$rsid
rownames(loc4_eANG)=loc4_eANG$rsid
rownames(loc4_aANG)=loc4_aANG$rsid
rownames(loc4_eID2)=loc4_eID2$rsid
rownames(loc4_aID2)=loc4_aID2$rsid
rownames(loc4_ePTPRJ)=loc4_ePTPRJ$rsid
rownames(loc4_aPTPRJ)=loc4_aPTPRJ$rsid
rownames(loc4_eTENM4)=loc4_eTENM4$rsid
rownames(loc4_aTENM4)=loc4_aTENM4$rsid
rownames(loc4_eEPHB2)=loc4_eEPHB2$rsid
rownames(loc4_aEPHB2)=loc4_aEPHB2$rsid

# Check out some relevant LocusCompare plots before picking the which to write to file
locuscompare(in_fn1=bmi4,in_fn2=loc4_eEPHB2,title1 = "BMI GWAS", title2 = "EPHB2 trans-eQTL",snp = "rs12422552") # Top GWAS SNP
locuscompare(in_fn1=bmi4,in_fn2=loc4_aEPHB2,title1 = "BMI GWAS", title2 = "EPHB2 trans-aQTL",snp = "rs12422552") # Top GWAS SNP
locuscompare(in_fn1=bmi4,in_fn2=loc4_eANG,title1 = "BMI GWAS", title2 = "ANG trans-eQTL",snp = "rs12422552") # Top GWAS SNP
locuscompare(in_fn1=bmi4,in_fn2=loc4_aANG,title1 = "BMI GWAS", title2 = "ANG trans-aQTL",snp = "rs12422552") # Top GWAS SNP
locuscompare(in_fn1=bmi4,in_fn2=loc4_eID2,title1 = "BMI GWAS", title2 = "ID2 trans-eQTL",snp = "rs12422552") # Top GWAS SNP
locuscompare(in_fn1=bmi4,in_fn2=loc4_aID2,title1 = "BMI GWAS", title2 = "ID2 trans-aQTL",snp = "rs12422552") # Top GWAS SNP
locuscompare(in_fn1=bmi4,in_fn2=loc4_ePTPRJ,title1 = "BMI GWAS", title2 = "PTPRJ trans-eQTL",snp = "rs12422552") # Top GWAS SNP
locuscompare(in_fn1=bmi4,in_fn2=loc4_aPTPRJ,title1 = "BMI GWAS", title2 = "PTPRJ trans-aQTL",snp = "rs12422552") # Top GWAS SNP
locuscompare(in_fn1=bmi4,in_fn2=loc4_eTENM4,title1 = "BMI GWAS", title2 = "TENM4 trans-eQTL",snp = "rs12422552") # Top GWAS SNP
locuscompare(in_fn1=bmi4,in_fn2=loc4_aTENM4,title1 = "BMI GWAS", title2 = "TENM4 trans-aQTL",snp = "rs12422552") # Top GWAS SNP

# Let's write some to PDFs
pdf("./BMI/rs12422552-ANG_eQTL_and_BMI_12p13.1_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi4,in_fn2=loc4_eANG,title1 = "BMI GWAS", title2 = "ANG cis-eQTL",snp = "rs12422552")
dev.off()
pdf("./BMI/rs12422552-ANG_aQTL_and_BMI_12p13.1_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi4,in_fn2=loc4_aANG,title1 = "BMI GWAS", title2 = "ANG cis-aQTL",snp = "rs12422552")
dev.off()
pdf("./BMI/rs12422552-ID2_eQTL_and_BMI_12p13.1_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi4,in_fn2=loc4_eID2,title1 = "BMI GWAS", title2 = "ID2 cis-eQTL",snp = "rs12422552")
dev.off()
pdf("./BMI/rs12422552-ID2_aQTL_and_BMI_12p13.1_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi4,in_fn2=loc4_aID2,title1 = "BMI GWAS", title2 = "ID2 cis-aQTL",snp = "rs12422552")
dev.off()
pdf("./BMI/rs12422552-PTPRJ_eQTL_and_BMI_12p13.1_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi4,in_fn2=loc4_ePTPRJ,title1 = "BMI GWAS", title2 = "PTPRJ cis-eQTL",snp = "rs12422552")
dev.off()
pdf("./BMI/rs12422552-PTPRJ_aQTL_and_BMI_12p13.1_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi4,in_fn2=loc4_aPTPRJ,title1 = "BMI GWAS", title2 = "PTPRJ cis-aQTL",snp = "rs12422552")
dev.off()
pdf("./BMI/rs12422552-TENM4_eQTL_and_BMI_12p13.1_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi4,in_fn2=loc4_eTENM4,title1 = "BMI GWAS", title2 = "TENM4 cis-eQTL",snp = "rs12422552")
dev.off()
pdf("./BMI/rs12422552-TENM4_aQTL_and_BMI_12p13.1_LocusCompare.pdf",width = 10)
locuscompare(in_fn1=bmi4,in_fn2=loc4_aTENM4,title1 = "BMI GWAS", title2 = "TENM4 cis-aQTL",snp = "rs12422552")
dev.off()

### Let's switch to pulling out data for Cytoscape network visualizations
# Read in data
bmi_pairColoc=read.table("./BMI/Pairwise_HyPrColoc_between_BMI_and_each_QTL_for_select_loci.txt",sep = "\t",header = T)
t2d_pairColoc=read.table("./T2D/Pairwise_HyPrColoc_between_BMIadjT2D_and_each_QTL_for_select_loci.txt",sep = "\t",header = T)
hdl_pairColoc=read.table("./HDL/Pairwise_HyPrColoc_between_HDL_and_each_QTL_for_select_loci.txt",sep = "\t",header = T)
triG_pairColoc=read.table("./Triglycerides/Pairwise_HyPrColoc_between_TriG_and_each_QTL_for_select_loci.txt",sep = "\t",header = T)
ephb2_pairColoc=read.table("./BMI/Pairwise_HyPrColoc_between_EPHB2_aQTL_and_each_other_QTL_for_1p36.txt",sep = "\t",header = T)
bmi_mrs=read.table("./BMI/Eurobats_adipose_time-matched_BMI_MRs_from_RF_modeling.txt",header = F)
homair_mrs=read.table("./HOMA-IR/Eurobats_adipose_time-matched_HOMA-IR_MRs_from_RF_modeling.txt",header = F)
hdl_mrs=read.table("./HDL/Eurobats_adipose_time-matched_HDL_MRs_from_RF_modeling.txt",header = F)
triG_mrs=read.table("./Triglycerides/Eurobats_adipose_time-matched_Triglycerides_MRs_from_RF_modeling.txt",header = F)
interactome=read.table("../Adipose expression data/FINAL_logTPMs_and_activities/Eurobats_adipose_900boots_regulon_with_LINC-PINT.txt",sep = "\t",header = T)
tpm=read.table("../Adipose expression data/FINAL_logTPMs_and_activities/Filtered_Eurobats_adipose_qnorm_INT_logTPMs_for_all_expressed_genes.txt",
               sep = "\t",header = T,row.names = 1)
vip=read.table("../Adipose expression data/FINAL_logTPMs_and_activities/Filtered_Eurobats_adipose_unnormalized_activities_from_logTPM_for_4213_regulators.txt",
               sep = "\t",header = T,row.names = 1)
phenos=read.table("../Eurobats phenotypes/Amendment_time-matched_phenotypes_E886_02082019_with_HOMA.txt",sep="\t",header = T,row.names = 1)
filt_pheno=phenos[na.omit(match(colnames(vip),rownames(phenos))),]
all(colnames(vip)==colnames(tpm)) # TRUE
all(colnames(vip)==rownames(filt_pheno)) # TRUE
sig_bmi=filt_bmi[filt_bmi$P<=5E-8,]
sig_t2d=filt_t2d[filt_t2d$Pvalue<=5E-8,]
sig_hdl=filt_hdl[filt_hdl$P<=5E-8,]
sig_triG=filt_triG[filt_triG$P<=5E-8,]

# Grab relevant sub-interactomes
bmi_MRregs=interactome[interactome$Target %in% bmi_mrs[,1],]
bmi_MRMR=bmi_MRregs[bmi_MRregs$Regulator %in% bmi_mrs[,1],]
homair_MRregs=interactome[interactome$Target %in% homair_mrs[,1],]
homair_MRMR=homair_MRregs[homair_MRregs$Regulator %in% homair_mrs[,1],]
hdl_MRregs=interactome[interactome$Target %in% hdl_mrs[,1],]
hdl_MRMR=hdl_MRregs[hdl_MRregs$Regulator %in% hdl_mrs[,1],]
triG_MRregs=interactome[interactome$Target %in% triG_mrs[,1],]
triG_MRMR=triG_MRregs[triG_MRregs$Regulator %in% triG_mrs[,1],]

# 1p36
# Grab interactions between EPHB2 and MRs in adipose interactome
EPHB2mrs=bmi_MRregs[bmi_MRregs$Regulator=="EPHB2",]
mrsEPHB2=interactome[(interactome$Regulator %in% bmi_mrs[,1]) & (interactome$Target=="EPHB2"),]
interactome1p36=rbind(bmi_MRMR,EPHB2mrs,mrsEPHB2)

# Grab pairwise colocalizations with PP>0.5 between BMI and QTLs and EPHB2 aQTL and trans-QTLs
bmi_pairColoc1p36=bmi_pairColoc[bmi_pairColoc$locus=="1p36.1" & bmi_pairColoc$posterior_prob>0.5,c(2,3,5)]
bmi_pairColoc1p36$traits=gsub("BMI, ","",bmi_pairColoc1p36$traits)
bmi_pairColoc1p36=cbind("trait1"=rep("BMI",dim(bmi_pairColoc1p36)[1]),bmi_pairColoc1p36)
bmi_e_pairColoc1=bmi_pairColoc1p36[grepl("-e_",bmi_pairColoc1p36$traits),]
bmi_a_pairColoc1=bmi_pairColoc1p36[grepl("-a_",bmi_pairColoc1p36$traits),]
bmi_e_pairColoc1$traits=gsub(".*_","",bmi_e_pairColoc1$traits)
bmi_a_pairColoc1$traits=gsub(".*_","",bmi_a_pairColoc1$traits)
bmi_netPair1p36=rbind(bmi_e_pairColoc1,bmi_a_pairColoc1)
bmi_netPair1p36=bmi_netPair1p36[!duplicated(bmi_netPair1p36$traits),-c(3,4)]
bmi_netPair1p36$eQTL_PP=bmi_e_pairColoc1[match(bmi_netPair1p36$traits,bmi_e_pairColoc1$traits),3]
bmi_netPair1p36$eQTL_SNP=bmi_e_pairColoc1[match(bmi_netPair1p36$traits,bmi_e_pairColoc1$traits),4]
bmi_netPair1p36$aQTL_PP=bmi_a_pairColoc1[match(bmi_netPair1p36$traits,bmi_a_pairColoc1$traits),3]
bmi_netPair1p36$aQTL_SNP=bmi_a_pairColoc1[match(bmi_netPair1p36$traits,bmi_a_pairColoc1$traits),4]

ephb2Coloc1=ephb2_pairColoc[!is.na(ephb2_pairColoc$candidate_snp),c(2,3,5)]
ephb2Coloc1=ephb2Coloc1[ephb2Coloc1$posterior_prob>0.5,]
ephb2Coloc1$traits=gsub("cis-a_EPHB2, ","",ephb2Coloc1$traits)
ephb2Coloc1=cbind("trait1"=rep("EPHB2",dim(ephb2Coloc1)[1]),ephb2Coloc1)
e_ephb2Coloc1=ephb2Coloc1[grepl("-e_",ephb2Coloc1$traits),]
a_ephb2Coloc1=ephb2Coloc1[grepl("-a_",ephb2Coloc1$traits),]
e_ephb2Coloc1$traits=gsub(".*_","",e_ephb2Coloc1$traits)
a_ephb2Coloc1$traits=gsub(".*_","",a_ephb2Coloc1$traits)
netEPHB2=rbind(e_ephb2Coloc1,a_ephb2Coloc1)
netEPHB2=netEPHB2[!duplicated(netEPHB2$traits),-c(3,4)]
netEPHB2$eQTL_PP=e_ephb2Coloc1[match(netEPHB2$traits,e_ephb2Coloc1$traits),3]
netEPHB2$eQTL_SNP=e_ephb2Coloc1[match(netEPHB2$traits,e_ephb2Coloc1$traits),4]
netEPHB2$aQTL_PP=a_ephb2Coloc1[match(netEPHB2$traits,a_ephb2Coloc1$traits),3]
netEPHB2$aQTL_SNP=a_ephb2Coloc1[match(netEPHB2$traits,a_ephb2Coloc1$traits),4]

colocNet1p36=rbind(bmi_netPair1p36,netEPHB2)
colocNet1p36[is.na(colocNet1p36)]=0

# Grab the -log10(Pmin) and betas for the eQTLs and aQTLs among BMI GWAS significant SNPs at the 1p36.1 locus
# Actually, though I originally made networks with nodes shaded according to their minP QTLs, I've since decided
# that for the manuscript I need to stick to a single SNP for all QTLs to avoid allele switching issues and to
# facilitate discussion in the manuscript. Therefore, for this locus I will focus on rs4654828 since it tends to
# be among the top SNPs for EPHB2 cis-aQTL and all BMI MR trans-aQTLs. However, since subsequent lines of code
# refer to the variables as min_cisE1, etc. I will keep that naming even though it's not an adequate description.
cisE1=filt_cis_eqtl[filt_cis_eqtl$chr==1,]
cisE1=cisE1[cisE1$snps %in% sig_bmi$SNP,]
cisE1=cisE1[order(cisE1$pvalue),]
#min_cisE1=cisE1[!duplicated(cisE1$gene),]
min_cisE1=cisE1[cisE1$snps=="rs4654828",]

cisA1=filt_cis_aqtl[filt_cis_aqtl$chr==1,]
cisA1=cisA1[cisA1$snps %in% sig_bmi$SNP,]
cisA1=cisA1[order(cisA1$pvalue),]
#min_cisA1=cisA1[!duplicated(cisA1$gene),]
min_cisA1=cisA1[cisA1$snps=="rs4654828",]

transE1=filt_trans_bmi_eqtl[filt_trans_bmi_eqtl$chr==1,]
transE1=transE1[transE1$snps %in% sig_bmi$SNP,]
transE1=transE1[order(transE1$pvalue),]
#min_transE1=transE1[!duplicated(transE1$gene),]
min_transE1=transE1[transE1$snps=="rs4654828",]

transA1=filt_trans_bmi_aqtl[filt_trans_bmi_aqtl$chr==1,]
transA1=transA1[transA1$snps %in% sig_bmi$SNP,]
transA1=transA1[order(transA1$pvalue),]
#min_transA1=transA1[!duplicated(transA1$gene),]
min_transA1=transA1[transA1$snps=="rs4654828",]

# Make node tables for 1p36 networks
inter_nodes1p36=data.frame("Node"=as.character(unique(interactome1p36$Regulator)),"BMI_exp_cor"=rep(0,length(unique(interactome1p36$Regulator))),
                           "BMI_act_cor"=rep(0,length(unique(interactome1p36$Regulator))),"rs4654828_eQTL_Beta"=rep(0,length(unique(interactome1p36$Regulator))),
                           "rs4654828_eQTL_logP"=rep(0,length(unique(interactome1p36$Regulator))),"rs4654828_aQTL_Beta"=rep(0,length(unique(interactome1p36$Regulator))),
                           "rs4654828_aQTL_logP"=rep(0,length(unique(interactome1p36$Regulator))))
for(i in 1:dim(inter_nodes1p36)[1]){
  inter_nodes1p36$BMI_exp_cor[i]=cor(as.numeric(tpm[as.character(inter_nodes1p36$Node[i]),]),filt_pheno$BMI)
  inter_nodes1p36$BMI_act_cor[i]=cor(as.numeric(vip[as.character(inter_nodes1p36$Node[i]),]),filt_pheno$BMI)
  inter_nodes1p36$rs4654828_eQTL_Beta[i]=ifelse(inter_nodes1p36$Node[i] %in% min_transE1$gene,
                                           min_transE1[min_transE1$gene==as.character(inter_nodes1p36$Node[i]),"beta"],
                                           0)
  inter_nodes1p36$rs4654828_eQTL_logP[i]=ifelse(inter_nodes1p36$Node[i] %in% min_transE1$gene,
                                           -log10(min_transE1[min_transE1$gene==as.character(inter_nodes1p36$Node[i]),"pvalue"]),
                                           0)
  inter_nodes1p36$rs4654828_aQTL_Beta[i]=ifelse(inter_nodes1p36$Node[i] %in% min_transA1$gene,
                                           min_transA1[min_transA1$gene==as.character(inter_nodes1p36$Node[i]),"beta"],
                                           0)
  inter_nodes1p36$rs4654828_aQTL_logP[i]=ifelse(inter_nodes1p36$Node[i] %in% min_transA1$gene,
                                           -log10(min_transA1[min_transA1$gene==as.character(inter_nodes1p36$Node[i]),"pvalue"]),
                                           0)
}

coloc_nodes1p36=data.frame("Node"=as.character(unique(colocNet1p36$traits)),"BMI_exp_cor"=rep(0,length(unique(colocNet1p36$traits))),
                           "BMI_act_cor"=rep(0,length(unique(colocNet1p36$traits))),"rs4654828_eQTL_Beta"=rep(0,length(unique(colocNet1p36$traits))),
                           "rs4654828_eQTL_logP"=rep(0,length(unique(colocNet1p36$traits))),"rs4654828_aQTL_Beta"=rep(0,length(unique(colocNet1p36$traits))),
                           "rs4654828_aQTL_logP"=rep(0,length(unique(colocNet1p36$traits))))
for(i in 1:dim(coloc_nodes1p36)[1]){
  coloc_nodes1p36$BMI_exp_cor[i]=cor(as.numeric(tpm[as.character(coloc_nodes1p36$Node[i]),]),filt_pheno$BMI)
  coloc_nodes1p36$BMI_act_cor[i]=cor(as.numeric(vip[as.character(coloc_nodes1p36$Node[i]),]),filt_pheno$BMI)
  coloc_nodes1p36$rs4654828_eQTL_Beta[i]=ifelse(coloc_nodes1p36$Node[i] %in% min_transE1$gene,
                                           min_transE1[min_transE1$gene==as.character(coloc_nodes1p36$Node[i]),"beta"],
                                           0)
  coloc_nodes1p36$rs4654828_eQTL_logP[i]=ifelse(coloc_nodes1p36$Node[i] %in% min_transE1$gene,
                                           -log10(min_transE1[min_transE1$gene==as.character(coloc_nodes1p36$Node[i]),"pvalue"]),
                                           0)
  coloc_nodes1p36$rs4654828_aQTL_Beta[i]=ifelse(coloc_nodes1p36$Node[i] %in% min_transA1$gene,
                                           min_transA1[min_transA1$gene==as.character(coloc_nodes1p36$Node[i]),"beta"],
                                           0)
  coloc_nodes1p36$rs4654828_aQTL_logP[i]=ifelse(coloc_nodes1p36$Node[i] %in% min_transA1$gene,
                                           -log10(min_transA1[min_transA1$gene==as.character(coloc_nodes1p36$Node[i]),"pvalue"]),
                                           0)
}

# Since EPHB2 is the only cis gene here, I'll just deal with it manually
inter_nodes1p36[inter_nodes1p36$Node=="EPHB2","rs4654828_eQTL_Beta"]=min_cisE1[min_cisE1$gene=="EPHB2","beta"]
inter_nodes1p36[inter_nodes1p36$Node=="EPHB2","rs4654828_eQTL_logP"]=-log10(min_cisE1[min_cisE1$gene=="EPHB2","pvalue"])
inter_nodes1p36[inter_nodes1p36$Node=="EPHB2","rs4654828_aQTL_Beta"]=min_cisA1[min_cisA1$gene=="EPHB2","beta"]
inter_nodes1p36[inter_nodes1p36$Node=="EPHB2","rs4654828_aQTL_logP"]=-log10(min_cisA1[min_cisA1$gene=="EPHB2","pvalue"])
coloc_nodes1p36[coloc_nodes1p36$Node=="EPHB2","rs4654828_eQTL_Beta"]=min_cisE1[min_cisE1$gene=="EPHB2","beta"]
coloc_nodes1p36[coloc_nodes1p36$Node=="EPHB2","rs4654828_eQTL_logP"]=-log10(min_cisE1[min_cisE1$gene=="EPHB2","pvalue"])
coloc_nodes1p36[coloc_nodes1p36$Node=="EPHB2","rs4654828_aQTL_Beta"]=min_cisA1[min_cisA1$gene=="EPHB2","beta"]
coloc_nodes1p36[coloc_nodes1p36$Node=="EPHB2","rs4654828_aQTL_logP"]=-log10(min_cisA1[min_cisA1$gene=="EPHB2","pvalue"])

# I think it may be more convenient to merge the networks into one and then just change which attributes I visualize in Cytoscape
# Start with 2 temporary columns concatinating the regulator-target and target-regulator for easier matching.
interactome1p36$temp1=paste(interactome1p36$Regulator,interactome1p36$Target)
interactome1p36$temp2=paste(interactome1p36$Target,interactome1p36$Regulator)
colocNet1p36$temp1=paste(colocNet1p36$trait1,colocNet1p36$traits)
colocNet1p36$temp2=paste(colocNet1p36$traits,colocNet1p36$trait1)

# Then grab colocalization data for gene pairs in interactome
temp=as.data.frame(matrix(nrow = dim(interactome1p36)[1],ncol = 4))
for(i in 1:dim(interactome1p36)[1]){
  temp[i,1:4]=colocNet1p36[ifelse(is.na(match(interactome1p36$temp1[i],colocNet1p36$temp1)),
                                  match(interactome1p36$temp1[i],colocNet1p36$temp2),
                                  match(interactome1p36$temp1[i],colocNet1p36$temp1)),3:6]
}
temp[is.na(temp)]=0

# Then combine with the BMI colocalizations
colnames(temp)=colnames(colocNet1p36)[3:6]
temp=rbind(temp,colocNet1p36[colocNet1p36$trait1=="BMI",3:6])

# Then add rows for BMI-Gene connections with 0 for MoA and likelihood
full1p36=interactome1p36[,1:4]
temp2=colocNet1p36[colocNet1p36$trait1=="BMI",1:4]
colnames(temp2)=colnames(interactome1p36)[1:4]
temp2[,3:4]=0
full1p36=rbind(full1p36,temp2)

# Finally, combine the colocalization columns with the interactome columns
full1p36=cbind(full1p36,temp)

# The nodes data also needs to be combined and duplicate rows removed
full1p36_nodes=rbind(inter_nodes1p36,coloc_nodes1p36)
full1p36_nodes=full1p36_nodes[!duplicated(full1p36_nodes$Node),]

# Write networks and node data to file for Cytoscape visualizations
write.table(interactome1p36,"Chr1p36_EPHB2_and_BMI_MRs_interactome.txt",sep = "\t",quote = F,row.names = F)
write.table(inter_nodes1p36,"Chr1p36_EPHB2_and_BMI_MRs_interactome_node_info.txt",sep = "\t",quote = F,row.names = F)
write.table(colocNet1p36,"Chr1p36_BMI_EPHB2_and_BMI_MRs_pairwise_colocalization_network.txt",sep = "\t",quote = F,row.names = F)
write.table(coloc_nodes1p36,"Chr1p36_BMI_EPHB2_and_BMI_MRs_pairwise_colocalization_network_node_info.txt",sep = "\t",quote = F,row.names = F)
write.table(full1p36,"Chr1p36_EPHB2_and_BMI_MRs_interactome_and_pairwise_colocalization.txt",sep = "\t",quote = F,row.names = F)
write.table(full1p36_nodes,"Chr1p36_EPHB2_and_BMI_MRs_interactome_and_pairwise_colocalization_node_info.txt",sep = "\t",quote = F,row.names = F)

# 7q32
# Grab interactions between LINC-PINT, KLF14 and MRs in adipose interactome
linc_bmi_mrs=bmi_MRregs[bmi_MRregs$Regulator=="LINC-PINT",]
linc_homair_mrs=homair_MRregs[homair_MRregs$Regulator=="LINC-PINT",]
linc_hdl_mrs=hdl_MRregs[hdl_MRregs$Regulator=="LINC-PINT",]
linc_triG_mrs=triG_MRregs[triG_MRregs$Regulator=="LINC-PINT",]
bmi_mrsLINC=interactome[(interactome$Regulator %in% bmi_mrs[,1]) & (interactome$Target=="LINC-PINT"),]
homair_mrsLINC=interactome[(interactome$Regulator %in% homair_mrs[,1]) & (interactome$Target=="LINC-PINT"),]
hdl_mrsLINC=interactome[(interactome$Regulator %in% hdl_mrs[,1]) & (interactome$Target=="LINC-PINT"),]
triG_mrsLINC=interactome[(interactome$Regulator %in% triG_mrs[,1]) & (interactome$Target=="LINC-PINT"),]

klf14_bmi_mrs=bmi_MRregs[bmi_MRregs$Regulator=="KLF14",]
klf14_homair_mrs=homair_MRregs[homair_MRregs$Regulator=="KLF14",]
klf14_hdl_mrs=hdl_MRregs[hdl_MRregs$Regulator=="KLF14",]
klf14_triG_mrs=triG_MRregs[triG_MRregs$Regulator=="KLF14",]
bmi_mrsKLF14=interactome[(interactome$Regulator %in% bmi_mrs[,1]) & (interactome$Target=="KLF14"),]
homair_mrsKLF14=interactome[(interactome$Regulator %in% homair_mrs[,1]) & (interactome$Target=="KLF14"),]
hdl_mrsKLF14=interactome[(interactome$Regulator %in% hdl_mrs[,1]) & (interactome$Target=="KLF14"),]
triG_mrsKLF14=interactome[(interactome$Regulator %in% triG_mrs[,1]) & (interactome$Target=="KLF14"),]

bmi_interactome7q32=rbind(bmi_MRMR,linc_bmi_mrs,bmi_mrsLINC,klf14_bmi_mrs,bmi_mrsKLF14)
homair_interactome7q32=rbind(homair_MRMR,linc_homair_mrs,homair_mrsLINC,klf14_homair_mrs,homair_mrsKLF14)
hdl_interactome7q32=rbind(hdl_MRMR,linc_hdl_mrs,hdl_mrsLINC,klf14_hdl_mrs,hdl_mrsKLF14)
triG_interactome7q32=rbind(triG_MRMR,linc_triG_mrs,triG_mrsLINC,klf14_triG_mrs,triG_mrsKLF14)

# Grab pairwise colocalizations with PP>0.5 between each GWAS and QTLs. I did not run pairwise colocalization analyses for LINC-PINT or KLF14 yet.
bmi_pairColoc7q32=bmi_pairColoc[bmi_pairColoc$locus=="7q32" & bmi_pairColoc$posterior_prob>0.5,c(2,3,5)]
bmi_pairColoc7q32$traits=gsub("BMI, ","",bmi_pairColoc7q32$traits)
bmi_pairColoc7q32=cbind("trait1"=rep("BMI",dim(bmi_pairColoc7q32)[1]),bmi_pairColoc7q32)
bmi_e_pairColoc1=bmi_pairColoc7q32[grepl("-e_",bmi_pairColoc7q32$traits),]
bmi_a_pairColoc1=bmi_pairColoc7q32[grepl("-a_",bmi_pairColoc7q32$traits),]
bmi_e_pairColoc1$traits=gsub(".*_","",bmi_e_pairColoc1$traits)
bmi_a_pairColoc1$traits=gsub(".*_","",bmi_a_pairColoc1$traits)
bmi_netPair7q32=rbind(bmi_e_pairColoc1,bmi_a_pairColoc1)
bmi_netPair7q32=bmi_netPair7q32[!duplicated(bmi_netPair7q32$traits),-c(3,4)]
bmi_netPair7q32$eQTL_PP=bmi_e_pairColoc1[match(bmi_netPair7q32$traits,bmi_e_pairColoc1$traits),3]
bmi_netPair7q32$eQTL_SNP=bmi_e_pairColoc1[match(bmi_netPair7q32$traits,bmi_e_pairColoc1$traits),4]
bmi_netPair7q32$aQTL_PP=bmi_a_pairColoc1[match(bmi_netPair7q32$traits,bmi_a_pairColoc1$traits),3]
bmi_netPair7q32$aQTL_SNP=bmi_a_pairColoc1[match(bmi_netPair7q32$traits,bmi_a_pairColoc1$traits),4]
bmi_netPair7q32[is.na(bmi_netPair7q32)]=0

t2d_pairColoc7q32=t2d_pairColoc[t2d_pairColoc$locus=="7q32" & t2d_pairColoc$posterior_prob>0.5,c(2,3,5)]
t2d_pairColoc7q32$traits=gsub("T2D, ","",t2d_pairColoc7q32$traits)
t2d_pairColoc7q32=cbind("trait1"=rep("T2D",dim(t2d_pairColoc7q32)[1]),t2d_pairColoc7q32)
t2d_e_pairColoc1=t2d_pairColoc7q32[grepl("-e_",t2d_pairColoc7q32$traits),]
t2d_a_pairColoc1=t2d_pairColoc7q32[grepl("-a_",t2d_pairColoc7q32$traits),]
t2d_e_pairColoc1$traits=gsub(".*_","",t2d_e_pairColoc1$traits)
t2d_a_pairColoc1$traits=gsub(".*_","",t2d_a_pairColoc1$traits)
t2d_netPair7q32=rbind(t2d_e_pairColoc1,t2d_a_pairColoc1)
t2d_netPair7q32=t2d_netPair7q32[!duplicated(t2d_netPair7q32$traits),-c(3,4)]
t2d_netPair7q32$eQTL_PP=t2d_e_pairColoc1[match(t2d_netPair7q32$traits,t2d_e_pairColoc1$traits),3]
t2d_netPair7q32$eQTL_SNP=t2d_e_pairColoc1[match(t2d_netPair7q32$traits,t2d_e_pairColoc1$traits),4]
t2d_netPair7q32$aQTL_PP=t2d_a_pairColoc1[match(t2d_netPair7q32$traits,t2d_a_pairColoc1$traits),3]
t2d_netPair7q32$aQTL_SNP=t2d_a_pairColoc1[match(t2d_netPair7q32$traits,t2d_a_pairColoc1$traits),4]
t2d_netPair7q32[is.na(t2d_netPair7q32)]=0

hdl_pairColoc7q32=hdl_pairColoc[hdl_pairColoc$locus=="7q32" & hdl_pairColoc$posterior_prob>0.5,c(2,3,5)]
hdl_pairColoc7q32$traits=gsub("HDL, ","",hdl_pairColoc7q32$traits)
hdl_pairColoc7q32=cbind("trait1"=rep("HDL",dim(hdl_pairColoc7q32)[1]),hdl_pairColoc7q32)
hdl_e_pairColoc1=hdl_pairColoc7q32[grepl("-e_",hdl_pairColoc7q32$traits),]
hdl_a_pairColoc1=hdl_pairColoc7q32[grepl("-a_",hdl_pairColoc7q32$traits),]
hdl_e_pairColoc1$traits=gsub(".*_","",hdl_e_pairColoc1$traits)
hdl_a_pairColoc1$traits=gsub(".*_","",hdl_a_pairColoc1$traits)
hdl_netPair7q32=rbind(hdl_e_pairColoc1,hdl_a_pairColoc1)
hdl_netPair7q32=hdl_netPair7q32[!duplicated(hdl_netPair7q32$traits),-c(3,4)]
hdl_netPair7q32$eQTL_PP=hdl_e_pairColoc1[match(hdl_netPair7q32$traits,hdl_e_pairColoc1$traits),3]
hdl_netPair7q32$eQTL_SNP=hdl_e_pairColoc1[match(hdl_netPair7q32$traits,hdl_e_pairColoc1$traits),4]
hdl_netPair7q32$aQTL_PP=hdl_a_pairColoc1[match(hdl_netPair7q32$traits,hdl_a_pairColoc1$traits),3]
hdl_netPair7q32$aQTL_SNP=hdl_a_pairColoc1[match(hdl_netPair7q32$traits,hdl_a_pairColoc1$traits),4]
hdl_netPair7q32[is.na(hdl_netPair7q32)]=0

triG_pairColoc7q32=triG_pairColoc[triG_pairColoc$locus=="7q32" & triG_pairColoc$posterior_prob>0.5,c(2,3,5)]
triG_pairColoc7q32$traits=gsub("TriG, ","",triG_pairColoc7q32$traits)
triG_pairColoc7q32=cbind("trait1"=rep("TriG",dim(triG_pairColoc7q32)[1]),triG_pairColoc7q32)
triG_e_pairColoc1=triG_pairColoc7q32[grepl("-e_",triG_pairColoc7q32$traits),]
triG_a_pairColoc1=triG_pairColoc7q32[grepl("-a_",triG_pairColoc7q32$traits),]
triG_e_pairColoc1$traits=gsub(".*_","",triG_e_pairColoc1$traits)
triG_a_pairColoc1$traits=gsub(".*_","",triG_a_pairColoc1$traits)
triG_netPair7q32=rbind(triG_e_pairColoc1,triG_a_pairColoc1)
triG_netPair7q32=triG_netPair7q32[!duplicated(triG_netPair7q32$traits),-c(3,4)]
triG_netPair7q32$eQTL_PP=triG_e_pairColoc1[match(triG_netPair7q32$traits,triG_e_pairColoc1$traits),3]
triG_netPair7q32$eQTL_SNP=triG_e_pairColoc1[match(triG_netPair7q32$traits,triG_e_pairColoc1$traits),4]
triG_netPair7q32$aQTL_PP=triG_a_pairColoc1[match(triG_netPair7q32$traits,triG_a_pairColoc1$traits),3]
triG_netPair7q32$aQTL_SNP=triG_a_pairColoc1[match(triG_netPair7q32$traits,triG_a_pairColoc1$traits),4]
triG_netPair7q32[is.na(triG_netPair7q32)]=0

# Grab the -log10(Pmin) and betas for the eQTLs and aQTLs among GWAS significant SNPs at the 7q32 locus
bmi_cisE1=filt_cis_eqtl[filt_cis_eqtl$chr==7,]
bmi_cisE1=bmi_cisE1[bmi_cisE1$snps %in% sig_bmi$SNP,]
bmi_cisE1=bmi_cisE1[order(bmi_cisE1$pvalue),]
min_bmi_cisE1=bmi_cisE1[!duplicated(bmi_cisE1$gene),]

bmi_cisA1=filt_cis_aqtl[filt_cis_aqtl$chr==7,]
bmi_cisA1=bmi_cisA1[bmi_cisA1$snps %in% sig_bmi$SNP,]
bmi_cisA1=bmi_cisA1[order(bmi_cisA1$pvalue),]
min_bmi_cisA1=bmi_cisA1[!duplicated(bmi_cisA1$gene),]

bmi_transE1=filt_trans_bmi_eqtl[filt_trans_bmi_eqtl$chr==7,]
bmi_transE1=bmi_transE1[bmi_transE1$snps %in% sig_bmi$SNP,]
bmi_transE1=bmi_transE1[order(bmi_transE1$pvalue),]
min_bmi_transE1=bmi_transE1[!duplicated(bmi_transE1$gene),]

bmi_transA1=filt_trans_bmi_aqtl[filt_trans_bmi_aqtl$chr==7,]
bmi_transA1=bmi_transA1[bmi_transA1$snps %in% sig_bmi$SNP,]
bmi_transA1=bmi_transA1[order(bmi_transA1$pvalue),]
min_bmi_transA1=bmi_transA1[!duplicated(bmi_transA1$gene),]

t2d_cisE1=filt_cis_eqtl[filt_cis_eqtl$chr==7,]
t2d_cisE1=t2d_cisE1[t2d_cisE1$snps %in% sig_t2d$rsID,]
t2d_cisE1=t2d_cisE1[order(t2d_cisE1$pvalue),]
min_t2d_cisE1=t2d_cisE1[!duplicated(t2d_cisE1$gene),]

t2d_cisA1=filt_cis_aqtl[filt_cis_aqtl$chr==7,]
t2d_cisA1=t2d_cisA1[t2d_cisA1$snps %in% sig_t2d$rsID,]
t2d_cisA1=t2d_cisA1[order(t2d_cisA1$pvalue),]
min_t2d_cisA1=t2d_cisA1[!duplicated(t2d_cisA1$gene),]

t2d_transE1=filt_trans_t2d_eqtl[filt_trans_t2d_eqtl$chr==7,]
t2d_transE1=t2d_transE1[t2d_transE1$snps %in% sig_t2d$rsID,]
t2d_transE1=t2d_transE1[order(t2d_transE1$pvalue),]
min_t2d_transE1=t2d_transE1[!duplicated(t2d_transE1$gene),]

t2d_transA1=filt_trans_t2d_aqtl[filt_trans_t2d_aqtl$chr==7,]
t2d_transA1=t2d_transA1[t2d_transA1$snps %in% sig_t2d$rsID,]
t2d_transA1=t2d_transA1[order(t2d_transA1$pvalue),]
min_t2d_transA1=t2d_transA1[!duplicated(t2d_transA1$gene),]

hdl_cisE1=filt_cis_eqtl[filt_cis_eqtl$chr==7,]
hdl_cisE1=hdl_cisE1[hdl_cisE1$snps %in% sig_hdl$SNP,]
hdl_cisE1=hdl_cisE1[order(hdl_cisE1$pvalue),]
min_hdl_cisE1=hdl_cisE1[!duplicated(hdl_cisE1$gene),]

hdl_cisA1=filt_cis_aqtl[filt_cis_aqtl$chr==7,]
hdl_cisA1=hdl_cisA1[hdl_cisA1$snps %in% sig_hdl$SNP,]
hdl_cisA1=hdl_cisA1[order(hdl_cisA1$pvalue),]
min_hdl_cisA1=hdl_cisA1[!duplicated(hdl_cisA1$gene),]

hdl_transE1=filt_trans_hdl_eqtl[filt_trans_hdl_eqtl$chr==7,]
hdl_transE1=hdl_transE1[hdl_transE1$snps %in% sig_hdl$SNP,]
hdl_transE1=hdl_transE1[order(hdl_transE1$pvalue),]
min_hdl_transE1=hdl_transE1[!duplicated(hdl_transE1$gene),]

hdl_transA1=filt_trans_hdl_aqtl[filt_trans_hdl_aqtl$chr==7,]
hdl_transA1=hdl_transA1[hdl_transA1$snps %in% sig_hdl$SNP,]
hdl_transA1=hdl_transA1[order(hdl_transA1$pvalue),]
min_hdl_transA1=hdl_transA1[!duplicated(hdl_transA1$gene),]

triG_cisE1=filt_cis_eqtl[filt_cis_eqtl$chr==7,]
triG_cisE1=triG_cisE1[triG_cisE1$snps %in% sig_triG$SNP,]
triG_cisE1=triG_cisE1[order(triG_cisE1$pvalue),]
min_triG_cisE1=triG_cisE1[!duplicated(triG_cisE1$gene),]

triG_cisA1=filt_cis_aqtl[filt_cis_aqtl$chr==7,]
triG_cisA1=triG_cisA1[triG_cisA1$snps %in% sig_triG$SNP,]
triG_cisA1=triG_cisA1[order(triG_cisA1$pvalue),]
min_triG_cisA1=triG_cisA1[!duplicated(triG_cisA1$gene),]

triG_transE1=filt_trans_triG_eqtl[filt_trans_triG_eqtl$chr==7,]
triG_transE1=triG_transE1[triG_transE1$snps %in% sig_triG$SNP,]
triG_transE1=triG_transE1[order(triG_transE1$pvalue),]
min_triG_transE1=triG_transE1[!duplicated(triG_transE1$gene),]

triG_transA1=filt_trans_triG_aqtl[filt_trans_triG_aqtl$chr==7,]
triG_transA1=triG_transA1[triG_transA1$snps %in% sig_triG$SNP,]
triG_transA1=triG_transA1[order(triG_transA1$pvalue),]
min_triG_transA1=triG_transA1[!duplicated(triG_transA1$gene),]

# Make node tables for 7q32 networks for each GWAS. Note that some GWAS (T2D and TriG) failed to have their MRs connect at all with LINC-PINT or KLF14, 
# so I manually added those to the node lists when needed.

# BMI
bmi_inter_nodes7q32=data.frame("Node"=as.character(unique(bmi_interactome7q32$Regulator)),"BMI_exp_cor"=rep(0,length(unique(bmi_interactome7q32$Regulator))),
                               "BMI_act_cor"=rep(0,length(unique(bmi_interactome7q32$Regulator))),"Best_eQTL_Beta"=rep(0,length(unique(bmi_interactome7q32$Regulator))),
                               "Best_eQTL_logP"=rep(0,length(unique(bmi_interactome7q32$Regulator))),"Best_aQTL_Beta"=rep(0,length(unique(bmi_interactome7q32$Regulator))),
                               "Best_aQTL_logP"=rep(0,length(unique(bmi_interactome7q32$Regulator))))
bmi_coloc_nodes7q32=data.frame("Node"=as.character(unique(bmi_netPair7q32$traits)),"BMI_exp_cor"=rep(0,length(unique(bmi_netPair7q32$traits))),
                               "BMI_act_cor"=rep(0,length(unique(bmi_netPair7q32$traits))),"Best_eQTL_Beta"=rep(0,length(unique(bmi_netPair7q32$traits))),
                               "Best_eQTL_logP"=rep(0,length(unique(bmi_netPair7q32$traits))),"Best_aQTL_Beta"=rep(0,length(unique(bmi_netPair7q32$traits))),
                               "Best_aQTL_logP"=rep(0,length(unique(bmi_netPair7q32$traits))))

for(i in 1:dim(bmi_inter_nodes7q32)[1]){
  bmi_inter_nodes7q32$BMI_exp_cor[i]=cor(as.numeric(tpm[as.character(bmi_inter_nodes7q32$Node[i]),]),filt_pheno$BMI)
  bmi_inter_nodes7q32$BMI_act_cor[i]=cor(as.numeric(vip[as.character(bmi_inter_nodes7q32$Node[i]),]),filt_pheno$BMI)
  bmi_inter_nodes7q32$Best_eQTL_Beta[i]=ifelse(bmi_inter_nodes7q32$Node[i] %in% min_bmi_transE1$gene,
                                           min_bmi_transE1[min_bmi_transE1$gene==as.character(bmi_inter_nodes7q32$Node[i]),"beta"],
                                           0)
  bmi_inter_nodes7q32$Best_eQTL_logP[i]=ifelse(bmi_inter_nodes7q32$Node[i] %in% min_bmi_transE1$gene,
                                           -log10(min_bmi_transE1[min_bmi_transE1$gene==as.character(bmi_inter_nodes7q32$Node[i]),"pvalue"]),
                                           0)
  bmi_inter_nodes7q32$Best_aQTL_Beta[i]=ifelse(bmi_inter_nodes7q32$Node[i] %in% min_bmi_transA1$gene,
                                           min_bmi_transA1[min_bmi_transA1$gene==as.character(bmi_inter_nodes7q32$Node[i]),"beta"],
                                           0)
  bmi_inter_nodes7q32$Best_aQTL_logP[i]=ifelse(bmi_inter_nodes7q32$Node[i] %in% min_bmi_transA1$gene,
                                           -log10(min_bmi_transA1[min_bmi_transA1$gene==as.character(bmi_inter_nodes7q32$Node[i]),"pvalue"]),
                                           0)
}

for(i in 1:dim(bmi_coloc_nodes7q32)[1]){
  bmi_coloc_nodes7q32$BMI_exp_cor[i]=cor(as.numeric(tpm[as.character(bmi_coloc_nodes7q32$Node[i]),]),filt_pheno$BMI)
  bmi_coloc_nodes7q32$BMI_act_cor[i]=cor(as.numeric(vip[as.character(bmi_coloc_nodes7q32$Node[i]),]),filt_pheno$BMI)
  bmi_coloc_nodes7q32$Best_eQTL_Beta[i]=ifelse(bmi_coloc_nodes7q32$Node[i] %in% min_bmi_transE1$gene,
                                               min_bmi_transE1[min_bmi_transE1$gene==as.character(bmi_coloc_nodes7q32$Node[i]),"beta"],
                                               0)
  bmi_coloc_nodes7q32$Best_eQTL_logP[i]=ifelse(bmi_coloc_nodes7q32$Node[i] %in% min_bmi_transE1$gene,
                                               -log10(min_bmi_transE1[min_bmi_transE1$gene==as.character(bmi_coloc_nodes7q32$Node[i]),"pvalue"]),
                                               0)
  bmi_coloc_nodes7q32$Best_aQTL_Beta[i]=ifelse(bmi_coloc_nodes7q32$Node[i] %in% min_bmi_transA1$gene,
                                               min_bmi_transA1[min_bmi_transA1$gene==as.character(bmi_coloc_nodes7q32$Node[i]),"beta"],
                                               0)
  bmi_coloc_nodes7q32$Best_aQTL_logP[i]=ifelse(bmi_coloc_nodes7q32$Node[i] %in% min_bmi_transA1$gene,
                                               -log10(min_bmi_transA1[min_bmi_transA1$gene==as.character(bmi_coloc_nodes7q32$Node[i]),"pvalue"]),
                                               0)
}

# T2D
t2d_inter_nodes7q32=data.frame("Node"=c(as.character(unique(homair_interactome7q32$Regulator)),"LINC-PINT"),"HOMA.IR_exp_cor"=rep(0,length(unique(homair_interactome7q32$Regulator))+1),
                               "HOMA.IR_act_cor"=rep(0,length(unique(homair_interactome7q32$Regulator))+1),"Best_eQTL_Beta"=rep(0,length(unique(homair_interactome7q32$Regulator))+1),
                               "Best_eQTL_logP"=rep(0,length(unique(homair_interactome7q32$Regulator))+1),"Best_aQTL_Beta"=rep(0,length(unique(homair_interactome7q32$Regulator))+1),
                               "Best_aQTL_logP"=rep(0,length(unique(homair_interactome7q32$Regulator))+1))
t2d_coloc_nodes7q32=data.frame("Node"=as.character(unique(t2d_netPair7q32$traits)),"HOMA.IR_exp_cor"=rep(0,length(unique(t2d_netPair7q32$traits))),
                               "HOMA.IR_act_cor"=rep(0,length(unique(t2d_netPair7q32$traits))),"Best_eQTL_Beta"=rep(0,length(unique(t2d_netPair7q32$traits))),
                               "Best_eQTL_logP"=rep(0,length(unique(t2d_netPair7q32$traits))),"Best_aQTL_Beta"=rep(0,length(unique(t2d_netPair7q32$traits))),
                               "Best_aQTL_logP"=rep(0,length(unique(t2d_netPair7q32$traits))))

for(i in 1:dim(t2d_inter_nodes7q32)[1]){
  noNA_samples=rownames(filt_pheno)[!is.na(filt_pheno$HOMA.IR)]
  t2d_inter_nodes7q32$HOMA.IR_exp_cor[i]=cor(as.numeric(tpm[as.character(t2d_inter_nodes7q32$Node[i]),noNA_samples]),filt_pheno[noNA_samples,"HOMA.IR"])
  t2d_inter_nodes7q32$HOMA.IR_act_cor[i]=cor(as.numeric(vip[as.character(t2d_inter_nodes7q32$Node[i]),noNA_samples]),filt_pheno[noNA_samples,"HOMA.IR"])
  t2d_inter_nodes7q32$Best_eQTL_Beta[i]=ifelse(t2d_inter_nodes7q32$Node[i] %in% min_t2d_transE1$gene,
                                               min_t2d_transE1[min_t2d_transE1$gene==as.character(t2d_inter_nodes7q32$Node[i]),"beta"],
                                               0)
  t2d_inter_nodes7q32$Best_eQTL_logP[i]=ifelse(t2d_inter_nodes7q32$Node[i] %in% min_t2d_transE1$gene,
                                               -log10(min_t2d_transE1[min_t2d_transE1$gene==as.character(t2d_inter_nodes7q32$Node[i]),"pvalue"]),
                                               0)
  t2d_inter_nodes7q32$Best_aQTL_Beta[i]=ifelse(t2d_inter_nodes7q32$Node[i] %in% min_t2d_transA1$gene,
                                               min_t2d_transA1[min_t2d_transA1$gene==as.character(t2d_inter_nodes7q32$Node[i]),"beta"],
                                               0)
  t2d_inter_nodes7q32$Best_aQTL_logP[i]=ifelse(t2d_inter_nodes7q32$Node[i] %in% min_t2d_transA1$gene,
                                               -log10(min_t2d_transA1[min_t2d_transA1$gene==as.character(t2d_inter_nodes7q32$Node[i]),"pvalue"]),
                                               0)
}

for(i in 1:dim(t2d_coloc_nodes7q32)[1]){
  noNA_samples=rownames(filt_pheno)[!is.na(filt_pheno$HOMA.IR)]
  t2d_coloc_nodes7q32$HOMA.IR_exp_cor[i]=cor(as.numeric(tpm[as.character(t2d_coloc_nodes7q32$Node[i]),noNA_samples]),filt_pheno[noNA_samples,"HOMA.IR"])
  t2d_coloc_nodes7q32$HOMA.IR_act_cor[i]=cor(as.numeric(vip[as.character(t2d_coloc_nodes7q32$Node[i]),noNA_samples]),filt_pheno[noNA_samples,"HOMA.IR"])
  t2d_coloc_nodes7q32$Best_eQTL_Beta[i]=ifelse(t2d_coloc_nodes7q32$Node[i] %in% min_t2d_transE1$gene,
                                               min_t2d_transE1[min_t2d_transE1$gene==as.character(t2d_coloc_nodes7q32$Node[i]),"beta"],
                                               0)
  t2d_coloc_nodes7q32$Best_eQTL_logP[i]=ifelse(t2d_coloc_nodes7q32$Node[i] %in% min_t2d_transE1$gene,
                                               -log10(min_t2d_transE1[min_t2d_transE1$gene==as.character(t2d_coloc_nodes7q32$Node[i]),"pvalue"]),
                                               0)
  t2d_coloc_nodes7q32$Best_aQTL_Beta[i]=ifelse(t2d_coloc_nodes7q32$Node[i] %in% min_t2d_transA1$gene,
                                               min_t2d_transA1[min_t2d_transA1$gene==as.character(t2d_coloc_nodes7q32$Node[i]),"beta"],
                                               0)
  t2d_coloc_nodes7q32$Best_aQTL_logP[i]=ifelse(t2d_coloc_nodes7q32$Node[i] %in% min_t2d_transA1$gene,
                                               -log10(min_t2d_transA1[min_t2d_transA1$gene==as.character(t2d_coloc_nodes7q32$Node[i]),"pvalue"]),
                                               0)
}

# HDL
hdl_inter_nodes7q32=data.frame("Node"=as.character(unique(hdl_interactome7q32$Regulator)),"HDL_exp_cor"=rep(0,length(unique(hdl_interactome7q32$Regulator))),
                               "HDL_act_cor"=rep(0,length(unique(hdl_interactome7q32$Regulator))),"Best_eQTL_Beta"=rep(0,length(unique(hdl_interactome7q32$Regulator))),
                               "Best_eQTL_logP"=rep(0,length(unique(hdl_interactome7q32$Regulator))),"Best_aQTL_Beta"=rep(0,length(unique(hdl_interactome7q32$Regulator))),
                               "Best_aQTL_logP"=rep(0,length(unique(hdl_interactome7q32$Regulator))))
hdl_coloc_nodes7q32=data.frame("Node"=as.character(unique(hdl_netPair7q32$traits)),"HDL_exp_cor"=rep(0,length(unique(hdl_netPair7q32$traits))),
                               "HDL_act_cor"=rep(0,length(unique(hdl_netPair7q32$traits))),"Best_eQTL_Beta"=rep(0,length(unique(hdl_netPair7q32$traits))),
                               "Best_eQTL_logP"=rep(0,length(unique(hdl_netPair7q32$traits))),"Best_aQTL_Beta"=rep(0,length(unique(hdl_netPair7q32$traits))),
                               "Best_aQTL_logP"=rep(0,length(unique(hdl_netPair7q32$traits))))

for(i in 1:dim(hdl_inter_nodes7q32)[1]){
  noNA_samples=rownames(filt_pheno)[!is.na(filt_pheno$HDLcholesterol)]
  hdl_inter_nodes7q32$HDL_exp_cor[i]=cor(as.numeric(tpm[as.character(hdl_inter_nodes7q32$Node[i]),noNA_samples]),filt_pheno[noNA_samples,"HDLcholesterol"])
  hdl_inter_nodes7q32$HDL_act_cor[i]=cor(as.numeric(vip[as.character(hdl_inter_nodes7q32$Node[i]),noNA_samples]),filt_pheno[noNA_samples,"HDLcholesterol"])
  hdl_inter_nodes7q32$Best_eQTL_Beta[i]=ifelse(hdl_inter_nodes7q32$Node[i] %in% min_hdl_transE1$gene,
                                               min_hdl_transE1[min_hdl_transE1$gene==as.character(hdl_inter_nodes7q32$Node[i]),"beta"],
                                               0)
  hdl_inter_nodes7q32$Best_eQTL_logP[i]=ifelse(hdl_inter_nodes7q32$Node[i] %in% min_hdl_transE1$gene,
                                               -log10(min_hdl_transE1[min_hdl_transE1$gene==as.character(hdl_inter_nodes7q32$Node[i]),"pvalue"]),
                                               0)
  hdl_inter_nodes7q32$Best_aQTL_Beta[i]=ifelse(hdl_inter_nodes7q32$Node[i] %in% min_hdl_transA1$gene,
                                               min_hdl_transA1[min_hdl_transA1$gene==as.character(hdl_inter_nodes7q32$Node[i]),"beta"],
                                               0)
  hdl_inter_nodes7q32$Best_aQTL_logP[i]=ifelse(hdl_inter_nodes7q32$Node[i] %in% min_hdl_transA1$gene,
                                               -log10(min_hdl_transA1[min_hdl_transA1$gene==as.character(hdl_inter_nodes7q32$Node[i]),"pvalue"]),
                                               0)
}

for(i in 1:dim(hdl_coloc_nodes7q32)[1]){
  noNA_samples=rownames(filt_pheno)[!is.na(filt_pheno$HDLcholesterol)]
  hdl_coloc_nodes7q32$HDL_exp_cor[i]=cor(as.numeric(tpm[as.character(hdl_coloc_nodes7q32$Node[i]),noNA_samples]),filt_pheno[noNA_samples,"HDLcholesterol"])
  hdl_coloc_nodes7q32$HDL_act_cor[i]=cor(as.numeric(vip[as.character(hdl_coloc_nodes7q32$Node[i]),noNA_samples]),filt_pheno[noNA_samples,"HDLcholesterol"])
  hdl_coloc_nodes7q32$Best_eQTL_Beta[i]=ifelse(hdl_coloc_nodes7q32$Node[i] %in% min_hdl_transE1$gene,
                                               min_hdl_transE1[min_hdl_transE1$gene==as.character(hdl_coloc_nodes7q32$Node[i]),"beta"],
                                               0)
  hdl_coloc_nodes7q32$Best_eQTL_logP[i]=ifelse(hdl_coloc_nodes7q32$Node[i] %in% min_hdl_transE1$gene,
                                               -log10(min_hdl_transE1[min_hdl_transE1$gene==as.character(hdl_coloc_nodes7q32$Node[i]),"pvalue"]),
                                               0)
  hdl_coloc_nodes7q32$Best_aQTL_Beta[i]=ifelse(hdl_coloc_nodes7q32$Node[i] %in% min_hdl_transA1$gene,
                                               min_hdl_transA1[min_hdl_transA1$gene==as.character(hdl_coloc_nodes7q32$Node[i]),"beta"],
                                               0)
  hdl_coloc_nodes7q32$Best_aQTL_logP[i]=ifelse(hdl_coloc_nodes7q32$Node[i] %in% min_hdl_transA1$gene,
                                               -log10(min_hdl_transA1[min_hdl_transA1$gene==as.character(hdl_coloc_nodes7q32$Node[i]),"pvalue"]),
                                               0)
}

# TriG
triG_inter_nodes7q32=data.frame("Node"=c(as.character(unique(triG_interactome7q32$Regulator)),"LINC-PINT","KLF14"),"TriG_exp_cor"=rep(0,length(unique(triG_interactome7q32$Regulator))+2),
                                "TriG_act_cor"=rep(0,length(unique(triG_interactome7q32$Regulator))+2),"Best_eQTL_Beta"=rep(0,length(unique(triG_interactome7q32$Regulator))+2),
                                "Best_eQTL_logP"=rep(0,length(unique(triG_interactome7q32$Regulator))+2),"Best_aQTL_Beta"=rep(0,length(unique(triG_interactome7q32$Regulator))+2),
                                "Best_aQTL_logP"=rep(0,length(unique(triG_interactome7q32$Regulator))+2))
triG_coloc_nodes7q32=data.frame("Node"=as.character(unique(triG_netPair7q32$traits)),"TriG_exp_cor"=rep(0,length(unique(triG_netPair7q32$traits))),
                                "TriG_act_cor"=rep(0,length(unique(triG_netPair7q32$traits))),"Best_eQTL_Beta"=rep(0,length(unique(triG_netPair7q32$traits))),
                                "Best_eQTL_logP"=rep(0,length(unique(triG_netPair7q32$traits))),"Best_aQTL_Beta"=rep(0,length(unique(triG_netPair7q32$traits))),
                                "Best_aQTL_logP"=rep(0,length(unique(triG_netPair7q32$traits))))

for(i in 1:dim(triG_inter_nodes7q32)[1]){
  noNA_samples=rownames(filt_pheno)[!is.na(filt_pheno$TotalTriglycerides)]
  triG_inter_nodes7q32$TriG_exp_cor[i]=cor(as.numeric(tpm[as.character(triG_inter_nodes7q32$Node[i]),noNA_samples]),filt_pheno[noNA_samples,"TotalTriglycerides"])
  triG_inter_nodes7q32$TriG_act_cor[i]=cor(as.numeric(vip[as.character(triG_inter_nodes7q32$Node[i]),noNA_samples]),filt_pheno[noNA_samples,"TotalTriglycerides"])
  triG_inter_nodes7q32$Best_eQTL_Beta[i]=ifelse(triG_inter_nodes7q32$Node[i] %in% min_triG_transE1$gene,
                                               min_triG_transE1[min_triG_transE1$gene==as.character(triG_inter_nodes7q32$Node[i]),"beta"],
                                               0)
  triG_inter_nodes7q32$Best_eQTL_logP[i]=ifelse(triG_inter_nodes7q32$Node[i] %in% min_triG_transE1$gene,
                                               -log10(min_triG_transE1[min_triG_transE1$gene==as.character(triG_inter_nodes7q32$Node[i]),"pvalue"]),
                                               0)
  triG_inter_nodes7q32$Best_aQTL_Beta[i]=ifelse(triG_inter_nodes7q32$Node[i] %in% min_triG_transA1$gene,
                                               min_triG_transA1[min_triG_transA1$gene==as.character(triG_inter_nodes7q32$Node[i]),"beta"],
                                               0)
  triG_inter_nodes7q32$Best_aQTL_logP[i]=ifelse(triG_inter_nodes7q32$Node[i] %in% min_triG_transA1$gene,
                                               -log10(min_triG_transA1[min_triG_transA1$gene==as.character(triG_inter_nodes7q32$Node[i]),"pvalue"]),
                                               0)
}

for(i in 1:dim(triG_coloc_nodes7q32)[1]){
  noNA_samples=rownames(filt_pheno)[!is.na(filt_pheno$TotalTriglycerides)]
  triG_coloc_nodes7q32$TriG_exp_cor[i]=cor(as.numeric(tpm[as.character(triG_coloc_nodes7q32$Node[i]),noNA_samples]),filt_pheno[noNA_samples,"TotalTriglycerides"])
  triG_coloc_nodes7q32$TriG_act_cor[i]=cor(as.numeric(vip[as.character(triG_coloc_nodes7q32$Node[i]),noNA_samples]),filt_pheno[noNA_samples,"TotalTriglycerides"])
  triG_coloc_nodes7q32$Best_eQTL_Beta[i]=ifelse(triG_coloc_nodes7q32$Node[i] %in% min_triG_transE1$gene,
                                               min_triG_transE1[min_triG_transE1$gene==as.character(triG_coloc_nodes7q32$Node[i]),"beta"],
                                               0)
  triG_coloc_nodes7q32$Best_eQTL_logP[i]=ifelse(triG_coloc_nodes7q32$Node[i] %in% min_triG_transE1$gene,
                                               -log10(min_triG_transE1[min_triG_transE1$gene==as.character(triG_coloc_nodes7q32$Node[i]),"pvalue"]),
                                               0)
  triG_coloc_nodes7q32$Best_aQTL_Beta[i]=ifelse(triG_coloc_nodes7q32$Node[i] %in% min_triG_transA1$gene,
                                               min_triG_transA1[min_triG_transA1$gene==as.character(triG_coloc_nodes7q32$Node[i]),"beta"],
                                               0)
  triG_coloc_nodes7q32$Best_aQTL_logP[i]=ifelse(triG_coloc_nodes7q32$Node[i] %in% min_triG_transA1$gene,
                                               -log10(min_triG_transA1[min_triG_transA1$gene==as.character(triG_coloc_nodes7q32$Node[i]),"pvalue"]),
                                               0)
}

# I think it may be more convenient to merge the networks into one and then just change which attributes I visualize in Cytoscape
# Start with 2 temporary columns concatinating the regulator-target and target-regulator for easier matching.

# BMI
bmi_interactome7q32$temp1=paste(bmi_interactome7q32$Regulator,bmi_interactome7q32$Target)
bmi_interactome7q32$temp2=paste(bmi_interactome7q32$Target,bmi_interactome7q32$Regulator)
bmi_netPair7q32$temp1=paste(bmi_netPair7q32$trait1,bmi_netPair7q32$traits)
bmi_netPair7q32$temp2=paste(bmi_netPair7q32$traits,bmi_netPair7q32$trait1)

# Then grab colocalization data for gene pairs in interactome
temp=as.data.frame(matrix(nrow = dim(bmi_interactome7q32)[1],ncol = 4))
for(i in 1:dim(bmi_interactome7q32)[1]){
  temp[i,1:4]=bmi_netPair7q32[ifelse(is.na(match(bmi_interactome7q32$temp1[i],bmi_netPair7q32$temp1)),
                                  match(bmi_interactome7q32$temp1[i],bmi_netPair7q32$temp2),
                                  match(bmi_interactome7q32$temp1[i],bmi_netPair7q32$temp1)),3:6]
}
temp[is.na(temp)]=0

# Then combine with the BMI colocalizations
colnames(temp)=colnames(bmi_netPair7q32)[3:6]
temp=rbind(temp,bmi_netPair7q32[bmi_netPair7q32$trait1=="BMI",3:6])

# Then add rows for BMI-Gene connections with 0 for MoA and likelihood
bmi_full7q32=bmi_interactome7q32[,1:4]
temp2=bmi_netPair7q32[bmi_netPair7q32$trait1=="BMI",1:4]
colnames(temp2)=colnames(bmi_interactome7q32)[1:4]
temp2[,3:4]=0
bmi_full7q32=rbind(bmi_full7q32,temp2)

# Finally, combine the colocalization columns with the interactome columns
bmi_full7q32=cbind(bmi_full7q32,temp)

# The nodes data also needs to be combined and duplicate rows removed
bmi_full7q32_nodes=rbind(bmi_inter_nodes7q32,bmi_coloc_nodes7q32)
bmi_full7q32_nodes=bmi_full7q32_nodes[!duplicated(bmi_full7q32_nodes$Node),]

# T2D
homair_interactome7q32$temp1=paste(homair_interactome7q32$Regulator,homair_interactome7q32$Target)
homair_interactome7q32$temp2=paste(homair_interactome7q32$Target,homair_interactome7q32$Regulator)
t2d_netPair7q32$temp1=paste(t2d_netPair7q32$trait1,t2d_netPair7q32$traits)
t2d_netPair7q32$temp2=paste(t2d_netPair7q32$traits,t2d_netPair7q32$trait1)

# Then grab colocalization data for gene pairs in interactome
temp=as.data.frame(matrix(nrow = dim(homair_interactome7q32)[1],ncol = 4))
for(i in 1:dim(homair_interactome7q32)[1]){
  temp[i,1:4]=t2d_netPair7q32[ifelse(is.na(match(homair_interactome7q32$temp1[i],t2d_netPair7q32$temp1)),
                                     match(homair_interactome7q32$temp1[i],t2d_netPair7q32$temp2),
                                     match(homair_interactome7q32$temp1[i],t2d_netPair7q32$temp1)),3:6]
}
temp[is.na(temp)]=0

# Then combine with the T2D colocalizations
colnames(temp)=colnames(t2d_netPair7q32)[3:6]
temp=rbind(temp,t2d_netPair7q32[t2d_netPair7q32$trait1=="T2D",3:6])

# Then add rows for T2D-Gene connections with 0 for MoA and likelihood
t2d_full7q32=homair_interactome7q32[,1:4]
temp2=t2d_netPair7q32[t2d_netPair7q32$trait1=="T2D",1:4]
colnames(temp2)=colnames(homair_interactome7q32)[1:4]
temp2[,3:4]=0
t2d_full7q32=rbind(t2d_full7q32,temp2)

# Finally, combine the colocalization columns with the interactome columns
t2d_full7q32=cbind(t2d_full7q32,temp)

# The nodes data also needs to be combined and duplicate rows removed
t2d_full7q32_nodes=rbind(t2d_inter_nodes7q32,t2d_coloc_nodes7q32)
t2d_full7q32_nodes=t2d_full7q32_nodes[!duplicated(t2d_full7q32_nodes$Node),]

# HDL
hdl_interactome7q32$temp1=paste(hdl_interactome7q32$Regulator,hdl_interactome7q32$Target)
hdl_interactome7q32$temp2=paste(hdl_interactome7q32$Target,hdl_interactome7q32$Regulator)
hdl_netPair7q32$temp1=paste(hdl_netPair7q32$trait1,hdl_netPair7q32$traits)
hdl_netPair7q32$temp2=paste(hdl_netPair7q32$traits,hdl_netPair7q32$trait1)

# Then grab colocalization data for gene pairs in interactome
temp=as.data.frame(matrix(nrow = dim(hdl_interactome7q32)[1],ncol = 4))
for(i in 1:dim(hdl_interactome7q32)[1]){
  temp[i,1:4]=hdl_netPair7q32[ifelse(is.na(match(hdl_interactome7q32$temp1[i],hdl_netPair7q32$temp1)),
                                     match(hdl_interactome7q32$temp1[i],hdl_netPair7q32$temp2),
                                     match(hdl_interactome7q32$temp1[i],hdl_netPair7q32$temp1)),3:6]
}
temp[is.na(temp)]=0

# Then combine with the HDL colocalizations
colnames(temp)=colnames(hdl_netPair7q32)[3:6]
temp=rbind(temp,hdl_netPair7q32[hdl_netPair7q32$trait1=="HDL",3:6])

# Then add rows for HDL-Gene connections with 0 for MoA and likelihood
hdl_full7q32=hdl_interactome7q32[,1:4]
temp2=hdl_netPair7q32[hdl_netPair7q32$trait1=="HDL",1:4]
colnames(temp2)=colnames(hdl_interactome7q32)[1:4]
temp2[,3:4]=0
hdl_full7q32=rbind(hdl_full7q32,temp2)

# Finally, combine the colocalization columns with the interactome columns
hdl_full7q32=cbind(hdl_full7q32,temp)

# The nodes data also needs to be combined and duplicate rows removed
hdl_full7q32_nodes=rbind(hdl_inter_nodes7q32,hdl_coloc_nodes7q32)
hdl_full7q32_nodes=hdl_full7q32_nodes[!duplicated(hdl_full7q32_nodes$Node),]

# TriG
triG_interactome7q32$temp1=paste(triG_interactome7q32$Regulator,triG_interactome7q32$Target)
triG_interactome7q32$temp2=paste(triG_interactome7q32$Target,triG_interactome7q32$Regulator)
triG_netPair7q32$temp1=paste(triG_netPair7q32$trait1,triG_netPair7q32$traits)
triG_netPair7q32$temp2=paste(triG_netPair7q32$traits,triG_netPair7q32$trait1)

# Then grab colocalization data for gene pairs in interactome
temp=as.data.frame(matrix(nrow = dim(triG_interactome7q32)[1],ncol = 4))
for(i in 1:dim(triG_interactome7q32)[1]){
  temp[i,1:4]=triG_netPair7q32[ifelse(is.na(match(triG_interactome7q32$temp1[i],triG_netPair7q32$temp1)),
                                     match(triG_interactome7q32$temp1[i],triG_netPair7q32$temp2),
                                     match(triG_interactome7q32$temp1[i],triG_netPair7q32$temp1)),3:6]
}
temp[is.na(temp)]=0

# Then combine with the TriG colocalizations
colnames(temp)=colnames(triG_netPair7q32)[3:6]
temp=rbind(temp,triG_netPair7q32[triG_netPair7q32$trait1=="TriG",3:6])

# Then add rows for TriG-Gene connections with 0 for MoA and likelihood
triG_full7q32=triG_interactome7q32[,1:4]
temp2=triG_netPair7q32[triG_netPair7q32$trait1=="TriG",1:4]
colnames(temp2)=colnames(triG_interactome7q32)[1:4]
temp2[,3:4]=0
triG_full7q32=rbind(triG_full7q32,temp2)

# Finally, combine the colocalization columns with the interactome columns
triG_full7q32=cbind(triG_full7q32,temp)

# The nodes data also needs to be combined and duplicate rows removed
triG_full7q32_nodes=rbind(triG_inter_nodes7q32,triG_coloc_nodes7q32)
triG_full7q32_nodes=triG_full7q32_nodes[!duplicated(triG_full7q32_nodes$Node),]

# Since LINC-PINT, KLF14 and AC016831.7 is the only cis gene here, I'll just deal with them manually
bmi_full7q32_nodes[bmi_full7q32_nodes$Node=="LINC-PINT","Best_eQTL_Beta"]=min_bmi_cisE1[min_bmi_cisE1$gene=="LINC-PINT","beta"]
bmi_full7q32_nodes[bmi_full7q32_nodes$Node=="LINC-PINT","Best_eQTL_logP"]=-log10(min_bmi_cisE1[min_bmi_cisE1$gene=="LINC-PINT","pvalue"])
bmi_full7q32_nodes[bmi_full7q32_nodes$Node=="KLF14","Best_eQTL_Beta"]=min_bmi_cisE1[min_bmi_cisE1$gene=="KLF14","beta"]
bmi_full7q32_nodes[bmi_full7q32_nodes$Node=="KLF14","Best_eQTL_logP"]=-log10(min_bmi_cisE1[min_bmi_cisE1$gene=="KLF14","pvalue"])
bmi_full7q32_nodes[bmi_full7q32_nodes$Node=="AC016831.7","Best_eQTL_Beta"]=min_bmi_cisE1[min_bmi_cisE1$gene=="AC016831.7","beta"]
bmi_full7q32_nodes[bmi_full7q32_nodes$Node=="AC016831.7","Best_eQTL_logP"]=-log10(min_bmi_cisE1[min_bmi_cisE1$gene=="AC016831.7","pvalue"])
t2d_full7q32_nodes[t2d_full7q32_nodes$Node=="LINC-PINT","Best_eQTL_Beta"]=min_t2d_cisE1[min_t2d_cisE1$gene=="LINC-PINT","beta"]
t2d_full7q32_nodes[t2d_full7q32_nodes$Node=="LINC-PINT","Best_eQTL_logP"]=-log10(min_t2d_cisE1[min_t2d_cisE1$gene=="LINC-PINT","pvalue"])
t2d_full7q32_nodes[t2d_full7q32_nodes$Node=="KLF14","Best_eQTL_Beta"]=min_t2d_cisE1[min_t2d_cisE1$gene=="KLF14","beta"]
t2d_full7q32_nodes[t2d_full7q32_nodes$Node=="KLF14","Best_eQTL_logP"]=-log10(min_t2d_cisE1[min_t2d_cisE1$gene=="KLF14","pvalue"])
t2d_full7q32_nodes[t2d_full7q32_nodes$Node=="AC016831.7","Best_eQTL_Beta"]=min_t2d_cisE1[min_t2d_cisE1$gene=="AC016831.7","beta"]
t2d_full7q32_nodes[t2d_full7q32_nodes$Node=="AC016831.7","Best_eQTL_logP"]=-log10(min_t2d_cisE1[min_t2d_cisE1$gene=="AC016831.7","pvalue"])
hdl_full7q32_nodes[hdl_full7q32_nodes$Node=="LINC-PINT","Best_eQTL_Beta"]=min_hdl_cisE1[min_hdl_cisE1$gene=="LINC-PINT","beta"]
hdl_full7q32_nodes[hdl_full7q32_nodes$Node=="LINC-PINT","Best_eQTL_logP"]=-log10(min_hdl_cisE1[min_hdl_cisE1$gene=="LINC-PINT","pvalue"])
hdl_full7q32_nodes[hdl_full7q32_nodes$Node=="KLF14","Best_eQTL_Beta"]=min_hdl_cisE1[min_hdl_cisE1$gene=="KLF14","beta"]
hdl_full7q32_nodes[hdl_full7q32_nodes$Node=="KLF14","Best_eQTL_logP"]=-log10(min_hdl_cisE1[min_hdl_cisE1$gene=="KLF14","pvalue"])
hdl_full7q32_nodes[hdl_full7q32_nodes$Node=="AC016831.7","Best_eQTL_Beta"]=min_hdl_cisE1[min_hdl_cisE1$gene=="AC016831.7","beta"]
hdl_full7q32_nodes[hdl_full7q32_nodes$Node=="AC016831.7","Best_eQTL_logP"]=-log10(min_hdl_cisE1[min_hdl_cisE1$gene=="AC016831.7","pvalue"])
triG_full7q32_nodes[triG_full7q32_nodes$Node=="LINC-PINT","Best_eQTL_Beta"]=min_triG_cisE1[min_triG_cisE1$gene=="LINC-PINT","beta"]
triG_full7q32_nodes[triG_full7q32_nodes$Node=="LINC-PINT","Best_eQTL_logP"]=-log10(min_triG_cisE1[min_triG_cisE1$gene=="LINC-PINT","pvalue"])
triG_full7q32_nodes[triG_full7q32_nodes$Node=="KLF14","Best_eQTL_Beta"]=min_triG_cisE1[min_triG_cisE1$gene=="KLF14","beta"]
triG_full7q32_nodes[triG_full7q32_nodes$Node=="KLF14","Best_eQTL_logP"]=-log10(min_triG_cisE1[min_triG_cisE1$gene=="KLF14","pvalue"])
triG_full7q32_nodes[triG_full7q32_nodes$Node=="AC016831.7","Best_eQTL_Beta"]=min_triG_cisE1[min_triG_cisE1$gene=="AC016831.7","beta"]
triG_full7q32_nodes[triG_full7q32_nodes$Node=="AC016831.7","Best_eQTL_logP"]=-log10(min_triG_cisE1[min_triG_cisE1$gene=="AC016831.7","pvalue"])

# Final touches by replacing NA with 0
bmi_full7q32_nodes[is.na(bmi_full7q32_nodes)]=0
t2d_full7q32_nodes[is.na(t2d_full7q32_nodes)]=0
hdl_full7q32_nodes[is.na(hdl_full7q32_nodes)]=0
triG_full7q32_nodes[is.na(triG_full7q32_nodes)]=0

# Write networks and node data to file for Cytoscape visualizations
write.table(bmi_full7q32,"./BMI/Chr7q32_cis-Genes_and_BMI_MRs_interactome_and_pairwise_colocalization.txt",sep = "\t",quote = F,row.names = F)
write.table(bmi_full7q32_nodes,"./BMI/Chr7q32_cis-Genes_and_BMI_MRs_interactome_and_pairwise_colocalization_node_info.txt",sep = "\t",quote = F,row.names = F)
write.table(t2d_full7q32,"./T2D/Chr7q32_cis-Genes_and_HOMA-IR_MRs_interactome_and_pairwise_colocalization.txt",sep = "\t",quote = F,row.names = F)
write.table(t2d_full7q32_nodes,"./T2D/Chr7q32_cis-Genes_and_HOMA-IR_MRs_interactome_and_pairwise_colocalization_node_info.txt",sep = "\t",quote = F,row.names = F)
write.table(hdl_full7q32,"./HDL/Chr7q32_cis-Genes_and_HDL_MRs_interactome_and_pairwise_colocalization.txt",sep = "\t",quote = F,row.names = F)
write.table(hdl_full7q32_nodes,"./HDL/Chr7q32_cis-Genes_and_HDL_MRs_interactome_and_pairwise_colocalization_node_info.txt",sep = "\t",quote = F,row.names = F)
write.table(triG_full7q32,"./Triglycerides/Chr7q32_cis-Genes_and_Triglycerides_MRs_interactome_and_pairwise_colocalization.txt",sep = "\t",quote = F,row.names = F)
write.table(triG_full7q32_nodes,"./Triglycerides/Chr7q32_cis-Genes_and_Triglycerides_MRs_interactome_and_pairwise_colocalization_node_info.txt",sep = "\t",quote = F,row.names = F)

# 12p13.1
interactome12p13=bmi_MRMR

# Grab pairwise colocalizations with PP>0.5 between BMI and QTLs and EPHB2 aQTL and trans-QTLs
bmi_pairColoc12p13=bmi_pairColoc[bmi_pairColoc$locus=="12p13.1" & bmi_pairColoc$posterior_prob>0.5,c(2,3,5)]
bmi_pairColoc12p13$traits=gsub("BMI, ","",bmi_pairColoc12p13$traits)
bmi_pairColoc12p13=cbind("trait1"=rep("BMI",dim(bmi_pairColoc12p13)[1]),bmi_pairColoc12p13)
bmi_e_pairColoc12p13=bmi_pairColoc12p13[grepl("-e_",bmi_pairColoc12p13$traits),]
bmi_a_pairColoc12p13=bmi_pairColoc12p13[grepl("-a_",bmi_pairColoc12p13$traits),]
bmi_e_pairColoc12p13$traits=gsub(".*_","",bmi_e_pairColoc12p13$traits)
bmi_a_pairColoc12p13$traits=gsub(".*_","",bmi_a_pairColoc12p13$traits)
bmi_netPair12p13=rbind(bmi_e_pairColoc12p13,bmi_a_pairColoc12p13)
bmi_netPair12p13=bmi_netPair12p13[!duplicated(bmi_netPair12p13$traits),-c(3,4)]
bmi_netPair12p13$eQTL_PP=bmi_e_pairColoc12p13[match(bmi_netPair12p13$traits,bmi_e_pairColoc12p13$traits),3]
bmi_netPair12p13$eQTL_SNP=bmi_e_pairColoc12p13[match(bmi_netPair12p13$traits,bmi_e_pairColoc12p13$traits),4]
bmi_netPair12p13$aQTL_PP=bmi_a_pairColoc12p13[match(bmi_netPair12p13$traits,bmi_a_pairColoc12p13$traits),3]
bmi_netPair12p13$aQTL_SNP=bmi_a_pairColoc12p13[match(bmi_netPair12p13$traits,bmi_a_pairColoc12p13$traits),4]

colocNet12p13=bmi_netPair12p13
colocNet12p13[is.na(colocNet12p13)]=0

# Grab the -log10(Pmin) and betas for the eQTLs and aQTLs among BMI GWAS significant SNPs at the 12p13.1 locus
transE4=filt_trans_bmi_eqtl[filt_trans_bmi_eqtl$chr==12 & filt_trans_bmi_eqtl$position>13900000 & filt_trans_bmi_eqtl$position<15000000,]
transE4=transE4[transE4$snps %in% sig_bmi$SNP,]
transE4=transE4[order(transE4$pvalue),]
min_transE4=transE4[!duplicated(transE4$gene),]

transA4=filt_trans_bmi_aqtl[filt_trans_bmi_aqtl$chr==12 & filt_trans_bmi_aqtl$position>13900000 & filt_trans_bmi_aqtl$position<15000000,]
transA4=transA4[transA4$snps %in% sig_bmi$SNP,]
transA4=transA4[order(transA4$pvalue),]
min_transA4=transA4[!duplicated(transA4$gene),]

# Make node tables for 12p13 networks
inter_nodes12p13=data.frame("Node"=as.character(unique(interactome12p13$Regulator)),"BMI_exp_cor"=rep(0,length(unique(interactome12p13$Regulator))),
                           "BMI_act_cor"=rep(0,length(unique(interactome12p13$Regulator))),"Best_eQTL_Beta"=rep(0,length(unique(interactome12p13$Regulator))),
                           "Best_eQTL_logP"=rep(0,length(unique(interactome12p13$Regulator))),"Best_aQTL_Beta"=rep(0,length(unique(interactome12p13$Regulator))),
                           "Best_aQTL_logP"=rep(0,length(unique(interactome12p13$Regulator))))
for(i in 1:dim(inter_nodes12p13)[1]){
  inter_nodes12p13$BMI_exp_cor[i]=cor(as.numeric(tpm[as.character(inter_nodes12p13$Node[i]),]),filt_pheno$BMI)
  inter_nodes12p13$BMI_act_cor[i]=cor(as.numeric(vip[as.character(inter_nodes12p13$Node[i]),]),filt_pheno$BMI)
  inter_nodes12p13$Best_eQTL_Beta[i]=ifelse(inter_nodes12p13$Node[i] %in% min_transE4$gene,
                                           min_transE4[min_transE4$gene==as.character(inter_nodes12p13$Node[i]),"beta"],
                                           0)
  inter_nodes12p13$Best_eQTL_logP[i]=ifelse(inter_nodes12p13$Node[i] %in% min_transE4$gene,
                                           -log10(min_transE4[min_transE4$gene==as.character(inter_nodes12p13$Node[i]),"pvalue"]),
                                           0)
  inter_nodes12p13$Best_aQTL_Beta[i]=ifelse(inter_nodes12p13$Node[i] %in% min_transA4$gene,
                                           min_transA4[min_transA4$gene==as.character(inter_nodes12p13$Node[i]),"beta"],
                                           0)
  inter_nodes12p13$Best_aQTL_logP[i]=ifelse(inter_nodes12p13$Node[i] %in% min_transA4$gene,
                                           -log10(min_transA4[min_transA4$gene==as.character(inter_nodes12p13$Node[i]),"pvalue"]),
                                           0)
}

coloc_nodes12p13=data.frame("Node"=as.character(unique(colocNet12p13$traits)),"BMI_exp_cor"=rep(0,length(unique(colocNet12p13$traits))),
                           "BMI_act_cor"=rep(0,length(unique(colocNet12p13$traits))),"Best_eQTL_Beta"=rep(0,length(unique(colocNet12p13$traits))),
                           "Best_eQTL_logP"=rep(0,length(unique(colocNet12p13$traits))),"Best_aQTL_Beta"=rep(0,length(unique(colocNet12p13$traits))),
                           "Best_aQTL_logP"=rep(0,length(unique(colocNet12p13$traits))))
for(i in 1:dim(coloc_nodes12p13)[1]){
  coloc_nodes12p13$BMI_exp_cor[i]=cor(as.numeric(tpm[as.character(coloc_nodes12p13$Node[i]),]),filt_pheno$BMI)
  coloc_nodes12p13$BMI_act_cor[i]=cor(as.numeric(vip[as.character(coloc_nodes12p13$Node[i]),]),filt_pheno$BMI)
  coloc_nodes12p13$Best_eQTL_Beta[i]=ifelse(coloc_nodes12p13$Node[i] %in% min_transE4$gene,
                                           min_transE4[min_transE4$gene==as.character(coloc_nodes12p13$Node[i]),"beta"],
                                           0)
  coloc_nodes12p13$Best_eQTL_logP[i]=ifelse(coloc_nodes12p13$Node[i] %in% min_transE4$gene,
                                           -log10(min_transE4[min_transE4$gene==as.character(coloc_nodes12p13$Node[i]),"pvalue"]),
                                           0)
  coloc_nodes12p13$Best_aQTL_Beta[i]=ifelse(coloc_nodes12p13$Node[i] %in% min_transA4$gene,
                                           min_transA4[min_transA4$gene==as.character(coloc_nodes12p13$Node[i]),"beta"],
                                           0)
  coloc_nodes12p13$Best_aQTL_logP[i]=ifelse(coloc_nodes12p13$Node[i] %in% min_transA4$gene,
                                           -log10(min_transA4[min_transA4$gene==as.character(coloc_nodes12p13$Node[i]),"pvalue"]),
                                           0)
}

# I think it may be more convenient to merge the networks into one and then just change which attributes I visualize in Cytoscape
# Start with 2 temporary columns concatinating the regulator-target and target-regulator for easier matching.
interactome12p13$temp1=paste(interactome12p13$Regulator,interactome12p13$Target)
interactome12p13$temp2=paste(interactome12p13$Target,interactome12p13$Regulator)
colocNet12p13$temp1=paste(colocNet12p13$trait1,colocNet12p13$traits)
colocNet12p13$temp2=paste(colocNet12p13$traits,colocNet12p13$trait1)

# Then grab colocalization data for gene pairs in interactome
temp=as.data.frame(matrix(nrow = dim(interactome12p13)[1],ncol = 4))
for(i in 1:dim(interactome12p13)[1]){
  temp[i,1:4]=colocNet12p13[ifelse(is.na(match(interactome12p13$temp1[i],colocNet12p13$temp1)),
                                  match(interactome12p13$temp1[i],colocNet12p13$temp2),
                                  match(interactome12p13$temp1[i],colocNet12p13$temp1)),3:6]
}
temp[is.na(temp)]=0

# Then combine with the BMI colocalizations
colnames(temp)=colnames(colocNet12p13)[3:6]
temp=rbind(temp,colocNet12p13[colocNet12p13$trait1=="BMI",3:6])

# Then add rows for BMI-Gene connections with 0 for MoA and likelihood
full12p13=interactome12p13[,1:4]
temp2=colocNet12p13[colocNet12p13$trait1=="BMI",1:4]
colnames(temp2)=colnames(interactome12p13)[1:4]
temp2[,3:4]=0
full12p13=rbind(full12p13,temp2)

# Finally, combine the colocalization columns with the interactome columns
full12p13=cbind(full12p13,temp)

# The nodes data also needs to be combined and duplicate rows removed
full12p13_nodes=rbind(inter_nodes12p13,coloc_nodes12p13)
full12p13_nodes=full12p13_nodes[!duplicated(full12p13_nodes$Node),]

# Write networks and node data to file for Cytoscape visualizations
write.table(interactome12p13,"Chr12p13_BMI_MRs_interactome.txt",sep = "\t",quote = F,row.names = F)
write.table(inter_nodes12p13,"Chr12p13_BMI_MRs_interactome_node_info.txt",sep = "\t",quote = F,row.names = F)
write.table(colocNet12p13,"Chr12p13_BMI_and_BMI_MRs_pairwise_colocalization_network.txt",sep = "\t",quote = F,row.names = F)
write.table(coloc_nodes12p13,"Chr12p13_BMI_and_BMI_MRs_pairwise_colocalization_network_node_info.txt",sep = "\t",quote = F,row.names = F)
write.table(full12p13,"Chr12p13_BMI_MRs_interactome_and_pairwise_colocalization.txt",sep = "\t",quote = F,row.names = F)
write.table(full12p13_nodes,"Chr12p13_BMI_MRs_interactome_and_pairwise_colocalization_node_info.txt",sep = "\t",quote = F,row.names = F)

  










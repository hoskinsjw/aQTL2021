### This script is for generating violin plots for select eQTLs and aQTLs for BMI loci

library(ggplot2)
library(RColorBrewer)
library(snowfall)
library(foreach)
library(doParallel)
library(Hmisc)

setwd("YOUR WORKING DIRECTORY")

# Read in the LD matrices, and the GWAS and QTL data. The file locations are relative to your working directory, so adjust accordingly.
exp=read.table("./Filtered_Eurobats_adipose_qnorm_INT_logTPMs_for_4213_regulators.txt",
               sep = "\t",header = T,row.names = 1)
vip=read.table("./Filtered_Eurobats_adipose_unnormalized_activities_from_logTPM_for_4213_regulators.txt",
               sep = "\t",header = T,row.names = 1)
bmi1p36=read.table("./Eurobats_chr1p36.1_genotypes_for_colocalizations.dosage",
               sep = "\t",header = T,row.names = 1)
bmi7q32=read.table("./Eurobats_chr7q32_genotypes_for_colocalizations.dosage",
                   sep = "\t",header = T,row.names = 1)
bmi12p13=read.table("./Eurobats_chr12p13.1_genotypes_for_colocalizations.dosage",
                   sep = "\t",header = T,row.names = 1)
all(colnames(exp)==colnames(vip)) # TRUE
all(colnames(bmi1p36)==colnames(vip)) # TRUE

# Grab relevant SNPs and genes
filt_1p36=as.data.frame(round(t(bmi1p36[c("rs6692586","rs4654828"),]))) # Need to round imputed dosages to nearest whole number for plots
filt_7q32=as.data.frame(round(t(bmi7q32[c("rs972283","rs287621","rs738134"),]))) # Need to round imputed dosages to nearest whole number for plots
filt_12p13=as.data.frame(round(t(bmi12p13[c("rs12422552"),]))) # Need to round imputed dosages to nearest whole number for plots
filt_gwas=cbind(filt_1p36,filt_7q32,filt_12p13)
filt_exp=as.data.frame(t(exp[c("EPHB2","DAPK2","GNA14","LASP1","RASSF4","LINC-PINT","KLF14","TBX4","NR2F1","AGT","ID2","PTPRJ"),]))
filt_vip=as.data.frame(t(vip[c("EPHB2","DAPK2","GNA14","LASP1","RASSF4","LINC-PINT","KLF14","TBX4","NR2F1","AGT","ID2","PTPRJ"),]))

# Z-scale expression and activity values
cl=makeCluster(6)
registerDoParallel(cl)
exp_z=list()
exp_z=foreach(i=1:dim(filt_exp)[2]) %dopar%
  as.numeric(scale(as.numeric(filt_exp[,i])))
exp_z=as.data.frame(structure(exp_z, row.names = c(NA, -length(exp_z[[1]])), class = "data.frame")) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(exp_z)=colnames(filt_exp)
rownames(exp_z)=rownames(filt_exp)

cl=makeCluster(6)
registerDoParallel(cl)
vip_z=list()
vip_z=foreach(i=1:dim(filt_vip)[2]) %dopar%
  as.numeric(scale(as.numeric(filt_vip[,i])))
vip_z=as.data.frame(structure(vip_z, row.names = c(NA, -length(vip_z[[1]])), class = "data.frame")) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(vip_z)=colnames(filt_vip)
rownames(vip_z)=rownames(filt_vip)

filt_exp=cbind(exp_z,filt_gwas,"Type"=rep("Expression",dim(filt_gwas)[1]))
filt_vip=cbind(vip_z,filt_gwas,"Type"=rep("Activity",dim(filt_gwas)[1]))
both=rbind(filt_exp,filt_vip)

# Make sure the dosage is consistently minor allele dosage
table(both$rs6692586) # major allelic dosage
table(both$rs4654828) # major allelic dosage
table(both$rs972283) # minor allelic dosage
table(both$rs287621) # major allelic dosage
table(both$rs738134) # minor allelic dosage
table(both$rs12422552) # minor allelic dosage
both$rs6692586=abs(both$rs6692586-2)
both$rs4654828=abs(both$rs4654828-2)
both$rs287621=abs(both$rs287621-2)

# Need to convert the allelic dosages to factors for the violin plot
both$rs6692586=as.factor(both$rs6692586)
both$rs4654828=as.factor(both$rs4654828)
both$rs972283=as.factor(both$rs972283)
both$rs287621=as.factor(both$rs287621)
both$rs738134=as.factor(both$rs738134)
both$rs12422552=as.factor(both$rs12422552)

# The hyphen in LINC-PINT is causing trouble in calling its column. Let's change it to a period.
colnames(both)[6]="LINC.PINT"

# Make the violin plots
ephb2=ggplot(both,aes(x=rs4654828,y=EPHB2,fill=Type)) + 
  geom_violin(trim=F,position = position_dodge(0.8),width = 0.7) +
  scale_fill_manual(breaks = c("Expression","Activity"),values = c("#0099FF","#FF6666")) +
  labs(x="rs4654828\n(minor allelic dosage)",y="Scaled EPHB2 expression or activity") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic()

dapk2=ggplot(both,aes(x=rs4654828,y=DAPK2,fill=Type)) + 
  geom_violin(trim=F,position = position_dodge(0.8),width = 0.7) +
  scale_fill_manual(breaks = c("Expression","Activity"),values = c("#0099FF","#FF6666")) +
  labs(x="rs4654828\n(minor allelic dosage)",y="Scaled DAPK2 expression or activity") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic()

gna14=ggplot(both,aes(x=rs4654828,y=GNA14,fill=Type)) + 
  geom_violin(trim=F,position = position_dodge(0.8),width = 0.7) +
  scale_fill_manual(breaks = c("Expression","Activity"),values = c("#0099FF","#FF6666")) +
  labs(x="rs4654828\n(minor allelic dosage)",y="Scaled GNA14 expression or activity") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic()

lasp1=ggplot(both,aes(x=rs4654828,y=LASP1,fill=Type)) + 
  geom_violin(trim=F,position = position_dodge(0.8),width = 0.7) +
  scale_fill_manual(breaks = c("Expression","Activity"),values = c("#0099FF","#FF6666")) +
  labs(x="rs4654828\n(minor allelic dosage)",y="Scaled LASP1 expression or activity") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic()

rassf4=ggplot(both,aes(x=rs4654828,y=RASSF4,fill=Type)) + 
  geom_violin(trim=F,position = position_dodge(0.8),width = 0.7) +
  scale_fill_manual(breaks = c("Expression","Activity"),values = c("#0099FF","#FF6666")) +
  labs(x="rs4654828\n(minor allelic dosage)",y="Scaled RASSF4 expression or activity") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic()

linc=ggplot(both,aes(x=rs972283,y=both$LINC.PINT,fill=Type)) + 
  geom_violin(trim=F,position = position_dodge(0.8),width = 0.7) +
  scale_fill_manual(breaks = c("Expression","Activity"),values = c("#0099FF","#FF6666")) +
  labs(x="rs972283\n(minor allelic dosage)",y="Scaled LINC-PINT expression or activity") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic()

klf14=ggplot(both,aes(x=rs972283,y=KLF14,fill=Type)) + 
  geom_violin(trim=F,position = position_dodge(0.8),width = 0.7) +
  scale_fill_manual(breaks = c("Expression","Activity"),values = c("#0099FF","#FF6666")) +
  labs(x="rs972283\n(minor allelic dosage)",y="Scaled KLF14 expression or activity") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic()

tbx4=ggplot(both,aes(x=rs972283,y=TBX4,fill=Type)) + 
  geom_violin(trim=F,position = position_dodge(0.8),width = 0.7) +
  scale_fill_manual(breaks = c("Expression","Activity"),values = c("#0099FF","#FF6666")) +
  labs(x="rs972283\n(minor allelic dosage)",y="Scaled TBX4 expression or activity") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic()
  
nr2f1=ggplot(both,aes(x=rs972283,y=NR2F1,fill=Type)) + 
  geom_violin(trim=F,position = position_dodge(0.8),width = 0.7) +
  scale_fill_manual(breaks = c("Expression","Activity"),values = c("#0099FF","#FF6666")) +
  labs(x="rs972283\n(minor allelic dosage)",y="Scaled NR2F1 expression or activity") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic()

agt=ggplot(both,aes(x=rs972283,y=AGT,fill=Type)) + 
  geom_violin(trim=F,position = position_dodge(0.8),width = 0.7) +
  scale_fill_manual(breaks = c("Expression","Activity"),values = c("#0099FF","#FF6666")) +
  labs(x="rs972283\n(minor allelic dosage)",y="Scaled AGT expression or activity") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic()

id2=ggplot(both,aes(x=rs12422552,y=ID2,fill=Type)) + 
  geom_violin(trim=F,position = position_dodge(0.8),width = 0.7) +
  scale_fill_manual(breaks = c("Expression","Activity"),values = c("#0099FF","#FF6666")) +
  labs(x="rs12422552\n(minor allelic dosage)",y="Scaled ID2 expression or activity") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic()

ptprj=ggplot(both,aes(x=rs12422552,y=PTPRJ,fill=Type)) + 
  geom_violin(trim=F,position = position_dodge(0.8),width = 0.7) +
  scale_fill_manual(breaks = c("Expression","Activity"),values = c("#0099FF","#FF6666")) +
  labs(x="rs12422552\n(minor allelic dosage)",y="Scaled PTPRJ expression or activity") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic()

# Write the plots to files
pdf("rs4654828-EPHB2_violin_plot.pdf",width = 8)
ephb2
dev.off()

pdf("rs4654828-DAPK2_violin_plot.pdf",width = 8)
dapk2
dev.off()

pdf("rs4654828-GNA14_violin_plot.pdf",width = 8)
gna14
dev.off()

pdf("rs4654828-LASP1_violin_plot.pdf",width = 8)
lasp1
dev.off()

pdf("rs4654828-RASSF4_violin_plot.pdf",width = 8)
rassf4
dev.off()

pdf("rs972283-LINC-PINT_violin_plot.pdf",width = 8)
linc
dev.off()

pdf("rs972283-KLF14_violin_plot.pdf",width = 8)
klf14
dev.off()

pdf("rs972283-TBX4_violin_plot.pdf",width = 8)
tbx4
dev.off()

pdf("rs972283-NR2F1_violin_plot.pdf",width = 8)
nr2f1
dev.off()

pdf("rs972283-AGT_violin_plot.pdf",width = 8)
agt
dev.off()

pdf("rs12422552-ID2_violin_plot.pdf",width = 8)
id2
dev.off()

pdf("rs12422552-PTPRJ_violin_plot.pdf",width = 8)
ptprj
dev.off()




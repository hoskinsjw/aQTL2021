library(RColorBrewer)
library(ggplot2)
library(ggpolypath)
library(venn)

setwd("YOUR WORKING DIRECTORY")

# Read in the LD matrices, and the GWAS and QTL data. The file locations are relative to your working directory, so adjust accordingly.
mrs=read.table("./MR_comparison_table.txt",sep = "\t",header = T,row.names = 1)
mrs_exp=read.table("./MRs_BY_EXPRESSION_comparison_table.txt",sep = "\t",header = T,row.names = 1)
exp=read.table("./Filtered_Eurobats_adipose_qnorm_INT_logTPMs_for_4213_regulators.txt",
               sep = "\t",header = T,row.names = 1)
vip=read.table("./Filtered_Eurobats_adipose_unnormalized_activities_from_logTPM_for_4213_regulators.txt",
               sep = "\t",header = T,row.names = 1)

### Make a numeric/logical version of mrs
# For MRs by activity
mrs_logic=mrs
for(i in 1:dim(mrs_logic)[2]){
  mrs_logic[,i]=sapply(mrs_logic[,i],function(x) x=="TRUE")
  mrs_logic[,i]=as.numeric(mrs_logic[,i])
}

# For MRs by expression
mrs_exp_logic=mrs_exp
for(i in 1:dim(mrs_exp_logic)[2]){
  mrs_exp_logic[,i]=sapply(mrs_exp_logic[,i],function(x) x=="TRUE")
  mrs_exp_logic[,i]=as.numeric(mrs_exp_logic[,i])
}

### Filter down to just the MRs
# For MRs by activity
mr_exp=exp[na.omit(match(rownames(mrs_logic),rownames(exp))),]
mr_vip=vip[na.omit(match(rownames(mrs_logic),rownames(vip))),]

# For MRs by expression
mr_by_exp_exp=exp[na.omit(match(rownames(mrs_exp_logic),rownames(exp))),]
mr_by_exp_vip=vip[na.omit(match(rownames(mrs_exp_logic),rownames(vip))),]

### Make correlation matrices
# For MRs by activity
exp_exp=cor(t(mr_exp),t(mr_exp))
vip_vip=cor(t(mr_vip),t(mr_vip))
exp_vip=cor(t(mr_exp),t(mr_vip)) 
# Note that for exp_vip each correlation coefficient is between the expression values of the regulator indicted by the rowname
# and the activity scores of the regulator indicated by colname. Therefore, the correlation matrix is not symmetrical
# about the diagonal since, for instance, the correlation between DMRT2 expression and BMP3 activity will not be 
# identical to the correlation between BMP3 expression and DMRT2 activity.

# For MRs by expression
exp_exp_by_exp=cor(t(mr_by_exp_exp),t(mr_by_exp_exp))
vip_vip_by_exp=cor(t(mr_by_exp_vip),t(mr_by_exp_vip))
exp_vip_by_exp=cor(t(mr_by_exp_exp),t(mr_by_exp_vip))
# Note that for exp_vip each correlation coefficient is between the expression values of the regulator indicted by the rowname
# and the activity scores of the regulator indicated by colname. Therefore, the correlation matrix is not symmetrical
# about the diagonal since, for instance, the correlation between DMRT2 expression and BMP3 activity will not be 
# identical to the correlation between BMP3 expression and DMRT2 activity.

# For MRs by activity versus MRs by expression
exp_by_exp_exp_by_vip=cor(t(mr_by_exp_exp),t(mr_exp))
vip_by_exp_vip_by_vip=cor(t(mr_by_exp_vip),t(mr_vip))
exp_by_exp_vip_by_vip=cor(t(mr_by_exp_exp),t(mr_vip))
exp_by_vip_vip_by_exp=cor(t(mr_exp),t(mr_by_exp_vip))

# Expression-activity correlation for matching vs non-matching MRs and MRs vs non-MRs
# For MRs by expression
match_by_exp=c()
for(i in 1:dim(mrs_exp)[1]){
  match_by_exp[i]=exp_vip_by_exp[i,i]
  names(match_by_exp)[i]=rownames(exp_vip_by_exp)[i]
}

no_match_by_exp=c()
for(i in 1:dim(mrs_exp)[1]){
  no_match_by_exp=append(no_match_by_exp,exp_vip_by_exp[-i,i])
}

# For MRs by activity
match_by_vip=c()
for(i in 1:dim(mrs)[1]){
  match_by_vip[i]=exp_vip[i,i]
  names(match_by_vip)[i]=rownames(exp_vip)[i]
}

no_match_by_vip=c()
for(i in 1:dim(mrs)[1]){
  no_match_by_vip=append(no_match_by_vip,exp_vip[-i,i])
}

### Let's check the distributions of correlation coefficients
par(mfrow=c(3,2))
plot(density(exp_exp_by_exp),xlim=c(-1.2,1.2),ylim=c(0,1.2),main="Exp vs Exp (MRs by expression)")
plot(density(vip_vip_by_exp),xlim=c(-1.2,1.2),ylim=c(0,1.2),main="Act vs Act (MRs by expression)")
plot(density(exp_exp),xlim=c(-1.2,1.2),ylim=c(0,1.2),main="Exp vs Exp (MRs by activity)")
plot(density(vip_vip),xlim=c(-1.2,1.2),ylim=c(0,1.2),main="Act vs Act (MRs by activity)")
plot(density(exp_by_exp_exp_by_vip),xlim=c(-1.2,1.2),ylim=c(0,1.2),main="Exp vs Exp (MRs by expression vs activity)")
plot(density(vip_by_exp_vip_by_vip),xlim=c(-1.2,1.2),ylim=c(0,1.2),main="Act vs Act (MRs by expression vs activity)")

par(mfrow=c(2,2))
plot(density(exp_vip_by_exp),xlim=c(-1.2,1.2),ylim=c(0,1.2),main="Exp vs Act (MRs by expression)")
plot(density(exp_vip),xlim=c(-1.2,1.2),ylim=c(0,1.2),main="Exp vs Act (MRs by activity)")
plot(density(exp_by_exp_vip_by_vip),xlim=c(-1.2,1.2),ylim=c(0,1.2),main="Exp vs Act (MRs by expression vs activity, respectively)")
plot(density(exp_by_vip_vip_by_exp),xlim=c(-1.2,1.2),ylim=c(0,1.2),main="Act vs Exp (MRs by expression vs activity, respectively)")

par(mfrow=c(2,2))
plot(density(match_by_exp),xlim=c(-1.2,1.2),ylim=c(0,4.5),main="Exp vs Act for matching MRs by expression")
plot(density(match_by_vip),xlim=c(-1.2,1.2),ylim=c(0,4.5),main="Exp vs Act for matching MRs by activity")
plot(density(no_match_by_exp),xlim=c(-1.2,1.2),ylim=c(0,4.5),main="Exp vs Act for non-matching MRs by expression")
plot(density(no_match_by_vip),xlim=c(-1.2,1.2),ylim=c(0,4.5),main="Exp vs Act for non-matching MRs by activity")

# Let's make some better quality images for the Exp vs Act matching and unmatching for MRs by activity
matched=ggplot(as.data.frame(match_by_vip),aes(match_by_vip)) +
  geom_density() + coord_cartesian(xlim = c(0, 1), ylim = c(0, 3)) + labs(x="Pearson correlation coefficient",y="Density") + theme_classic()
matched

unmatched=ggplot(as.data.frame(no_match_by_vip),aes(no_match_by_vip)) +
  geom_density() + coord_cartesian(xlim = c(-1, 1), ylim = c(0, 3)) + labs(x="Pearson correlation coefficient",y="Density") + theme_classic()
unmatched

pdf("Corr_distribution_for_matched_expression_and_activity_for_MRs.pdf",width = 8)
matched
dev.off()

pdf("Corr_distribution_for_unmatched_expression_and_activity_for_MRs.pdf",width = 8)
unmatched
dev.off()


### Grab some summary stats
summary(match_by_exp)
summary(match_by_vip)
summary(no_match_by_exp)
summary(no_match_by_vip)

# The following will get the modes of the bimodal density distributions for the no_match plots
# For MRs by expression
which.max(density(no_match_by_exp)$y) # 130 is the index that will be used to find the x-coord of the highest peak
density(no_match_by_exp)$x[130] # Highest mode = -0.4892589
MaxY=max(density(no_match_by_exp)$y[density(no_match_by_exp)$x>0]) # This finds the maxY for the distribution for x > 0
which(density(no_match_by_exp)$y==MaxY) # 374 is the index that will be used to find the x-coord of the other peak
density(no_match_by_exp)$x[374] # Second mode = 0.4981385

# For MRs by activity
which.max(density(no_match_by_vip)$y) # 138 is the index that will be used to find the x-coord of the highest peak
density(no_match_by_vip)$x[138] # Highest mode = -0.4451985
MaxY=max(density(no_match_by_vip)$y[density(no_match_by_vip)$x>0]) # This finds the maxY for the distribution for x > 0
which(density(no_match_by_vip)$y==MaxY) # 359 is the index that will be used to find the x-coord of the other peak
density(no_match_by_vip)$x[359] # Second mode = 0.4332423

### Now let's work on clustering with heatmaps
# Set up the color palette I will use for the heatmaps
RdBu=rev(brewer.pal(11,"RdBu")) # rev() reverses the vector order so the palette goes from blue to red rather than from red to blue
pdf("Heatmap_color_palette.pdf")
display.brewer.pal(7,"RdBu")
dev.off()

# Make heatmaps from the correlation matrices for the MRs inferred by activities.
pdf("MR_expression_vs_expression_correlation_heatmap.pdf",width = 8,height = 8)
EEheat=heatmap(exp_exp,symm = T,scale = "none",col = RdBu,labRow = F, labCol = F)
dev.off()
pdf("MR_activity_vs_activity_correlation_heatmap.pdf",width = 8,height = 8)
AAheat=heatmap(vip_vip,symm = T,scale = "none",col = RdBu,labRow = F, labCol = F)
dev.off()
pdf("MR_expression_(row_regulator)_vs_activity_(column_regulator)_correlation_heatmap_with_symmetrical_regulator_clustering_by_columns.pdf",width = 8,height = 8)
EAheat=heatmap(exp_vip,symm = T,scale = "none",col = RdBu,labRow = F, labCol = F)
dev.off()
# Note that each correlation coefficient is between the expression values of the regulator indicated by row and the 
# activity scores of the regulator indicated by column.

# Make MR label graphs according to the heatmap clustering
EElab=mrs_logic[EEheat$rowInd,]
pdf("MR_expression_vs_expression_block_clusters.pdf")
EEgraph=heatmap(as.matrix(EElab),Rowv = NA, Colv = NA,scale = "none",col = c("white","black"),labRow = F,labCol = F)
dev.off()

AAlab=mrs_logic[AAheat$rowInd,]
pdf("MR_activity_vs_activity_block_clusters.pdf")
AAgraph=heatmap(as.matrix(AAlab),Rowv = NA, Colv = NA,scale = "none",col = c("white","black"),labRow = F,labCol = F)
dev.off()

# Since the clustering was performed by column and applied to the rows as well, due to the symm argument being set to 
# TRUE, we need to get the regulator order for the heatmap columns.
EAcols=mrs_logic[EAheat$colInd,]
pdf("MR_expression_vs_activity_block_clusters_by_columns.pdf")
EAgraph=heatmap(as.matrix(EAcols),Rowv = NA, Colv = NA,scale = "none",col = c("white","black"),labRow = F,labCol = F)
dev.off()

### Venn diagram
# After much playing around with different packages and online tools, I decided it would work best to manually make my 
# own in Powerpoint. However, I did get the basic idea for a good layout from the "venn" package, which features avocado 
# shaped categories arranged in a star pattern.

### Further comparison between MRs inferred from expression versus MRs inferred from activities.
# Based on visual inspection of the correlation coefficient distributions for Exp-Act for matching MRs by expression
# or by activity, it appears that the MRs inferred from expression tend to have higher correlation between their
# expression and activity profiles than for the MRs inferred from activities. One hypothesis for this is that the
# regulators that are best for predicting the given phenotype based on expression are more likely to have strong
# correlation between expression and activity. In other words, regulators whose expression values are well correlated
# with phenotype while their activities are not are less likely to serve as robust predictors by expression. Therefore, 
# the MR analysis by expression seems to bias the selection to regulators that have strongly correlated expression and
# activity. Conversely, MR inference by activity would be relatively agnostic of the regulator's own expression and
# so there is no selection against MRs with relatively weaker correlation between their expression and activity. If this
# hypothesis is correct, I would expect to see that the MRs in common between the inference by expression and by
# activity would be those with stronger correlation (e.g. r > ~0.7) between expression and activity. Let's test this.
all_mrs=c(rownames(mrs),rownames(mrs_exp)) # 596 MRs
all_mrs=unique(all_mrs) # 461 unique MRs
all_mrs=data.frame(row.names = all_mrs,"MR_by_exp"=(all_mrs %in% rownames(mrs_exp)),
                   "MR_by_act"=(all_mrs %in% rownames(mrs)))
sum(all_mrs$MR_by_exp==all_mrs$MR_by_act) # 135/461 MRs are in common.
for(i in 1:dim(all_mrs)[1]){
  if(all_mrs$MR_by_exp[i]){
    all_mrs$exp_act_cor[i]=match_by_exp[rownames(all_mrs)[i]]
  } else {
    all_mrs$exp_act_cor[i]=match_by_vip[rownames(all_mrs)[i]]
  }
}

summary(all_mrs$exp_act_cor) # All MRs
# Min.    1st Qu.  Median  Mean    3rd Qu.  Max. 
# 0.2007  0.6694   0.7588  0.7335  0.8293   0.9287

summary(match_by_exp) # MRs by expression
# Min.    1st Qu.  Median  Mean    3rd Qu.  Max. 
# 0.2917  0.7098   0.7895  0.7663  0.8460   0.9287

summary(match_by_vip) # MRs by activity
# Min.    1st Qu.  Median  Mean    3rd Qu.  Max. 
# 0.2007  0.6392   0.7316  0.7154  0.8141   0.9287

summary(all_mrs$exp_act_cor[(all_mrs$MR_by_exp & all_mrs$MR_by_act)]) # MRs by both expression and activity
# Min.    1st Qu.  Median  Mean    3rd Qu.  Max. 
# 0.5130  0.6972   0.7833  0.7684  0.8443   0.9287

summary(all_mrs$exp_act_cor[(all_mrs$MR_by_exp & !all_mrs$MR_by_act)]) # MRs by expression only
# Min.    1st Qu.  Median  Mean    3rd Qu.  Max. 
# 0.2917  0.7127  0.7917  0.7645  0.8458  0.9285

summary(all_mrs$exp_act_cor[(all_mrs$MR_by_act & !all_mrs$MR_by_exp)]) # MRs by activity only
# Min.    1st Qu.  Median  Mean    3rd Qu.  Max. 
# 0.2007  0.5856   0.6863  0.6695  0.7754   0.9202

# Plot these distributions
par(mfrow=c(2,3))
plot(density(all_mrs$exp_act_cor),xlim=c(0,1.2),ylim=c(0,4.5),main="Exp vs Act for all matching MRs")
plot(density(match_by_exp),xlim=c(0,1.2),ylim=c(0,4.5),main="Exp vs Act for matching MRs by expression")
plot(density(match_by_vip),xlim=c(0,1.2),ylim=c(0,4.5),main="Exp vs Act for matching MRs by activity")
plot(density(all_mrs$exp_act_cor[(all_mrs$MR_by_exp & all_mrs$MR_by_act)]),xlim=c(0,1.2),ylim=c(0,4.5),main="Exp vs Act for matching MRs by expression AND activity")
plot(density(all_mrs$exp_act_cor[(all_mrs$MR_by_exp & !all_mrs$MR_by_act)]),xlim=c(0,1.2),ylim=c(0,4.5),main="Exp vs Act for matching MRs by expression only")
plot(density(all_mrs$exp_act_cor[(all_mrs$MR_by_act & !all_mrs$MR_by_exp)]),xlim=c(0,1.2),ylim=c(0,4.5),main="Exp vs Act for matching MRs by activity only")

# My hypothesis holds up to this analysis. The distribution of correlation coefficients between MRs' expression and
# activity is notably right-shifted for MRs inferred by expression relative to MRs inferred by activity (that are more
# left shifted). Another potentially interesting test of this would be to check the rank of MRs (inferred by expression) 
# based on simple linear associations between expression and phenotype. What I'm wondering is whether regulators with
# stronger linear associations with the phenotype are excluded from the RF model because they have poorer correlation
# with their corresponding activities and are consequently less robust in predicting the phenotype.

# Read in Eurobats phenotype data
pheno=read.table("./Amendment_time-matched_phenotypes_E886_02082019_with_HOMA.txt",sep="\t",header = T,row.names = 1)
whr=read.table("./Eurobats_WHR_for_subjects_with_less_than_10percent_change_in_BMI.txt",sep="\t",header = T,row.names = 1)

# Filter to overlapping samples for phenotype and expression and activity data. Note that the WHR data has fewer samples,
# so the data needs to be separately filtered.
filt_pheno=pheno[na.omit(match(colnames(vip),rownames(pheno))),]
whr_samples=rownames(whr)[na.omit(match(colnames(vip),rownames(whr)))]
filt_whr=data.frame(row.names = whr_samples,"WHR"=whr[whr_samples,])
filt_vip=vip[,rownames(filt_pheno)]
filt_vip_whr=vip[,rownames(filt_whr)]
filt_exp=exp[,colnames(filt_vip)]
filt_exp_whr=exp[,rownames(filt_whr)]
all(colnames(filt_vip)==rownames(filt_pheno)) # TRUE
all(colnames(filt_vip)==colnames(filt_exp)) # TRUE
all(colnames(filt_vip_whr)==rownames(filt_whr)) # TRUE
all(colnames(filt_vip_whr)==colnames(filt_exp_whr)) # TRUE

# Identify the regulators whose expression or activities are best associated with the phenotypes
# For expression
e_bmi=list()
for(i in 1:dim(filt_exp)[1]){
  e_bmi[[i]]=summary(lm(filt_pheno$BMI~as.numeric(filt_exp[i,])))
}
names(e_bmi)=rownames(filt_exp)
sum_e_bmi=as.data.frame(matrix(nrow = length(e_bmi),ncol = 2))
for(i in 1:length(e_bmi)){
  sum_e_bmi[i,1]=coef(e_bmi[[i]])[2,1]
  sum_e_bmi[i,2]=coef(e_bmi[[i]])[2,4]
}
colnames(sum_e_bmi)=c("Beta","P")
rownames(sum_e_bmi)=names(e_bmi)
sum_e_bmi=sum_e_bmi[order(sum_e_bmi$P),]
sum_e_bmi$BonfP=sum_e_bmi$P*dim(sum_e_bmi)[1]

e_WHR=list()
for(i in 1:dim(filt_exp_whr)[1]){
  e_WHR[[i]]=summary(lm(filt_whr$WHR~as.numeric(filt_exp_whr[i,])))
}
names(e_WHR)=rownames(filt_exp_whr)
sum_e_WHR=as.data.frame(matrix(nrow = length(e_WHR),ncol = 2))
for(i in 1:length(e_WHR)){
  sum_e_WHR[i,1]=coef(e_WHR[[i]])[2,1]
  sum_e_WHR[i,2]=coef(e_WHR[[i]])[2,4]
}
colnames(sum_e_WHR)=c("Beta","P")
rownames(sum_e_WHR)=names(e_WHR)
sum_e_WHR=sum_e_WHR[order(sum_e_WHR$P),]
sum_e_WHR$BonferroniP=sum_e_WHR$P*dim(sum_e_WHR)[1]

e_HOMAIR=list()
for(i in 1:dim(filt_exp)[1]){
  e_HOMAIR[[i]]=summary(lm(filt_pheno$HOMA.IR~as.numeric(filt_exp[i,])))
}
names(e_HOMAIR)=rownames(filt_exp)
sum_e_HOMAIR=as.data.frame(matrix(nrow = length(e_HOMAIR),ncol = 2))
for(i in 1:length(e_HOMAIR)){
  sum_e_HOMAIR[i,1]=coef(e_HOMAIR[[i]])[2,1]
  sum_e_HOMAIR[i,2]=coef(e_HOMAIR[[i]])[2,4]
}
colnames(sum_e_HOMAIR)=c("Beta","P")
rownames(sum_e_HOMAIR)=names(e_HOMAIR)
sum_e_HOMAIR=sum_e_HOMAIR[order(sum_e_HOMAIR$P),]
sum_e_HOMAIR$BonferroniP=sum_e_HOMAIR$P*dim(sum_e_HOMAIR)[1]

e_HDL=list()
for(i in 1:dim(filt_exp)[1]){
  e_HDL[[i]]=summary(lm(filt_pheno$HDLcholesterol~as.numeric(filt_exp[i,])))
}
names(e_HDL)=rownames(filt_exp)
sum_e_HDL=as.data.frame(matrix(nrow = length(e_HDL),ncol = 2))
for(i in 1:length(e_HDL)){
  sum_e_HDL[i,1]=coef(e_HDL[[i]])[2,1]
  sum_e_HDL[i,2]=coef(e_HDL[[i]])[2,4]
}
colnames(sum_e_HDL)=c("Beta","P")
rownames(sum_e_HDL)=names(e_HDL)
sum_e_HDL=sum_e_HDL[order(sum_e_HDL$P),]
sum_e_HDL$BonferroniP=sum_e_HDL$P*dim(sum_e_HDL)[1]

e_TriG=list()
for(i in 1:dim(filt_exp)[1]){
  e_TriG[[i]]=summary(lm(filt_pheno$TotalTriglycerides~as.numeric(filt_exp[i,])))
}
names(e_TriG)=rownames(filt_exp)
sum_e_TriG=as.data.frame(matrix(nrow = length(e_TriG),ncol = 2))
for(i in 1:length(e_TriG)){
  sum_e_TriG[i,1]=coef(e_TriG[[i]])[2,1]
  sum_e_TriG[i,2]=coef(e_TriG[[i]])[2,4]
}
colnames(sum_e_TriG)=c("Beta","P")
rownames(sum_e_TriG)=names(e_TriG)
sum_e_TriG=sum_e_TriG[order(sum_e_TriG$P),]
sum_e_TriG$BonferroniP=sum_e_TriG$P*dim(sum_e_TriG)[1]

# For activities
a_bmi=list()
for(i in 1:dim(filt_vip)[1]){
  a_bmi[[i]]=summary(lm(filt_pheno$BMI~as.numeric(filt_vip[i,])))
}
names(a_bmi)=rownames(filt_vip)
sum_a_bmi=as.data.frame(matrix(nrow = length(a_bmi),ncol = 2))
for(i in 1:length(a_bmi)){
  sum_a_bmi[i,1]=coef(a_bmi[[i]])[2,1]
  sum_a_bmi[i,2]=coef(a_bmi[[i]])[2,4]
}
colnames(sum_a_bmi)=c("Beta","P")
rownames(sum_a_bmi)=names(a_bmi)
sum_a_bmi=sum_a_bmi[order(sum_a_bmi$P),]
sum_a_bmi$BonfP=sum_a_bmi$P*dim(sum_a_bmi)[1]

a_WHR=list()
for(i in 1:dim(filt_vip_whr)[1]){
  a_WHR[[i]]=summary(lm(filt_whr$WHR~as.numeric(filt_vip_whr[i,])))
}
names(a_WHR)=rownames(filt_vip_whr)
sum_a_WHR=as.data.frame(matrix(nrow = length(a_WHR),ncol = 2))
for(i in 1:length(a_WHR)){
  sum_a_WHR[i,1]=coef(a_WHR[[i]])[2,1]
  sum_a_WHR[i,2]=coef(a_WHR[[i]])[2,4]
}
colnames(sum_a_WHR)=c("Beta","P")
rownames(sum_a_WHR)=names(a_WHR)
sum_a_WHR=sum_a_WHR[order(sum_a_WHR$P),]
sum_a_WHR$BonferroniP=sum_a_WHR$P*dim(sum_a_WHR)[1]

a_HOMAIR=list()
for(i in 1:dim(filt_vip)[1]){
  a_HOMAIR[[i]]=summary(lm(filt_pheno$HOMA.IR~as.numeric(filt_vip[i,])))
}
names(a_HOMAIR)=rownames(filt_vip)
sum_a_HOMAIR=as.data.frame(matrix(nrow = length(a_HOMAIR),ncol = 2))
for(i in 1:length(a_HOMAIR)){
  sum_a_HOMAIR[i,1]=coef(a_HOMAIR[[i]])[2,1]
  sum_a_HOMAIR[i,2]=coef(a_HOMAIR[[i]])[2,4]
}
colnames(sum_a_HOMAIR)=c("Beta","P")
rownames(sum_a_HOMAIR)=names(a_HOMAIR)
sum_a_HOMAIR=sum_a_HOMAIR[order(sum_a_HOMAIR$P),]
sum_a_HOMAIR$BonferroniP=sum_a_HOMAIR$P*dim(sum_a_HOMAIR)[1]

a_HDL=list()
for(i in 1:dim(filt_vip)[1]){
  a_HDL[[i]]=summary(lm(filt_pheno$HDLcholesterol~as.numeric(filt_vip[i,])))
}
names(a_HDL)=rownames(filt_vip)
sum_a_HDL=as.data.frame(matrix(nrow = length(a_HDL),ncol = 2))
for(i in 1:length(a_HDL)){
  sum_a_HDL[i,1]=coef(a_HDL[[i]])[2,1]
  sum_a_HDL[i,2]=coef(a_HDL[[i]])[2,4]
}
colnames(sum_a_HDL)=c("Beta","P")
rownames(sum_a_HDL)=names(a_HDL)
sum_a_HDL=sum_a_HDL[order(sum_a_HDL$P),]
sum_a_HDL$BonferroniP=sum_a_HDL$P*dim(sum_a_HDL)[1]

a_TriG=list()
for(i in 1:dim(filt_vip)[1]){
  a_TriG[[i]]=summary(lm(filt_pheno$TotalTriglycerides~as.numeric(filt_vip[i,])))
}
names(a_TriG)=rownames(filt_vip)
sum_a_TriG=as.data.frame(matrix(nrow = length(a_TriG),ncol = 2))
for(i in 1:length(a_TriG)){
  sum_a_TriG[i,1]=coef(a_TriG[[i]])[2,1]
  sum_a_TriG[i,2]=coef(a_TriG[[i]])[2,4]
}
colnames(sum_a_TriG)=c("Beta","P")
rownames(sum_a_TriG)=names(a_TriG)
sum_a_TriG=sum_a_TriG[order(sum_a_TriG$P),]
sum_a_TriG$BonferroniP=sum_a_TriG$P*dim(sum_a_TriG)[1]

# Now let's grab the phenotypic MRs' rankings by linear regression and compare between MRs by expression vs activity
# For MRs by expression
mrs_exp$BMI_exp_lm_rank=match(rownames(mrs_exp),rownames(sum_e_bmi))
mrs_exp$WHR_exp_lm_rank=match(rownames(mrs_exp),rownames(sum_e_WHR))
mrs_exp$HOMAIR_exp_lm_rank=match(rownames(mrs_exp),rownames(sum_e_HOMAIR))
mrs_exp$HDL_exp_lm_rank=match(rownames(mrs_exp),rownames(sum_e_HDL))
mrs_exp$TriG_exp_lm_rank=match(rownames(mrs_exp),rownames(sum_e_TriG))

# For MRs by activity
mrs$BMI_act_lm_rank=match(rownames(mrs),rownames(sum_a_bmi))
mrs$WHR_act_lm_rank=match(rownames(mrs),rownames(sum_a_WHR))
mrs$HOMAIR_act_lm_rank=match(rownames(mrs),rownames(sum_a_HOMAIR))
mrs$HDL_act_lm_rank=match(rownames(mrs),rownames(sum_a_HDL))
mrs$TriG_act_lm_rank=match(rownames(mrs),rownames(sum_a_TriG))

# Let's check the summaries
summary(match(rownames(mrs_exp[mrs_exp$BMI_MR,]),rownames(sum_e_bmi)))
summary(match(rownames(mrs[mrs$BMI_MR,]),rownames(sum_a_bmi)))
summary(match(rownames(mrs_exp[mrs_exp$WHR_MR,]),rownames(sum_e_bmi)))
summary(match(rownames(mrs[mrs$WHR_MR,]),rownames(sum_a_bmi)))
summary(match(rownames(mrs_exp[mrs_exp$HOMA.IR_MR,]),rownames(sum_e_bmi)))
summary(match(rownames(mrs[mrs$HOMA.IR_MR,]),rownames(sum_a_bmi)))
summary(match(rownames(mrs_exp[mrs_exp$HDL_MR,]),rownames(sum_e_bmi)))
summary(match(rownames(mrs[mrs$HDL_MR,]),rownames(sum_a_bmi)))
summary(match(rownames(mrs_exp[mrs_exp$TriG_MR,]),rownames(sum_e_bmi)))
summary(match(rownames(mrs[mrs$TriG_MR,]),rownames(sum_a_bmi)))
# For all but TriG, the mean rank is quite a bit lower for the MRs by activity, though that is likely driven by
# the outliers, including the max that is much higher for the MRs by expression. Looking at the medians, which
# are more robust to outliers, we see the median is lower for the MRs by activity for BMI, WHR and HOMA-IR, but
# not really different for HDL or TriG. I guess the more direct way to test this idea would be to run regression
# analyses on RF importance (especially by expression) vs exp-act correlation. The prediction would be that
# regulators for which expression is both associated with the phenotype AND a strong proxy for activity are more 
# likely to be important. However, I'm not sure I have time to deal with this. I may come back to it if I finish
# the other stuff brought up by the reviewers.






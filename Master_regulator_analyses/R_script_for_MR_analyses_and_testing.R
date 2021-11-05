### This script is for the inference and testing of the MRs for the adipose cell states corresponding to the traits tested
### (i.e. BMI, WHR, HOMA-IR, HDL and Triglycerides).
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install(c("genefilter","preprocessCore","biomaRt","DESeq2","aracne.networks","viper","vsn"))

library(viper)
library(aracne.networks)
library(digest)
library(DESeq2)
library(hexbin)
library(grid)
library(gridExtra)
library(genefilter)
library(ggplot2)
library(preprocessCore)
library(RNOmni)
library(snowfall)
library(foreach)
library(doParallel)
library(vsn)
library(randomForest)

setwd("YOUR WORKING DIRECTORY")

# This is a re-write of the PCAplot() function present within DEseq that allows you to look at more than PC1 and PC2,
# and use other expression matrices than DESeq objects. Also, you can give it a vector of values by which to shade
# the data points, which shows whether or not the covariate is evenly distributed in the PC space. 
plotPCA4 <- function (object, pcs = c(1,2), ntop = 500, returnData = FALSE, colorVar=rep(1,dim(object)[2])){
  rv <- rowVars(object)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(object[select, ]),center = T,scale. = T)
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  # Note, I do not bother renaming the PC1 and PC2 columns.
  d <- data.frame(PC1 = pca$x[, as.numeric(pcs[1])], PC2 = pca$x[, as.numeric(pcs[2])], name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2")) + 
    geom_point(size = 3,aes(col=colorVar)) + xlab(paste0("PC",pcs[1],": ", round(percentVar[pcs[1]] * 100), "% variance")) + 
    ylab(paste0("PC",pcs[2],": ", round(percentVar[pcs[2]] * 100), "% variance")) + coord_fixed() + scale_color_gradientn(colours = c("blue","red"))
}

### Run each of the new MR analyses for BMI, WHR, HOMA-IR, HDL and Triglycerides.

# Read in Eurobats data. The file locations are relative to your working directory, so adjust accordingly.
euro_pheno=read.table("./Amendment_time-matched_phenotypes_E886_02082019_with_HOMA.txt",sep="\t",header = T,row.names = 1)
euro_whr=read.table("./Eurobats_WHR_for_subjects_with_less_than_10percent_change_in_BMI.txt",sep="\t",header = T,row.names = 1)
euro_exp=read.table("./Filtered_Eurobats_adipose_qnorm_INT_logTPMs_for_4213_regulators.txt",sep = "\t",header = T,row.names = 1)
euro_vip=read.table("./Filtered_Eurobats_adipose_unnormalized_activities_from_logTPM_for_4213_regulators.txt",sep = "\t",header = T,row.names = 1)

# BMI
# Filter to overlapping samples for BMI and VIPER data
filt_bmi=euro_pheno[na.omit(match(colnames(euro_vip),rownames(euro_pheno))),]
filt_vip=euro_vip[,rownames(filt_bmi)]
filt_exp=euro_exp[,colnames(filt_vip)]
all(colnames(filt_vip)==rownames(filt_bmi)) # TRUE
all(colnames(filt_vip)==colnames(filt_exp)) # TRUE

# Split data 70/30 for training/test sets.
RNGkind(sample.kind = "Rounding") # This is now necessary after R v3.6.0 for consistent results with pre-3.6.0 scripts because the default sampler changed.
set.seed(123)
rnd_samples=sample(colnames(filt_vip),size = 490)
train_exp=filt_exp[,rnd_samples]
train_vip=filt_vip[,rnd_samples]
train_bmi=filt_bmi[rnd_samples,]
test_exp=filt_exp[,!(colnames(filt_exp) %in% rnd_samples)]
test_vip=filt_vip[,!(colnames(filt_vip) %in% rnd_samples)]
test_bmi=filt_bmi[!(rownames(filt_bmi) %in% rnd_samples),]

# Identify the regulators whose expression or activities are best associated with BMI
# For expression
e_bmi=list()
for(i in 1:dim(train_exp)[1]){
  e_bmi[[i]]=summary(lm(train_bmi$BMI~as.numeric(train_exp[i,])))
}
names(e_bmi)=rownames(train_exp)
sum_e_bmi=as.data.frame(matrix(nrow = length(e_bmi),ncol = 2))
for(i in 1:length(e_bmi)){
  sum_e_bmi[i,1]=coef(e_bmi[[i]])[2,1]
  sum_e_bmi[i,2]=coef(e_bmi[[i]])[2,4]
}
colnames(sum_e_bmi)=c("Beta","P")
rownames(sum_e_bmi)=names(e_bmi)
sum_e_bmi=sum_e_bmi[order(sum_e_bmi$P),]
sum_e_bmi$BonfP=sum_e_bmi$P*dim(sum_e_bmi)[1]
sum(sum_e_bmi$BonfP<0.05) # 1673 regulators w/ Bonferroni adjusted P < 0.05

# For activities
a_bmi=list()
for(i in 1:dim(train_vip)[1]){
  a_bmi[[i]]=summary(lm(train_bmi$BMI~as.numeric(train_vip[i,])))
}
names(a_bmi)=rownames(train_vip)
sum_a_bmi=as.data.frame(matrix(nrow = length(a_bmi),ncol = 2))
for(i in 1:length(a_bmi)){
  sum_a_bmi[i,1]=coef(a_bmi[[i]])[2,1]
  sum_a_bmi[i,2]=coef(a_bmi[[i]])[2,4]
}
colnames(sum_a_bmi)=c("Beta","P")
rownames(sum_a_bmi)=names(a_bmi)
sum_a_bmi=sum_a_bmi[order(sum_a_bmi$P),]
sum_a_bmi$BonfP=sum_a_bmi$P*dim(sum_a_bmi)[1]
sum(sum_a_bmi$BonfP<0.05) # 2030 regulators w/ Bonferroni adjusted P < 0.05

# How many regulators are significantly associated with BMI at both the expression and activity levels?
sum(rownames(sum_e_bmi)[sum_e_bmi$BonfP<0.05] %in% rownames(sum_a_bmi)[sum_a_bmi$BonfP<0.05]) 
# 1408 (84% and 69% of regulators significant by expression and activity, respectively)

# Let's use the top 500 to test in our random forest, which correspond to thresholds of P = 7.65e-15 and P = 7.8e-23 for 
# expression and activities, respectively. Note: I tried running the rfcv analysis with all 2030 Bonf sig genes, but it took 
# WAY longer and the curve didn't plateau until nearly 500 regulators were used as features. This slower plateauing is probably 
# a result of the lower frequency for each regulator being picked during the bagging. Therefore, it takes more features to 
# increase the likelihood of including important regulators.
train_exp_mrs=train_exp[rownames(sum_e_bmi)[1:500],]
test_exp_mrs=test_exp[rownames(sum_e_bmi)[1:500],]
train_vip_mrs=train_vip[rownames(sum_a_bmi)[1:500],]
test_vip_mrs=test_vip[rownames(sum_a_bmi)[1:500],]

# Let's Z-transform the expression values and activities prior to RF modeling
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
train_z=list()
train_z=foreach(i=1:dim(train_exp_mrs)[1]) %dopar%
  as.numeric(scale(as.numeric(train_exp_mrs[i,])))
train_z=as.data.frame(t(structure(train_z, row.names = c(NA, -length(train_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(train_z)=colnames(train_exp_mrs)
rownames(train_z)=rownames(train_exp_mrs)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
test_z=list()
test_z=foreach(i=1:dim(test_exp_mrs)[1]) %dopar%
  as.numeric(scale(as.numeric(test_exp_mrs[i,])))
test_z=as.data.frame(t(structure(test_z, row.names = c(NA, -length(test_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(test_z)=colnames(test_exp_mrs)
rownames(test_z)=rownames(test_exp_mrs)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
train_zvip=list()
train_zvip=foreach(i=1:dim(train_vip_mrs)[1]) %dopar%
  as.numeric(scale(as.numeric(train_vip_mrs[i,])))
train_zvip=as.data.frame(t(structure(train_zvip, row.names = c(NA, -length(train_zvip[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(train_zvip)=colnames(train_vip_mrs)
rownames(train_zvip)=rownames(train_vip_mrs)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
test_zvip=list()
test_zvip=foreach(i=1:dim(test_vip_mrs)[1]) %dopar%
  as.numeric(scale(as.numeric(test_vip_mrs[i,])))
test_zvip=as.data.frame(t(structure(test_zvip, row.names = c(NA, -length(test_zvip[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(test_zvip)=colnames(test_vip_mrs)
rownames(test_zvip)=rownames(test_vip_mrs)

# Let's run cross-validation random forest to identify the number of features required to plateau the error when inputting 
# 500 initial candidate BMI MRs 
# Based on expression
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
cv_bmi_exp_rndFor=list()
cv_bmi_exp_rndFor=foreach(i=1:12,.packages = 'randomForest') %dopar% {
  set.seed(i)
  cv_bmi_exp_rndFor[i]=rfcv(trainx=t(train_z),trainy=train_bmi$BMI,scale = F,step = -5,cv.fold = 5) # Test how many features are required before MSE plateaus
}
stopCluster(cl)

par(mfrow=c(3,4))
with(cv_bmi_exp_rndFor[[1]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~50 features is the ideal balance between parsimony and accuracy
with(cv_bmi_exp_rndFor[[2]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~75 features is the ideal balance between parsimony and accuracy
with(cv_bmi_exp_rndFor[[3]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~50 features is the ideal balance between parsimony and accuracy
with(cv_bmi_exp_rndFor[[4]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~50 features is the ideal balance between parsimony and accuracy
with(cv_bmi_exp_rndFor[[5]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~45 features is the ideal balance between parsimony and accuracy
with(cv_bmi_exp_rndFor[[6]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~35 features is the ideal balance between parsimony and accuracy
with(cv_bmi_exp_rndFor[[7]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~30 features is the ideal balance between parsimony and accuracy
with(cv_bmi_exp_rndFor[[8]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~80 features is the ideal balance between parsimony and accuracy
with(cv_bmi_exp_rndFor[[9]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~85 features is the ideal balance between parsimony and accuracy
with(cv_bmi_exp_rndFor[[10]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~40 features is the ideal balance between parsimony and accuracy
with(cv_bmi_exp_rndFor[[11]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~70 features is the ideal balance between parsimony and accuracy
with(cv_bmi_exp_rndFor[[12]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~125 features is the ideal balance between parsimony and accuracy
# Let's average across all 12 different seeds for the 500 predictor rfcv() analysis with z-scaled expression
cv_bmi_exp_avg=rowMeans(cbind(cv_bmi_exp_rndFor[[1]]$error.cv,cv_bmi_exp_rndFor[[2]]$error.cv,cv_bmi_exp_rndFor[[3]]$error.cv,cv_bmi_exp_rndFor[[4]]$error.cv,
                          cv_bmi_exp_rndFor[[5]]$error.cv,cv_bmi_exp_rndFor[[6]]$error.cv,cv_bmi_exp_rndFor[[7]]$error.cv,cv_bmi_exp_rndFor[[8]]$error.cv,
                          cv_bmi_exp_rndFor[[9]]$error.cv,cv_bmi_exp_rndFor[[10]]$error.cv,cv_bmi_exp_rndFor[[11]]$error.cv,cv_bmi_exp_rndFor[[12]]$error.cv))
cv_bmi_exp_sd=rowSds(cbind(cv_bmi_exp_rndFor[[1]]$error.cv,cv_bmi_exp_rndFor[[2]]$error.cv,cv_bmi_exp_rndFor[[3]]$error.cv,cv_bmi_exp_rndFor[[4]]$error.cv,
                       cv_bmi_exp_rndFor[[5]]$error.cv,cv_bmi_exp_rndFor[[6]]$error.cv,cv_bmi_exp_rndFor[[7]]$error.cv,cv_bmi_exp_rndFor[[8]]$error.cv,
                       cv_bmi_exp_rndFor[[9]]$error.cv,cv_bmi_exp_rndFor[[10]]$error.cv,cv_bmi_exp_rndFor[[11]]$error.cv,cv_bmi_exp_rndFor[[12]]$error.cv))
pdf("./BMI/Eurobats_adipose_CVRF_analysis_of_BMI_regulators_by_expression.pdf",width = 9)
par(mfrow=c(1,1))
plot(x=as.numeric(names(cv_bmi_exp_avg)),y=cv_bmi_exp_avg,ylim=range(c(cv_bmi_exp_avg-cv_bmi_exp_sd,cv_bmi_exp_avg+cv_bmi_exp_sd)),pch=19,xlab="Feature count",ylab="Mean cross-validation error",main="Average expression rfcv() results (12 seeds)")
arrows(x0=as.numeric(names(cv_bmi_exp_avg)),y0=cv_bmi_exp_avg-cv_bmi_exp_sd,x1=as.numeric(names(cv_bmi_exp_avg)),y1=cv_bmi_exp_avg+cv_bmi_exp_sd,length = 0.05, angle = 90, code = 3) # This trick draws the error bars
dev.off()
# This analysis points to ~65 features being sufficient when based on expression. However, for the sake of comparing to the
# activity-based MRs I will use 100.

# Based on activity
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
cv_bmi_vip_rndFor=list()
cv_bmi_vip_rndFor=foreach(i=1:12,.packages = 'randomForest') %dopar% {
  set.seed(i)
  cv_bmi_vip_rndFor[i]=rfcv(trainx=t(train_zvip),trainy=train_bmi$BMI,scale = F,step = -5,cv.fold = 5) # Test how many features are required before MSE plateaus
}
stopCluster(cl)

par(mfrow=c(3,4))
with(cv_bmi_vip_rndFor[[1]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~20 features is the ideal balance between parsimony and accuracy/purity
with(cv_bmi_vip_rndFor[[2]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~145 features is the ideal balance between parsimony and accuracy/purity
with(cv_bmi_vip_rndFor[[3]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~200 features is the ideal balance between parsimony and accuracy/purity
with(cv_bmi_vip_rndFor[[4]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~100 features is the ideal balance between parsimony and accuracy/purity
with(cv_bmi_vip_rndFor[[5]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~170 features is the ideal balance between parsimony and accuracy/purity
with(cv_bmi_vip_rndFor[[6]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~25 features is the ideal balance between parsimony and accuracy/purity
with(cv_bmi_vip_rndFor[[7]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~65 features is the ideal balance between parsimony and accuracy/purity
with(cv_bmi_vip_rndFor[[8]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~85 features is the ideal balance between parsimony and accuracy/purity
with(cv_bmi_vip_rndFor[[9]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~160 features is the ideal balance between parsimony and accuracy/purity
with(cv_bmi_vip_rndFor[[10]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~230 features is the ideal balance between parsimony and accuracy/purity
with(cv_bmi_vip_rndFor[[11]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~90 features is the ideal balance between parsimony and accuracy/purity
with(cv_bmi_vip_rndFor[[12]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~95 features is the ideal balance between parsimony and accuracy/purity
# Let's average across all 12 different seeds for the 500 predictor rfcv() analysis with z-scaled activity
cv_bmi_vip_avg=rowMeans(cbind(cv_bmi_vip_rndFor[[1]]$error.cv,cv_bmi_vip_rndFor[[2]]$error.cv,cv_bmi_vip_rndFor[[3]]$error.cv,cv_bmi_vip_rndFor[[4]]$error.cv,
                              cv_bmi_vip_rndFor[[5]]$error.cv,cv_bmi_vip_rndFor[[6]]$error.cv,cv_bmi_vip_rndFor[[7]]$error.cv,cv_bmi_vip_rndFor[[8]]$error.cv,
                              cv_bmi_vip_rndFor[[9]]$error.cv,cv_bmi_vip_rndFor[[10]]$error.cv,cv_bmi_vip_rndFor[[11]]$error.cv,cv_bmi_vip_rndFor[[12]]$error.cv))
cv_bmi_vip_sd=rowSds(cbind(cv_bmi_vip_rndFor[[1]]$error.cv,cv_bmi_vip_rndFor[[2]]$error.cv,cv_bmi_vip_rndFor[[3]]$error.cv,cv_bmi_vip_rndFor[[4]]$error.cv,
                           cv_bmi_vip_rndFor[[5]]$error.cv,cv_bmi_vip_rndFor[[6]]$error.cv,cv_bmi_vip_rndFor[[7]]$error.cv,cv_bmi_vip_rndFor[[8]]$error.cv,
                           cv_bmi_vip_rndFor[[9]]$error.cv,cv_bmi_vip_rndFor[[10]]$error.cv,cv_bmi_vip_rndFor[[11]]$error.cv,cv_bmi_vip_rndFor[[12]]$error.cv))
pdf("./BMI/Eurobats_adipose_CVRF_analysis_of_BMI_regulators_by_activity.pdf",width = 9)
par(mfrow=c(1,1))
plot(x=as.numeric(names(cv_bmi_vip_avg)),y=cv_bmi_vip_avg,ylim=range(c(cv_bmi_vip_avg-cv_bmi_vip_sd,cv_bmi_vip_avg+cv_bmi_vip_sd)),pch=19,xlab="Feature count",ylab="Mean cross-validation error",main="Average activity rfcv() results (12 seeds)")
arrows(x0=as.numeric(names(cv_bmi_vip_avg)),y0=cv_bmi_vip_avg-cv_bmi_vip_sd,x1=as.numeric(names(cv_bmi_vip_avg)),y1=cv_bmi_vip_avg+cv_bmi_vip_sd,length = 0.05, angle = 90, code = 3) # This trick draws the error bars
dev.off()
# This analysis points to ~100 features being ideal.

# Now let's look at the importance ranking among the 500 regulators in a forest based on expression, using %IncMSE to choose the top 100.
set.seed(135)
bmi_exp_rndFor500=randomForest(x=t(train_z),y=train_bmi$BMI,ntree = 1000,keep.forest = T,importance = T)
import_bmi_exp=importance(bmi_exp_rndFor500,type = 1,scale = F)
import_bmi_exp=import_bmi_exp[order(import_bmi_exp,decreasing = T),,drop=F]
bmi_mrs_exp=rownames(import_bmi_exp)[1:100]

# Now let's look at the importance ranking among the 500 regulators in a forest based on activities, using %IncMSE to choose the top 100.
set.seed(135)
bmi_vip_rndFor500=randomForest(x=t(train_zvip),y=train_bmi$BMI,ntree = 1000,keep.forest = T,importance = T)
import_bmi_vip=importance(bmi_vip_rndFor500,type = 1,scale = F)
import_bmi_vip=import_bmi_vip[order(import_bmi_vip,decreasing = T),,drop=F]
bmi_mrs_vip=rownames(import_bmi_vip)[1:100]

# How many BMI MRs are in common between the expression and activity based RFs?
sum(bmi_mrs_exp %in% bmi_mrs_vip) # 34 putative BMI MRs are identical regardless of using expression or activities

# Let's train and test the final Z-scaled expression RF model
train_bmi_exp_mrs=train_zvip[na.omit(match(bmi_mrs_exp,rownames(train_z))),]
test_bmi_exp_mrs=test_zvip[na.omit(match(bmi_mrs_exp,rownames(test_z))),]
set.seed(314)
final_bmi_exp_rndFor=randomForest(x=t(train_bmi_exp_mrs),y=train_bmi$BMI,xtest = t(test_bmi_exp_mrs),ytest = test_bmi$BMI,ntree = 1000,keep.forest = T,importance = T)
train_lm_bmi_exp=lm(train_bmi$BMI~final_bmi_exp_rndFor$predicted) # Beta=0.97514, P=8.14e-73, r-squared=0.4863
test_lm_bmi_exp=lm(test_bmi$BMI~final_bmi_exp_rndFor$test$predicted) # Beta=1.04005, P=8.90e-39, r-squared=0.558

# Let's train and test the final Z-VIPER RF model
train_bmi_vip_mrs=train_zvip[na.omit(match(bmi_mrs_vip,rownames(train_zvip))),]
test_bmi_vip_mrs=test_zvip[na.omit(match(bmi_mrs_vip,rownames(test_zvip))),]
set.seed(314)
final_bmi_vip_rndFor=randomForest(x=t(train_bmi_vip_mrs),y=train_bmi$BMI,xtest = t(test_bmi_vip_mrs),ytest = test_bmi$BMI,ntree = 1000,keep.forest = T,importance = T)
train_lm_bmi_vip=lm(train_bmi$BMI~final_bmi_vip_rndFor$predicted) # Beta=1.00234, P=2.51e-79, r-squared=0.5168
test_lm_bmi_vip=lm(test_bmi$BMI~final_bmi_vip_rndFor$test$predicted) # Beta=1.02113, P=1.12e-39, r-squared=0.5668
# Based on the TwinsUK data alone, there is not much difference in the quality of the BMI predictions based on expression
# versus activities MR RF models.

# Let's make some hi-res plots
pdf("./BMI/Eurobats_adipose_BMI_MR_by_expression_importance_plots.pdf",width = 6,height = 10)
varImpPlot(final_bmi_exp_rndFor,type=1,scale = F,n.var = 100,cex = 0.5)
dev.off()

pdf("./BMI/Eurobats_adipose_BMI_MR_by_activities_importance_plots.pdf",width = 6,height = 10)
varImpPlot(final_bmi_vip_rndFor,type=1,scale = F,n.var = 100,cex = 0.5)
dev.off()

### WHR
# Filter to overlapping samples for WHR and VIPER data
filt_whr=data.frame("WHR"=euro_whr[na.omit(match(colnames(euro_vip),rownames(euro_whr))),],row.names = rownames(euro_whr)[na.omit(match(colnames(euro_vip),rownames(euro_whr)))])
filt_vip=euro_vip[,rownames(filt_whr)]
filt_exp=euro_exp[,colnames(filt_vip)]
all(colnames(filt_vip)==rownames(filt_whr)) # TRUE
all(colnames(filt_vip)==colnames(filt_exp)) # TRUE

# Split data 70/30 for training/test sets.
RNGkind(sample.kind = "default") # Note that the original WHR RF MR analysis was performed in R >3.6.0, so it used the new default sampling.
set.seed(123)
rnd_samples=sample(colnames(filt_vip),size = 271) # 70% of 388 samples
train_exp=filt_exp[,rnd_samples]
train_vip=filt_vip[,rnd_samples]
train_whr=filt_whr[rnd_samples,"WHR"]
test_exp=filt_exp[,!(colnames(filt_exp) %in% rnd_samples)]
test_vip=filt_vip[,!(colnames(filt_vip) %in% rnd_samples)]
test_whr=filt_whr[!(rownames(filt_whr) %in% rnd_samples),"WHR"]

# Identify the regulators whose expression or activities are best associated with WHR
# For expression
e_WHR=list()
for(i in 1:dim(train_exp)[1]){
  e_WHR[[i]]=summary(lm(train_whr~as.numeric(train_exp[i,])))
}
names(e_WHR)=rownames(train_exp)
sum_e_WHR=as.data.frame(matrix(nrow = length(e_WHR),ncol = 2))
for(i in 1:length(e_WHR)){
  sum_e_WHR[i,1]=coef(e_WHR[[i]])[2,1]
  sum_e_WHR[i,2]=coef(e_WHR[[i]])[2,4]
}
colnames(sum_e_WHR)=c("Beta","P")
rownames(sum_e_WHR)=names(e_WHR)
sum_e_WHR=sum_e_WHR[order(sum_e_WHR$P),]
sum_e_WHR$BonferroniP=sum_e_WHR$P*dim(sum_e_WHR)[1]
sum(sum_e_WHR$BonferroniP<0.05) # 609 regulators w/ Bonferroni adjusted P < 0.05
sum(sum_e_WHR$P<0.05) # 2196 regulators w/ P < 0.05

# For activities
a_WHR=list()
for(i in 1:dim(train_vip)[1]){
  a_WHR[[i]]=summary(lm(train_whr~as.numeric(train_vip[i,])))
}
names(a_WHR)=rownames(train_vip)
sum_a_WHR=as.data.frame(matrix(nrow = length(a_WHR),ncol = 2))
for(i in 1:length(a_WHR)){
  sum_a_WHR[i,1]=coef(a_WHR[[i]])[2,1]
  sum_a_WHR[i,2]=coef(a_WHR[[i]])[2,4]
}
colnames(sum_a_WHR)=c("Beta","P")
rownames(sum_a_WHR)=names(a_WHR)
sum_a_WHR=sum_a_WHR[order(sum_a_WHR$P),]
sum_a_WHR$BonferroniP=sum_a_WHR$P*dim(sum_a_WHR)[1]
sum(sum_a_WHR$BonferroniP<0.05) # 1130 regulators w/ Bonferroni adjusted P < 0.05
sum(sum_a_WHR$P<0.05) # 2516 regulators w/ P < 0.05

# Let's use the top 500 to train our random forest
train_exp_top=train_exp[rownames(sum_e_WHR)[1:500],]
test_exp_top=test_exp[rownames(sum_e_WHR)[1:500],]
train_vip_top=train_vip[rownames(sum_a_WHR)[1:500],]
test_vip_top=test_vip[rownames(sum_a_WHR)[1:500],]

# Let's Z-transform the activities prior to RF modeling
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
train_z=list()
train_z=foreach(i=1:dim(train_exp_top)[1]) %dopar%
  as.numeric(scale(as.numeric(train_exp_top[i,])))
train_z=as.data.frame(t(structure(train_z, row.names = c(NA, -length(train_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(train_z)=colnames(train_exp_top)
rownames(train_z)=rownames(train_exp_top)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
test_z=list()
test_z=foreach(i=1:dim(test_exp_top)[1]) %dopar%
  as.numeric(scale(as.numeric(test_exp_top[i,])))
test_z=as.data.frame(t(structure(test_z, row.names = c(NA, -length(test_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(test_z)=colnames(test_exp_top)
rownames(test_z)=rownames(test_exp_top)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
train_zvip=list()
train_zvip=foreach(i=1:dim(train_vip_top)[1]) %dopar%
  as.numeric(scale(as.numeric(train_vip_top[i,])))
train_zvip=as.data.frame(t(structure(train_zvip, row.names = c(NA, -length(train_zvip[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(train_zvip)=colnames(train_vip_top)
rownames(train_zvip)=rownames(train_vip_top)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
test_zvip=list()
test_zvip=foreach(i=1:dim(test_vip_top)[1]) %dopar%
  as.numeric(scale(as.numeric(test_vip_top[i,])))
test_zvip=as.data.frame(t(structure(test_zvip, row.names = c(NA, -length(test_zvip[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(test_zvip)=colnames(test_vip_top)
rownames(test_zvip)=rownames(test_vip_top)

# Let's run cross-validation random forest to identify the number of features required to plateau the error when inputing 500 candidate WHR MRs
# For expression
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
cv_whr_exp_rndFor=list()
cv_whr_exp_rndFor=foreach(i=1:12,.packages = 'randomForest') %dopar% {
  set.seed(i)
  cv_whr_exp_rndFor[i]=rfcv(trainx=t(train_z),trainy=train_whr,scale = F,step = -5,cv.fold = 5) # Test how many features are required before MSE plateaus
}
stopCluster(cl)

par(mfrow=c(3,4),mar=c(5,4,4,2)+0.1)
with(cv_whr_exp_rndFor[[1]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~95 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_exp_rndFor[[2]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~100 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_exp_rndFor[[3]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~90 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_exp_rndFor[[4]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~220 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_exp_rndFor[[5]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~105 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_exp_rndFor[[6]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~300 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_exp_rndFor[[7]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~120 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_exp_rndFor[[8]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~95 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_exp_rndFor[[9]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~300 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_exp_rndFor[[10]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~60 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_exp_rndFor[[11]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~300 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_exp_rndFor[[12]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~50 features is the ideal balance between parsimony and accuracy/purity
# Let's average across all 12 different seeds for the 500 predictor rfcv() analysis with zvip
cv_whr_exp_avg=rowMeans(cbind(cv_whr_exp_rndFor[[1]]$error.cv,cv_whr_exp_rndFor[[2]]$error.cv,cv_whr_exp_rndFor[[3]]$error.cv,cv_whr_exp_rndFor[[4]]$error.cv,
                          cv_whr_exp_rndFor[[5]]$error.cv,cv_whr_exp_rndFor[[6]]$error.cv,cv_whr_exp_rndFor[[7]]$error.cv,cv_whr_exp_rndFor[[8]]$error.cv,
                          cv_whr_exp_rndFor[[9]]$error.cv,cv_whr_exp_rndFor[[10]]$error.cv,cv_whr_exp_rndFor[[11]]$error.cv,cv_whr_exp_rndFor[[12]]$error.cv))
cv_whr_exp_sd=rowSds(cbind(cv_whr_exp_rndFor[[1]]$error.cv,cv_whr_exp_rndFor[[2]]$error.cv,cv_whr_exp_rndFor[[3]]$error.cv,cv_whr_exp_rndFor[[4]]$error.cv,
                       cv_whr_exp_rndFor[[5]]$error.cv,cv_whr_exp_rndFor[[6]]$error.cv,cv_whr_exp_rndFor[[7]]$error.cv,cv_whr_exp_rndFor[[8]]$error.cv,
                       cv_whr_exp_rndFor[[9]]$error.cv,cv_whr_exp_rndFor[[10]]$error.cv,cv_whr_exp_rndFor[[11]]$error.cv,cv_whr_exp_rndFor[[12]]$error.cv))
pdf("./WHR/Eurobats_adipose_CVRF_analysis_of_WHR_regulators_by_expression.pdf",width = 9)
par(mfrow=c(1,1))
plot(x=as.numeric(names(cv_whr_exp_avg)),y=cv_whr_exp_avg,ylim=range(c(cv_whr_exp_avg-cv_whr_exp_sd,cv_whr_exp_avg+cv_whr_exp_sd)),pch=19,xlab="Feature count",ylab="Mean cross-validation error",main="Average expression rfcv() results (12 seeds)")
arrows(x0=as.numeric(names(cv_whr_exp_avg)),y0=cv_whr_exp_avg-cv_whr_exp_sd,x1=as.numeric(names(cv_whr_exp_avg)),y1=cv_whr_exp_avg+cv_whr_exp_sd,length = 0.05, angle = 90, code = 3) # This trick draws the error bars
dev.off()
# Let's go with ~100 MRs based on expression.

# For activity
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
cv_whr_vip_rndFor=list()
cv_whr_vip_rndFor=foreach(i=1:12,.packages = 'randomForest') %dopar% {
  set.seed(i)
  cv_whr_vip_rndFor[i]=rfcv(trainx=t(train_zvip),trainy=train_whr,scale = F,step = -5,cv.fold = 5) # Test how many features are required before MSE plateaus
}
stopCluster(cl)

par(mfrow=c(3,4),mar=c(5,4,4,2)+0.1)
with(cv_whr_vip_rndFor[[1]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~20 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_vip_rndFor[[2]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~145 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_vip_rndFor[[3]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~200 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_vip_rndFor[[4]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~100 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_vip_rndFor[[5]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~170 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_vip_rndFor[[6]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~25 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_vip_rndFor[[7]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~65 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_vip_rndFor[[8]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~85 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_vip_rndFor[[9]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~160 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_vip_rndFor[[10]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~230 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_vip_rndFor[[11]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~90 features is the ideal balance between parsimony and accuracy/purity
with(cv_whr_vip_rndFor[[12]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~95 features is the ideal balance between parsimony and accuracy/purity
# Let's average across all 12 different seeds for the 500 predictor rfcv() analysis with zvip
cv_whr_vip_avg=rowMeans(cbind(cv_whr_vip_rndFor[[1]]$error.cv,cv_whr_vip_rndFor[[2]]$error.cv,cv_whr_vip_rndFor[[3]]$error.cv,cv_whr_vip_rndFor[[4]]$error.cv,
                              cv_whr_vip_rndFor[[5]]$error.cv,cv_whr_vip_rndFor[[6]]$error.cv,cv_whr_vip_rndFor[[7]]$error.cv,cv_whr_vip_rndFor[[8]]$error.cv,
                              cv_whr_vip_rndFor[[9]]$error.cv,cv_whr_vip_rndFor[[10]]$error.cv,cv_whr_vip_rndFor[[11]]$error.cv,cv_whr_vip_rndFor[[12]]$error.cv))
cv_whr_vip_sd=rowSds(cbind(cv_whr_vip_rndFor[[1]]$error.cv,cv_whr_vip_rndFor[[2]]$error.cv,cv_whr_vip_rndFor[[3]]$error.cv,cv_whr_vip_rndFor[[4]]$error.cv,
                           cv_whr_vip_rndFor[[5]]$error.cv,cv_whr_vip_rndFor[[6]]$error.cv,cv_whr_vip_rndFor[[7]]$error.cv,cv_whr_vip_rndFor[[8]]$error.cv,
                           cv_whr_vip_rndFor[[9]]$error.cv,cv_whr_vip_rndFor[[10]]$error.cv,cv_whr_vip_rndFor[[11]]$error.cv,cv_whr_vip_rndFor[[12]]$error.cv))
pdf("./WHR/Eurobats_adipose_CVRF_analysis_of_WHR_regulators_by_activity.pdf",width = 9)
par(mfrow=c(1,1))
plot(x=as.numeric(names(cv_whr_vip_avg)),y=cv_whr_vip_avg,ylim=range(c(cv_whr_vip_avg-cv_whr_vip_sd,cv_whr_vip_avg+cv_whr_vip_sd)),pch=19,xlab="Feature count",ylab="Mean cross-validation error",main="Average activity rfcv() results (12 seeds)")
arrows(x0=as.numeric(names(cv_whr_vip_avg)),y0=cv_whr_vip_avg-cv_whr_vip_sd,x1=as.numeric(names(cv_whr_vip_avg)),y1=cv_whr_vip_avg+cv_whr_vip_sd,length = 0.05, angle = 90, code = 3) # This trick draws the error bars
dev.off()
# Let's go with 100 MRs as before.

# Let's train a random forest with the 500 regulators to pick the top 100 (by importance type 1) as WHR MRs
# For expression
set.seed(246)
whr_exp_rndFor500=randomForest(x=t(train_exp_top),y=train_whr,ntree = 1000,keep.forest = T,importance = T)
import_whr500_exp=importance(whr_exp_rndFor500,type = 1,scale = F)
import_whr500_exp=import_whr500_exp[order(import_whr500_exp,decreasing = T),,drop=F]
train_whr_exp_mrs=train_exp_top[rownames(import_whr500_exp)[1:100],]
test_whr_exp_mrs=test_exp_top[rownames(import_whr500_exp)[1:100],]

# For activity
set.seed(246)
whr_vip_rndFor500=randomForest(x=t(train_vip_top),y=train_whr,ntree = 1000,keep.forest = T,importance = T)
import_whr500_vip=importance(whr_vip_rndFor500,type = 1,scale = F)
import_whr500_vip=import_whr500_vip[order(import_whr500_vip,decreasing = T),,drop=F]
train_whr_vip_mrs=train_vip_top[rownames(import_whr500_vip)[1:100],]
test_whr_vip_mrs=test_vip_top[rownames(import_whr500_vip)[1:100],]

# How many WHR MRs are in common between the expression and activity based RFs?
sum(rownames(train_whr_exp_mrs) %in% rownames(train_whr_vip_mrs)) # 28 putative WHR MRs are identical regardless of using expression or activities

# Let's train and test the final WHR random forest models
# For expression
set.seed(314)
final_whr_exp_rndFor=randomForest(x=t(train_whr_exp_mrs),y=train_whr,xtest = t(test_whr_exp_mrs),ytest = test_whr,ntree = 1000,keep.forest = T,importance = T)
train_lm_whr_exp=lm(train_whr~final_whr_exp_rndFor$predicted) # Beta=1.03694, P=8.96e-27, r-squared=0.3453
test_lm_whr_exp=lm(test_whr~final_whr_exp_rndFor$test$predicted) # Beta=0.94845, P=3.44e-8, r-squared=0.2268

# For activity
set.seed(314)
final_whr_vip_rndFor=randomForest(x=t(train_whr_vip_mrs),y=train_whr,xtest = t(test_whr_vip_mrs),ytest = test_whr,ntree = 1000,keep.forest = T,importance = T)
train_lm_whr_vip=lm(train_whr~final_whr_vip_rndFor$predicted) # Beta=0.9257, P=7.21e-23, r-squared=0.3004
test_lm_whr_vip=lm(test_whr~final_whr_vip_rndFor$test$predicted) # Beta=0.8991, P=4.79e-8, r-squared=0.2224
# Based on the TwinsUK data alone, the WHR predictions based on expression versus activities MR RF models are
# better in the training data but very comparable in the test data.

# Let's make some hi-res plots
pdf("./WHR/Eurobats_adipose_WHR_MR_by_expression_importance_plots.pdf",width = 6,height = 10)
varImpPlot(final_whr_exp_rndFor,type=1,scale = F,n.var = 100,cex = 0.5)
dev.off()

pdf("./WHR/Eurobats_adipose_WHR_MR_by_activities_importance_plots.pdf",width = 6,height = 10)
varImpPlot(final_whr_vip_rndFor,type=1,scale = F,n.var = 100,cex = 0.5)
dev.off()

### HOMA-IR
# Filter to overlapping samples for HOMA-IR and VIPER data
filt_phenos=euro_pheno[na.omit(match(colnames(euro_vip),rownames(euro_pheno))),]
filt_phenos=filt_phenos[!is.na(filt_phenos$HOMA.IR),] # 659 samples
filt_vip=euro_vip[,rownames(filt_phenos)]
filt_exp=euro_exp[,colnames(filt_vip)]
all(colnames(filt_vip)==rownames(filt_phenos)) # TRUE
all(colnames(filt_vip)==colnames(filt_exp)) # TRUE

# HOMA-IR was ln transformed in the GWAS paper, so I'm doing it here too.
filt_phenos$HOMA.IR=log(filt_phenos$HOMA.IR)

# Split data 70/30 for training/test sets. 
RNGkind(sample.kind = "default") # Note that the original WHR RF MR analysis was performed in R >3.6.0, so it used the new default sampling.
set.seed(123)
rnd_samples=sample(colnames(filt_vip),size = 461) # 70% of 659 samples
train_exp=filt_exp[,rnd_samples]
train_vip=filt_vip[,rnd_samples]
train_HOMAIR=filt_phenos[rnd_samples,"HOMA.IR"]
test_exp=filt_exp[,!(colnames(filt_exp) %in% rnd_samples)]
test_vip=filt_vip[,!(colnames(filt_vip) %in% rnd_samples)]
test_HOMAIR=filt_phenos[!(rownames(filt_phenos) %in% rnd_samples),"HOMA.IR"]

# Identify the regulators whose expression or activities are best associated with HOMA-IR
# By expression
e_HOMAIR=list()
for(i in 1:dim(train_exp)[1]){
  e_HOMAIR[[i]]=summary(lm(train_HOMAIR~as.numeric(train_exp[i,])))
}
names(e_HOMAIR)=rownames(train_exp)
sum_e_HOMAIR=as.data.frame(matrix(nrow = length(e_HOMAIR),ncol = 2))
for(i in 1:length(e_HOMAIR)){
  sum_e_HOMAIR[i,1]=coef(e_HOMAIR[[i]])[2,1]
  sum_e_HOMAIR[i,2]=coef(e_HOMAIR[[i]])[2,4]
}
colnames(sum_e_HOMAIR)=c("Beta","P")
rownames(sum_e_HOMAIR)=names(e_HOMAIR)
sum_e_HOMAIR=sum_e_HOMAIR[order(sum_e_HOMAIR$P),]
sum_e_HOMAIR$BonferroniP=sum_e_HOMAIR$P*dim(sum_e_HOMAIR)[1]
sum(sum_e_HOMAIR$BonferroniP<0.05) # 1440 regulators w/ Bonferroni adjusted P < 0.05
sum(sum_e_HOMAIR$P<0.05) # 2856 regulators w/ P < 0.05

# By activity
a_HOMAIR=list()
for(i in 1:dim(train_vip)[1]){
  a_HOMAIR[[i]]=summary(lm(train_HOMAIR~as.numeric(train_vip[i,])))
}
names(a_HOMAIR)=rownames(train_vip)
sum_a_HOMAIR=as.data.frame(matrix(nrow = length(a_HOMAIR),ncol = 2))
for(i in 1:length(a_HOMAIR)){
  sum_a_HOMAIR[i,1]=coef(a_HOMAIR[[i]])[2,1]
  sum_a_HOMAIR[i,2]=coef(a_HOMAIR[[i]])[2,4]
}
colnames(sum_a_HOMAIR)=c("Beta","P")
rownames(sum_a_HOMAIR)=names(a_HOMAIR)
sum_a_HOMAIR=sum_a_HOMAIR[order(sum_a_HOMAIR$P),]
sum_a_HOMAIR$BonferroniP=sum_a_HOMAIR$P*dim(sum_a_HOMAIR)[1]
sum(sum_a_HOMAIR$BonferroniP<0.05) # 1827 regulators w/ Bonferroni adjusted P < 0.05
sum(sum_a_HOMAIR$P<0.05) # 3093 regulators w/ P < 0.05

# Let's use the top 500 to train our random forest
train_homair_exp_top=train_exp[rownames(sum_e_HOMAIR)[1:500],]
train_homair_vip_top=train_vip[rownames(sum_a_HOMAIR)[1:500],]
test_homair_exp_top=test_exp[rownames(sum_e_HOMAIR)[1:500],]
test_homair_vip_top=test_vip[rownames(sum_a_HOMAIR)[1:500],]

# Let's Z-transform the activities prior to RF modeling
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
train_z=list()
train_z=foreach(i=1:dim(train_homair_exp_top)[1]) %dopar%
  as.numeric(scale(as.numeric(train_homair_exp_top[i,])))
train_z=as.data.frame(t(structure(train_z, row.names = c(NA, -length(train_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(train_z)=colnames(train_homair_exp_top)
rownames(train_z)=rownames(train_homair_exp_top)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
test_z=list()
test_z=foreach(i=1:dim(test_homair_exp_top)[1]) %dopar%
  as.numeric(scale(as.numeric(test_homair_exp_top[i,])))
test_z=as.data.frame(t(structure(test_z, row.names = c(NA, -length(test_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(test_z)=colnames(test_homair_exp_top)
rownames(test_z)=rownames(test_homair_exp_top)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
train_zvip=list()
train_zvip=foreach(i=1:dim(train_homair_vip_top)[1]) %dopar%
  as.numeric(scale(as.numeric(train_homair_vip_top[i,])))
train_zvip=as.data.frame(t(structure(train_zvip, row.names = c(NA, -length(train_zvip[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(train_zvip)=colnames(train_homair_vip_top)
rownames(train_zvip)=rownames(train_homair_vip_top)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
test_zvip=list()
test_zvip=foreach(i=1:dim(test_homair_vip_top)[1]) %dopar%
  as.numeric(scale(as.numeric(test_homair_vip_top[i,])))
test_zvip=as.data.frame(t(structure(test_zvip, row.names = c(NA, -length(test_zvip[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(test_zvip)=colnames(test_homair_vip_top)
rownames(test_zvip)=rownames(test_homair_vip_top)

# Let's run cross-validation random forest to identify the number of features required to plateau the error when inputing 500 candidate HOMA-IR MRs
# For expression
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
cv_homair_exp_rndFor=list()
cv_homair_exp_rndFor=foreach(i=1:12,.packages = 'randomForest') %dopar% {
  set.seed(i)
  cv_homair_exp_rndFor[i]=rfcv(trainx=t(train_z),trainy=train_HOMAIR,scale = F,step = -5,cv.fold = 5) # Test how many features are required before MSE plateaus
}
stopCluster(cl)

par(mfrow=c(3,4))
with(cv_homair_exp_rndFor[[1]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~20 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_exp_rndFor[[2]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~145 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_exp_rndFor[[3]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~200 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_exp_rndFor[[4]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~100 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_exp_rndFor[[5]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~170 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_exp_rndFor[[6]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~25 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_exp_rndFor[[7]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~65 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_exp_rndFor[[8]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~85 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_exp_rndFor[[9]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~160 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_exp_rndFor[[10]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~230 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_exp_rndFor[[11]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~90 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_exp_rndFor[[12]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~95 features is the ideal balance between parsimony and accuracy/purity
# Let's average across all 12 different seeds for the 500 predictor rfcv() analysis with zvip
cv_homair_exp_avg=rowMeans(cbind(cv_homair_exp_rndFor[[1]]$error.cv,cv_homair_exp_rndFor[[2]]$error.cv,cv_homair_exp_rndFor[[3]]$error.cv,cv_homair_exp_rndFor[[4]]$error.cv,
                                 cv_homair_exp_rndFor[[5]]$error.cv,cv_homair_exp_rndFor[[6]]$error.cv,cv_homair_exp_rndFor[[7]]$error.cv,cv_homair_exp_rndFor[[8]]$error.cv,
                                 cv_homair_exp_rndFor[[9]]$error.cv,cv_homair_exp_rndFor[[10]]$error.cv,cv_homair_exp_rndFor[[11]]$error.cv,cv_homair_exp_rndFor[[12]]$error.cv))
cv_homair_exp_sd=rowSds(cbind(cv_homair_exp_rndFor[[1]]$error.cv,cv_homair_exp_rndFor[[2]]$error.cv,cv_homair_exp_rndFor[[3]]$error.cv,cv_homair_exp_rndFor[[4]]$error.cv,
                              cv_homair_exp_rndFor[[5]]$error.cv,cv_homair_exp_rndFor[[6]]$error.cv,cv_homair_exp_rndFor[[7]]$error.cv,cv_homair_exp_rndFor[[8]]$error.cv,
                              cv_homair_exp_rndFor[[9]]$error.cv,cv_homair_exp_rndFor[[10]]$error.cv,cv_homair_exp_rndFor[[11]]$error.cv,cv_homair_exp_rndFor[[12]]$error.cv))
pdf("./HOMA-IR/Eurobats_adipose_CVRF_analysis_of_HOMA-IR_regulators_by_expression.pdf",width = 9)
par(mfrow=c(1,1))
plot(x=as.numeric(names(cv_homair_exp_avg)),y=cv_homair_exp_avg,ylim=range(c(cv_homair_exp_avg-cv_homair_exp_sd,cv_homair_exp_avg+cv_homair_exp_sd)),pch=19,xlab="Feature count",ylab="Mean cross-validation error",main="Average expression rfcv() results (12 seeds)")
arrows(x0=as.numeric(names(cv_homair_exp_avg)),y0=cv_homair_exp_avg-cv_homair_exp_sd,x1=as.numeric(names(cv_homair_exp_avg)),y1=cv_homair_exp_avg+cv_homair_exp_sd,length = 0.05, angle = 90, code = 3) # This trick draws the error bars
dev.off()
# 100 HOMA-IR MRs should work fine.

# For activity
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
cv_homair_vip_rndFor=list()
cv_homair_vip_rndFor=foreach(i=1:12,.packages = 'randomForest') %dopar% {
  set.seed(i)
  cv_homair_vip_rndFor[i]=rfcv(trainx=t(train_zvip),trainy=train_HOMAIR,scale = F,step = -5,cv.fold = 5) # Test how many features are required before MSE plateaus
}
stopCluster(cl)

par(mfrow=c(3,4))
with(cv_homair_vip_rndFor[[1]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~20 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_vip_rndFor[[2]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~145 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_vip_rndFor[[3]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~200 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_vip_rndFor[[4]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~100 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_vip_rndFor[[5]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~170 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_vip_rndFor[[6]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~25 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_vip_rndFor[[7]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~65 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_vip_rndFor[[8]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~85 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_vip_rndFor[[9]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~160 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_vip_rndFor[[10]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~230 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_vip_rndFor[[11]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~90 features is the ideal balance between parsimony and accuracy/purity
with(cv_homair_vip_rndFor[[12]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~95 features is the ideal balance between parsimony and accuracy/purity
# Let's average across all 12 different seeds for the 500 predictor rfcv() analysis with zvip
cv_homair_vip_avg=rowMeans(cbind(cv_homair_vip_rndFor[[1]]$error.cv,cv_homair_vip_rndFor[[2]]$error.cv,cv_homair_vip_rndFor[[3]]$error.cv,cv_homair_vip_rndFor[[4]]$error.cv,
                                 cv_homair_vip_rndFor[[5]]$error.cv,cv_homair_vip_rndFor[[6]]$error.cv,cv_homair_vip_rndFor[[7]]$error.cv,cv_homair_vip_rndFor[[8]]$error.cv,
                                 cv_homair_vip_rndFor[[9]]$error.cv,cv_homair_vip_rndFor[[10]]$error.cv,cv_homair_vip_rndFor[[11]]$error.cv,cv_homair_vip_rndFor[[12]]$error.cv))
cv_homair_vip_sd=rowSds(cbind(cv_homair_vip_rndFor[[1]]$error.cv,cv_homair_vip_rndFor[[2]]$error.cv,cv_homair_vip_rndFor[[3]]$error.cv,cv_homair_vip_rndFor[[4]]$error.cv,
                              cv_homair_vip_rndFor[[5]]$error.cv,cv_homair_vip_rndFor[[6]]$error.cv,cv_homair_vip_rndFor[[7]]$error.cv,cv_homair_vip_rndFor[[8]]$error.cv,
                              cv_homair_vip_rndFor[[9]]$error.cv,cv_homair_vip_rndFor[[10]]$error.cv,cv_homair_vip_rndFor[[11]]$error.cv,cv_homair_vip_rndFor[[12]]$error.cv))
pdf("./HOMA-IR/Eurobats_adipose_CVRF_analysis_of_HOMA-IR_regulators_by_activity.pdf",width = 9)
par(mfrow=c(1,1))
plot(x=as.numeric(names(cv_homair_vip_avg)),y=cv_homair_vip_avg,ylim=range(c(cv_homair_vip_avg-cv_homair_vip_sd,cv_homair_vip_avg+cv_homair_vip_sd)),pch=19,xlab="Feature count",ylab="Mean cross-validation error",main="Average activity rfcv() results (12 seeds)")
arrows(x0=as.numeric(names(cv_homair_vip_avg)),y0=cv_homair_vip_avg-cv_homair_vip_sd,x1=as.numeric(names(cv_homair_vip_avg)),y1=cv_homair_vip_avg+cv_homair_vip_sd,length = 0.05, angle = 90, code = 3) # This trick draws the error bars
dev.off()
# 100 HOMA-IR MRs should work fine.

# Let's train a random forest with the 500 regulators to pick the top 100 (by importance type 1) as HOMA-IR MRs
# For expression
set.seed(246)
homair_exp_rndFor500=randomForest(x=t(train_homair_exp_top),y=train_HOMAIR,ntree = 1000,keep.forest = T,importance = T)
par(mfrow=c(1,1))
plot(homair_exp_rndFor500)
import_homair500_exp=importance(homair_exp_rndFor500,type = 1,scale = F)
import_homair500_exp=import_homair500_exp[order(import_homair500_exp,decreasing = T),,drop=F]
train_homair_exp_mrs=train_homair_exp_top[rownames(import_homair500_exp)[1:100],]
test_homair_exp_mrs=test_homair_exp_top[rownames(import_homair500_exp)[1:100],]

# For activity
set.seed(246)
homair_vip_rndFor500=randomForest(x=t(train_homair_vip_top),y=train_HOMAIR,ntree = 1000,keep.forest = T,importance = T)
par(mfrow=c(1,1))
plot(homair_vip_rndFor500)
import_homair500_exp=importance(homair_vip_rndFor500,type = 1,scale = F)
import_homair500_exp=import_homair500_exp[order(import_homair500_exp,decreasing = T),,drop=F]
train_homair_vip_mrs=train_homair_vip_top[rownames(import_homair500_exp)[1:100],]
test_homair_vip_mrs=test_homair_vip_top[rownames(import_homair500_exp)[1:100],]

# How many HOMA-IR MRs are in common between the expression and activity based RFs?
sum(rownames(train_homair_exp_mrs) %in% rownames(train_homair_vip_mrs)) # 39 putative HOMA-IR MRs are identical regardless of using expression or activities

# Let's train and test a final HOMA-IR random forest
# For expression
set.seed(314)
final_homair_exp_rndFor=randomForest(x=t(train_homair_exp_mrs),y=train_HOMAIR,xtest = t(test_homair_exp_mrs),ytest = test_HOMAIR,ntree = 1000,keep.forest = T,importance = T)
train_lm_homair_exp=lm(train_HOMAIR~final_homair_exp_rndFor$predicted) # Beta=1.08591, P=2.40e-68, r-squared=0.4849
test_lm_homair_exp=lm(test_HOMAIR~final_homair_exp_rndFor$test$predicted) # Beta=1.157705, P=8.65e-34, r-squared=0.5256

# For activity
set.seed(314)
final_homair_vip_rndFor=randomForest(x=t(train_homair_vip_mrs),y=train_HOMAIR,xtest = t(test_homair_vip_mrs),ytest = test_HOMAIR,ntree = 1000,keep.forest = T,importance = T)
train_lm_homair_vip=lm(train_HOMAIR~final_homair_vip_rndFor$predicted) # Beta=1.0013, P=5.96e-61, r-squared=0.4454
test_lm_homair_vip=lm(test_HOMAIR~final_homair_vip_rndFor$test$predicted) # Beta=1.0793, P=5.02e-31, r-squared=0.4939

# Let's make some hi-res plots
pdf("./HOMA-IR/Eurobats_adipose_HOMA-IR_MRs_by_expression_importance_plots.pdf",width = 6,height = 10)
varImpPlot(final_homair_exp_rndFor,type=1,scale = F,n.var = 100,cex = 0.5)
dev.off()

pdf("./HOMA-IR/Eurobats_adipose_HOMA-IR_MRs_by_activity_importance_plots.pdf",width = 6,height = 10)
varImpPlot(final_homair_vip_rndFor,type=1,scale = F,n.var = 100,cex = 0.5)
dev.off()

### HDL
# Filter to overlapping samples for HDL and VIPER data
filt_phenos=euro_pheno[na.omit(match(colnames(euro_vip),rownames(euro_pheno))),]
filt_phenos=filt_phenos[!is.na(filt_phenos$HDLcholesterol),] # 698 samples
filt_vip=euro_vip[,rownames(filt_phenos)]
filt_exp=euro_exp[,colnames(filt_vip)]
all(colnames(filt_vip)==rownames(filt_phenos)) # TRUE
all(colnames(filt_vip)==colnames(filt_exp)) # TRUE

# Split data 70/30 for training/test sets. 
RNGkind(sample.kind = "default") # Note that the original HDL RF MR analysis was performed in R >3.6.0, so it used the new default sampling.
set.seed(123)
rnd_samples=sample(colnames(filt_vip),size = 488) # 70% of 697 samples
train_exp=filt_exp[,rnd_samples]
train_vip=filt_vip[,rnd_samples]
train_HDL=filt_phenos[rnd_samples,"HDLcholesterol"]
test_exp=filt_exp[,!(colnames(filt_exp) %in% rnd_samples)]
test_vip=filt_vip[,!(colnames(filt_vip) %in% rnd_samples)]
test_HDL=filt_phenos[!(rownames(filt_phenos) %in% rnd_samples),"HDLcholesterol"]

# Identify the regulators whose expression or activities are best associated with HDL
# For expression
e_HDL=list()
for(i in 1:dim(train_exp)[1]){
  e_HDL[[i]]=summary(lm(train_HDL~as.numeric(train_exp[i,])))
}
names(e_HDL)=rownames(train_exp)
sum_e_HDL=as.data.frame(matrix(nrow = length(e_HDL),ncol = 2))
for(i in 1:length(e_HDL)){
  sum_e_HDL[i,1]=coef(e_HDL[[i]])[2,1]
  sum_e_HDL[i,2]=coef(e_HDL[[i]])[2,4]
}
colnames(sum_e_HDL)=c("Beta","P")
rownames(sum_e_HDL)=names(e_HDL)
sum_e_HDL=sum_e_HDL[order(sum_e_HDL$P),]
sum_e_HDL$BonferroniP=sum_e_HDL$P*dim(sum_e_HDL)[1]
sum(sum_e_HDL$BonferroniP<0.05) # 1445 regulators w/ Bonferroni adjusted P < 0.05
sum(sum_e_HDL$P<0.05) # 2883 regulators w/ P < 0.05

# For activity
a_HDL=list()
for(i in 1:dim(train_vip)[1]){
  a_HDL[[i]]=summary(lm(train_HDL~as.numeric(train_vip[i,])))
}
names(a_HDL)=rownames(train_vip)
sum_a_HDL=as.data.frame(matrix(nrow = length(a_HDL),ncol = 2))
for(i in 1:length(a_HDL)){
  sum_a_HDL[i,1]=coef(a_HDL[[i]])[2,1]
  sum_a_HDL[i,2]=coef(a_HDL[[i]])[2,4]
}
colnames(sum_a_HDL)=c("Beta","P")
rownames(sum_a_HDL)=names(a_HDL)
sum_a_HDL=sum_a_HDL[order(sum_a_HDL$P),]
sum_a_HDL$BonferroniP=sum_a_HDL$P*dim(sum_a_HDL)[1]
sum(sum_a_HDL$BonferroniP<0.05) # 2081 regulators w/ Bonferroni adjusted P < 0.05
sum(sum_a_HDL$P<0.05) # 3303 regulators w/ P < 0.05

# Let's use the top 500 to train our random forest
train_hdl_exp_top=train_exp[rownames(sum_e_HDL)[1:500],]
train_hdl_vip_top=train_vip[rownames(sum_a_HDL)[1:500],]
test_hdl_exp_top=test_exp[rownames(sum_e_HDL)[1:500],]
test_hdl_vip_top=test_vip[rownames(sum_a_HDL)[1:500],]

# Let's Z-transform the activities prior to RF modeling
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
train_z=list()
train_z=foreach(i=1:dim(train_hdl_exp_top)[1]) %dopar%
  as.numeric(scale(as.numeric(train_hdl_exp_top[i,])))
train_z=as.data.frame(t(structure(train_z, row.names = c(NA, -length(train_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(train_z)=colnames(train_hdl_exp_top)
rownames(train_z)=rownames(train_hdl_exp_top)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
test_z=list()
test_z=foreach(i=1:dim(test_hdl_exp_top)[1]) %dopar%
  as.numeric(scale(as.numeric(test_hdl_exp_top[i,])))
test_z=as.data.frame(t(structure(test_z, row.names = c(NA, -length(test_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(test_z)=colnames(test_hdl_exp_top)
rownames(test_z)=rownames(test_hdl_exp_top)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
train_zvip=list()
train_zvip=foreach(i=1:dim(train_hdl_vip_top)[1]) %dopar%
  as.numeric(scale(as.numeric(train_hdl_vip_top[i,])))
train_zvip=as.data.frame(t(structure(train_zvip, row.names = c(NA, -length(train_zvip[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(train_zvip)=colnames(train_hdl_vip_top)
rownames(train_zvip)=rownames(train_hdl_vip_top)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
test_zvip=list()
test_zvip=foreach(i=1:dim(test_hdl_vip_top)[1]) %dopar%
  as.numeric(scale(as.numeric(test_hdl_vip_top[i,])))
test_zvip=as.data.frame(t(structure(test_zvip, row.names = c(NA, -length(test_zvip[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(test_zvip)=colnames(test_hdl_vip_top)
rownames(test_zvip)=rownames(test_hdl_vip_top)

# Let's run cross-validation random forest to identify the number of features required to plateau the error when inputing 500 candidate HDL MRs
# For expression
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
cv_hdl_exp_rndFor=list()
cv_hdl_exp_rndFor=foreach(i=1:12,.packages = 'randomForest') %dopar% {
  set.seed(i)
  cv_hdl_exp_rndFor[i]=rfcv(trainx=t(train_z),trainy=train_HDL,scale = F,step = -5,cv.fold = 5) # Test how many features are required before MSE plateaus
}
stopCluster(cl)

par(mfrow=c(3,4))
with(cv_hdl_exp_rndFor[[1]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~20 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_exp_rndFor[[2]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~145 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_exp_rndFor[[3]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~200 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_exp_rndFor[[4]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~100 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_exp_rndFor[[5]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~170 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_exp_rndFor[[6]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~25 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_exp_rndFor[[7]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~65 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_exp_rndFor[[8]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~85 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_exp_rndFor[[9]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~160 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_exp_rndFor[[10]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~230 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_exp_rndFor[[11]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~90 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_exp_rndFor[[12]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~95 features is the ideal balance between parsimony and accuracy/purity
# Let's average across all 12 different seeds for the 500 predictor rfcv() analysis with zvip
cv_hdl_exp_avg=rowMeans(cbind(cv_hdl_exp_rndFor[[1]]$error.cv,cv_hdl_exp_rndFor[[2]]$error.cv,cv_hdl_exp_rndFor[[3]]$error.cv,cv_hdl_exp_rndFor[[4]]$error.cv,
                             cv_hdl_exp_rndFor[[5]]$error.cv,cv_hdl_exp_rndFor[[6]]$error.cv,cv_hdl_exp_rndFor[[7]]$error.cv,cv_hdl_exp_rndFor[[8]]$error.cv,
                             cv_hdl_exp_rndFor[[9]]$error.cv,cv_hdl_exp_rndFor[[10]]$error.cv,cv_hdl_exp_rndFor[[11]]$error.cv,cv_hdl_exp_rndFor[[12]]$error.cv))
cv_hdl_exp_sd=rowSds(cbind(cv_hdl_exp_rndFor[[1]]$error.cv,cv_hdl_exp_rndFor[[2]]$error.cv,cv_hdl_exp_rndFor[[3]]$error.cv,cv_hdl_exp_rndFor[[4]]$error.cv,
                          cv_hdl_exp_rndFor[[5]]$error.cv,cv_hdl_exp_rndFor[[6]]$error.cv,cv_hdl_exp_rndFor[[7]]$error.cv,cv_hdl_exp_rndFor[[8]]$error.cv,
                          cv_hdl_exp_rndFor[[9]]$error.cv,cv_hdl_exp_rndFor[[10]]$error.cv,cv_hdl_exp_rndFor[[11]]$error.cv,cv_hdl_exp_rndFor[[12]]$error.cv))
pdf("./HDL/Eurobats_adipose_CVRF_analysis_of_HDL_regulators_by_expression.pdf",width = 9)
par(mfrow=c(1,1))
plot(x=as.numeric(names(cv_hdl_exp_avg)),y=cv_hdl_exp_avg,ylim=range(c(cv_hdl_exp_avg-cv_hdl_exp_sd,cv_hdl_exp_avg+cv_hdl_exp_sd)),pch=19,xlab="Feature count",ylab="Mean cross-validation error",main="Average expression rfcv() results (12 seeds)")
arrows(x0=as.numeric(names(cv_hdl_exp_avg)),y0=cv_hdl_exp_avg-cv_hdl_exp_sd,x1=as.numeric(names(cv_hdl_exp_avg)),y1=cv_hdl_exp_avg+cv_hdl_exp_sd,length = 0.05, angle = 90, code = 3) # This trick draws the error bars
dev.off()
# Again, 100 HDL MRs should work fine.

# For activity
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
cv_hdl_vip_rndFor=list()
cv_hdl_vip_rndFor=foreach(i=1:12,.packages = 'randomForest') %dopar% {
  set.seed(i)
  cv_hdl_vip_rndFor[i]=rfcv(trainx=t(train_zvip),trainy=train_HDL,scale = F,step = -5,cv.fold = 5) # Test how many features are required before MSE plateaus
}
stopCluster(cl)

par(mfrow=c(3,4))
with(cv_hdl_vip_rndFor[[1]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~20 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_vip_rndFor[[2]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~145 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_vip_rndFor[[3]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~200 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_vip_rndFor[[4]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~100 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_vip_rndFor[[5]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~170 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_vip_rndFor[[6]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~25 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_vip_rndFor[[7]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~65 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_vip_rndFor[[8]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~85 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_vip_rndFor[[9]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~160 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_vip_rndFor[[10]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~230 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_vip_rndFor[[11]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~90 features is the ideal balance between parsimony and accuracy/purity
with(cv_hdl_vip_rndFor[[12]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~95 features is the ideal balance between parsimony and accuracy/purity
# Let's average across all 12 different seeds for the 500 predictor rfcv() analysis with zvip
cv_hdl_vip_avg=rowMeans(cbind(cv_hdl_vip_rndFor[[1]]$error.cv,cv_hdl_vip_rndFor[[2]]$error.cv,cv_hdl_vip_rndFor[[3]]$error.cv,cv_hdl_vip_rndFor[[4]]$error.cv,
                              cv_hdl_vip_rndFor[[5]]$error.cv,cv_hdl_vip_rndFor[[6]]$error.cv,cv_hdl_vip_rndFor[[7]]$error.cv,cv_hdl_vip_rndFor[[8]]$error.cv,
                              cv_hdl_vip_rndFor[[9]]$error.cv,cv_hdl_vip_rndFor[[10]]$error.cv,cv_hdl_vip_rndFor[[11]]$error.cv,cv_hdl_vip_rndFor[[12]]$error.cv))
cv_hdl_vip_sd=rowSds(cbind(cv_hdl_vip_rndFor[[1]]$error.cv,cv_hdl_vip_rndFor[[2]]$error.cv,cv_hdl_vip_rndFor[[3]]$error.cv,cv_hdl_vip_rndFor[[4]]$error.cv,
                           cv_hdl_vip_rndFor[[5]]$error.cv,cv_hdl_vip_rndFor[[6]]$error.cv,cv_hdl_vip_rndFor[[7]]$error.cv,cv_hdl_vip_rndFor[[8]]$error.cv,
                           cv_hdl_vip_rndFor[[9]]$error.cv,cv_hdl_vip_rndFor[[10]]$error.cv,cv_hdl_vip_rndFor[[11]]$error.cv,cv_hdl_vip_rndFor[[12]]$error.cv))
pdf("./HDL/Eurobats_adipose_CVRF_analysis_of_HDL_regulators_by_activity.pdf",width = 9)
par(mfrow=c(1,1))
plot(x=as.numeric(names(cv_hdl_vip_avg)),y=cv_hdl_vip_avg,ylim=range(c(cv_hdl_vip_avg-cv_hdl_vip_sd,cv_hdl_vip_avg+cv_hdl_vip_sd)),pch=19,xlab="Feature count",ylab="Mean cross-validation error",main="Average activity rfcv() results (12 seeds)")
arrows(x0=as.numeric(names(cv_hdl_vip_avg)),y0=cv_hdl_vip_avg-cv_hdl_vip_sd,x1=as.numeric(names(cv_hdl_vip_avg)),y1=cv_hdl_vip_avg+cv_hdl_vip_sd,length = 0.05, angle = 90, code = 3) # This trick draws the error bars
dev.off()
# Again, 100 HDL MRs should work fine.

# Let's train a random forest with the 500 regulators to pick the top 100 (by importance type 1) as HDL MRs
# For expression
set.seed(246)
hdl_exp_rndFor500=randomForest(x=t(train_hdl_exp_top),y=train_HDL,ntree = 1000,keep.forest = T,importance = T)
par(mfrow=c(1,1))
plot(hdl_exp_rndFor500)
import_hdl500_exp=importance(hdl_exp_rndFor500,type = 1,scale = F)
import_hdl500_exp=import_hdl500_exp[order(import_hdl500_exp,decreasing = T),,drop=F]
train_hdl_exp_mrs=train_hdl_exp_top[rownames(import_hdl500_exp)[1:100],]
test_hdl_exp_mrs=test_hdl_exp_top[rownames(import_hdl500_exp)[1:100],]

# For activity
set.seed(246)
hdl_vip_rndFor500=randomForest(x=t(train_hdl_vip_top),y=train_HDL,ntree = 1000,keep.forest = T,importance = T)
par(mfrow=c(1,1))
plot(hdl_vip_rndFor500)
import_hdl500_vip=importance(hdl_vip_rndFor500,type = 1,scale = F)
import_hdl500_vip=import_hdl500_vip[order(import_hdl500_vip,decreasing = T),,drop=F]
train_hdl_vip_mrs=train_hdl_vip_top[rownames(import_hdl500_vip)[1:100],]
test_hdl_vip_mrs=test_hdl_vip_top[rownames(import_hdl500_vip)[1:100],]

# How many HDL MRs are in common between the expression and activity based RFs?
sum(rownames(train_hdl_exp_mrs) %in% rownames(train_hdl_vip_mrs)) # 19 putative HDL MRs are identical regardless of using expression or activities

# Let's train and test a final HDL random forest
# For expression
set.seed(314)
final_hdl_exp_rndFor=randomForest(x=t(train_hdl_exp_mrs),y=train_HDL,xtest = t(test_hdl_exp_mrs),ytest = test_HDL,ntree = 1000,keep.forest = T,importance = T)
train_lm_hdl_exp=lm(train_HDL~final_hdl_exp_rndFor$predicted) # Beta=1.06247, P=3.58e-42, r-squared=0.3157
test_lm_hdl_exp=lm(test_HDL~final_hdl_exp_rndFor$test$predicted) # Beta=1.03438, P=2.38e-17, r-squared=0.2891

# For activity
set.seed(314)
final_hdl_vip_rndFor=randomForest(x=t(train_hdl_vip_mrs),y=train_HDL,xtest = t(test_hdl_vip_mrs),ytest = test_HDL,ntree = 1000,keep.forest = T,importance = T)
train_lm_hdl_vip=lm(train_HDL~final_hdl_vip_rndFor$predicted) # Beta=0.9546, P=5.96e-36, r-squared=0.2744
test_lm_hdl_vip=lm(test_HDL~final_hdl_vip_rndFor$test$predicted) # Beta=0.9293, P=8.91e-16, r-squared=0.2642
# Very similar results to the old MR analysis with missing regulators, suggesting the analysis is pretty robust.

# Let's make some hi-res plots
pdf("./HDL/Eurobats_adipose_HDL_MR_by_expression_importance_plot.pdf",width = 6,height = 10)
varImpPlot(final_hdl_exp_rndFor,type=1,scale = F,n.var = 100,cex = 0.5)
dev.off()

pdf("./HDL/Eurobats_adipose_HDL_MR_by_activity_importance_plot.pdf",width = 6,height = 10)
varImpPlot(final_hdl_vip_rndFor,type=1,scale = F,n.var = 100,cex = 0.5)
dev.off()

### Triglycerides
# Filter to overlapping samples for Triglycerides and VIPER data
filt_phenos=euro_pheno[na.omit(match(colnames(euro_vip),rownames(euro_pheno))),]
filt_phenos=filt_phenos[!is.na(filt_phenos$TotalTriglycerides),] # 698 samples
filt_vip=euro_vip[,rownames(filt_phenos)]
filt_exp=euro_exp[,colnames(filt_vip)]
all(colnames(filt_vip)==rownames(filt_phenos)) # TRUE
all(colnames(filt_vip)==colnames(filt_exp)) # TRUE

# Split data 70/30 for training/test sets. 
RNGkind(sample.kind = "default") # Note that the original TriG RF MR analysis was performed in R >3.6.0, so it used the new default sampling.
set.seed(123)
rnd_samples=sample(colnames(filt_vip),size = 488) # 70% of 697 samples
train_exp=filt_exp[,rnd_samples]
train_vip=filt_vip[,rnd_samples]
train_TriG=filt_phenos[rnd_samples,"TotalTriglycerides"]
test_exp=filt_exp[,!(colnames(filt_exp) %in% rnd_samples)]
test_vip=filt_vip[,!(colnames(filt_vip) %in% rnd_samples)]
test_TriG=filt_phenos[!(rownames(filt_phenos) %in% rnd_samples),"TotalTriglycerides"]

# Identify the regulators whose expression or activities are best associated with TriG
# For expression
e_TriG=list()
for(i in 1:dim(train_exp)[1]){
  e_TriG[[i]]=summary(lm(train_TriG~as.numeric(train_exp[i,])))
}
names(e_TriG)=rownames(train_exp)
sum_e_TriG=as.data.frame(matrix(nrow = length(e_TriG),ncol = 2))
for(i in 1:length(e_TriG)){
  sum_e_TriG[i,1]=coef(e_TriG[[i]])[2,1]
  sum_e_TriG[i,2]=coef(e_TriG[[i]])[2,4]
}
colnames(sum_e_TriG)=c("Beta","P")
rownames(sum_e_TriG)=names(e_TriG)
sum_e_TriG=sum_e_TriG[order(sum_e_TriG$P),]
sum_e_TriG$BonferroniP=sum_e_TriG$P*dim(sum_e_TriG)[1]
sum(sum_e_TriG$BonferroniP<0.05) # 1113 regulators w/ Bonferroni adjusted P < 0.05
sum(sum_e_TriG$P<0.05) # 2672 regulators w/ P < 0.05

# For activity
a_TriG=list()
for(i in 1:dim(train_vip)[1]){
  a_TriG[[i]]=summary(lm(train_TriG~as.numeric(train_vip[i,])))
}
names(a_TriG)=rownames(train_vip)
sum_a_TriG=as.data.frame(matrix(nrow = length(a_TriG),ncol = 2))
for(i in 1:length(a_TriG)){
  sum_a_TriG[i,1]=coef(a_TriG[[i]])[2,1]
  sum_a_TriG[i,2]=coef(a_TriG[[i]])[2,4]
}
colnames(sum_a_TriG)=c("Beta","P")
rownames(sum_a_TriG)=names(a_TriG)
sum_a_TriG=sum_a_TriG[order(sum_a_TriG$P),]
sum_a_TriG$BonferroniP=sum_a_TriG$P*dim(sum_a_TriG)[1]
sum(sum_a_TriG$BonferroniP<0.05) # 1543 regulators w/ Bonferroni adjusted P < 0.05
sum(sum_a_TriG$P<0.05) # 2910 regulators w/ P < 0.05

# Let's use the top 500 to train our random forest
train_triG_exp_top=train_exp[rownames(sum_e_TriG)[1:500],]
train_triG_vip_top=train_vip[rownames(sum_a_TriG)[1:500],]
test_triG_exp_top=test_exp[rownames(sum_e_TriG)[1:500],]
test_triG_vip_top=test_vip[rownames(sum_a_TriG)[1:500],]

# Let's Z-transform the activities prior to RF modeling
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
train_z=list()
train_z=foreach(i=1:dim(train_triG_exp_top)[1]) %dopar%
  as.numeric(scale(as.numeric(train_triG_exp_top[i,])))
train_z=as.data.frame(t(structure(train_z, row.names = c(NA, -length(train_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(train_z)=colnames(train_triG_exp_top)
rownames(train_z)=rownames(train_triG_exp_top)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
test_z=list()
test_z=foreach(i=1:dim(test_triG_exp_top)[1]) %dopar%
  as.numeric(scale(as.numeric(test_triG_exp_top[i,])))
test_z=as.data.frame(t(structure(test_z, row.names = c(NA, -length(test_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(test_z)=colnames(test_triG_exp_top)
rownames(test_z)=rownames(test_triG_exp_top)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
train_zvip=list()
train_zvip=foreach(i=1:dim(train_triG_vip_top)[1]) %dopar%
  as.numeric(scale(as.numeric(train_triG_vip_top[i,])))
train_zvip=as.data.frame(t(structure(train_zvip, row.names = c(NA, -length(train_zvip[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(train_zvip)=colnames(train_triG_vip_top)
rownames(train_zvip)=rownames(train_triG_vip_top)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
test_zvip=list()
test_zvip=foreach(i=1:dim(test_triG_vip_top)[1]) %dopar%
  as.numeric(scale(as.numeric(test_triG_vip_top[i,])))
test_zvip=as.data.frame(t(structure(test_zvip, row.names = c(NA, -length(test_zvip[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(test_zvip)=colnames(test_triG_vip_top)
rownames(test_zvip)=rownames(test_triG_vip_top)

# Let's run cross-validation random forest to identify the number of features required to plateau the error when inputing 500 candidate Triglycerides MRs
# For expression
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
cv_triG_exp_rndFor=list()
cv_triG_exp_rndFor=foreach(i=1:12,.packages = 'randomForest') %dopar% {
  set.seed(i)
  cv_triG_exp_rndFor[i]=rfcv(trainx=t(train_z),trainy=train_TriG,scale = F,step = -5,cv.fold = 5) # Test how many features are required before MSE plateaus
}
stopCluster(cl)

par(mfrow=c(3,4))
with(cv_triG_exp_rndFor[[1]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~20 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_exp_rndFor[[2]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~145 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_exp_rndFor[[3]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~200 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_exp_rndFor[[4]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~100 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_exp_rndFor[[5]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~170 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_exp_rndFor[[6]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~25 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_exp_rndFor[[7]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~65 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_exp_rndFor[[8]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~85 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_exp_rndFor[[9]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~160 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_exp_rndFor[[10]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~230 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_exp_rndFor[[11]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~90 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_exp_rndFor[[12]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~95 features is the ideal balance between parsimony and accuracy/purity
# Let's average across all 12 different seeds for the 500 predictor rfcv() analysis with zvip
cv_triG_exp_avg=rowMeans(cbind(cv_triG_exp_rndFor[[1]]$error.cv,cv_triG_exp_rndFor[[2]]$error.cv,cv_triG_exp_rndFor[[3]]$error.cv,cv_triG_exp_rndFor[[4]]$error.cv,
                           cv_triG_exp_rndFor[[5]]$error.cv,cv_triG_exp_rndFor[[6]]$error.cv,cv_triG_exp_rndFor[[7]]$error.cv,cv_triG_exp_rndFor[[8]]$error.cv,
                           cv_triG_exp_rndFor[[9]]$error.cv,cv_triG_exp_rndFor[[10]]$error.cv,cv_triG_exp_rndFor[[11]]$error.cv,cv_triG_exp_rndFor[[12]]$error.cv))
cv_triG_exp_sd=rowSds(cbind(cv_triG_exp_rndFor[[1]]$error.cv,cv_triG_exp_rndFor[[2]]$error.cv,cv_triG_exp_rndFor[[3]]$error.cv,cv_triG_exp_rndFor[[4]]$error.cv,
                        cv_triG_exp_rndFor[[5]]$error.cv,cv_triG_exp_rndFor[[6]]$error.cv,cv_triG_exp_rndFor[[7]]$error.cv,cv_triG_exp_rndFor[[8]]$error.cv,
                        cv_triG_exp_rndFor[[9]]$error.cv,cv_triG_exp_rndFor[[10]]$error.cv,cv_triG_exp_rndFor[[11]]$error.cv,cv_triG_exp_rndFor[[12]]$error.cv))
pdf("./Triglycerides/Eurobats_adipose_CVRF_analysis_of_TriG_regulators_by_expression.pdf",width = 9)
par(mfrow=c(1,1))
plot(x=as.numeric(names(cv_triG_exp_avg)),y=cv_triG_exp_avg,ylim=range(c(cv_triG_exp_avg-cv_triG_exp_sd,cv_triG_exp_avg+cv_triG_exp_sd)),pch=19,xlab="Feature count",ylab="Mean cross-validation error",main="Average expression rfcv() results (12 seeds)")
arrows(x0=as.numeric(names(cv_triG_exp_avg)),y0=cv_triG_exp_avg-cv_triG_exp_sd,x1=as.numeric(names(cv_triG_exp_avg)),y1=cv_triG_exp_avg+cv_triG_exp_sd,length = 0.05, angle = 90, code = 3) # This trick draws the error bars
dev.off()
# In this case, ~60 regulators was the ideal feature count

# For activity
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
cv_triG_vip_rndFor=list()
cv_triG_vip_rndFor=foreach(i=1:12,.packages = 'randomForest') %dopar% {
  set.seed(i)
  cv_triG_vip_rndFor[i]=rfcv(trainx=t(train_zvip),trainy=train_TriG,scale = F,step = -5,cv.fold = 5) # Test how many features are required before MSE plateaus
}
stopCluster(cl)

par(mfrow=c(3,4))
with(cv_triG_vip_rndFor[[1]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~20 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_vip_rndFor[[2]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~145 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_vip_rndFor[[3]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~200 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_vip_rndFor[[4]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~100 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_vip_rndFor[[5]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~170 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_vip_rndFor[[6]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~25 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_vip_rndFor[[7]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~65 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_vip_rndFor[[8]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~85 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_vip_rndFor[[9]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~160 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_vip_rndFor[[10]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~230 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_vip_rndFor[[11]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~90 features is the ideal balance between parsimony and accuracy/purity
with(cv_triG_vip_rndFor[[12]],plot(n.var,error.cv,log="x",type="o",lwd=2)) # This suggests ~95 features is the ideal balance between parsimony and accuracy/purity
# Let's average across all 12 different seeds for the 500 predictor rfcv() analysis with zvip
cv_triG_vip_avg=rowMeans(cbind(cv_triG_vip_rndFor[[1]]$error.cv,cv_triG_vip_rndFor[[2]]$error.cv,cv_triG_vip_rndFor[[3]]$error.cv,cv_triG_vip_rndFor[[4]]$error.cv,
                          cv_triG_vip_rndFor[[5]]$error.cv,cv_triG_vip_rndFor[[6]]$error.cv,cv_triG_vip_rndFor[[7]]$error.cv,cv_triG_vip_rndFor[[8]]$error.cv,
                          cv_triG_vip_rndFor[[9]]$error.cv,cv_triG_vip_rndFor[[10]]$error.cv,cv_triG_vip_rndFor[[11]]$error.cv,cv_triG_vip_rndFor[[12]]$error.cv))
cv_triG_vip_sd=rowSds(cbind(cv_triG_vip_rndFor[[1]]$error.cv,cv_triG_vip_rndFor[[2]]$error.cv,cv_triG_vip_rndFor[[3]]$error.cv,cv_triG_vip_rndFor[[4]]$error.cv,
                       cv_triG_vip_rndFor[[5]]$error.cv,cv_triG_vip_rndFor[[6]]$error.cv,cv_triG_vip_rndFor[[7]]$error.cv,cv_triG_vip_rndFor[[8]]$error.cv,
                       cv_triG_vip_rndFor[[9]]$error.cv,cv_triG_vip_rndFor[[10]]$error.cv,cv_triG_vip_rndFor[[11]]$error.cv,cv_triG_vip_rndFor[[12]]$error.cv))
pdf("./Triglycerides/Eurobats_adipose_CVRF_analysis_of_TriG_regulators_by_activity.pdf",width = 9)
par(mfrow=c(1,1))
plot(x=as.numeric(names(cv_triG_vip_avg)),y=cv_triG_vip_avg,ylim=range(c(cv_triG_vip_avg-cv_triG_vip_sd,cv_triG_vip_avg+cv_triG_vip_sd)),pch=19,xlab="Feature count",ylab="Mean cross-validation error",main="Average activity rfcv() results (12 seeds)")
arrows(x0=as.numeric(names(cv_triG_vip_avg)),y0=cv_triG_vip_avg-cv_triG_vip_sd,x1=as.numeric(names(cv_triG_vip_avg)),y1=cv_triG_vip_avg+cv_triG_vip_sd,length = 0.05, angle = 90, code = 3) # This trick draws the error bars
dev.off()
# In this case, ~60 regulators was the ideal feature count.

# Let's train a random forest with the 500 regulators to pick the top 60 (by importance type 1) as TriG MRs
# For expression
set.seed(246)
triG_exp_rndFor500=randomForest(x=t(train_triG_exp_top),y=train_TriG,ntree = 1000,keep.forest = T,importance = T)
par(mfrow=c(1,1))
plot(triG_exp_rndFor500)
import_triG500_exp=importance(triG_exp_rndFor500,type = 1,scale = F)
import_triG500_exp=import_triG500_exp[order(import_triG500_exp,decreasing = T),,drop=F]
train_triG_exp_mrs=train_triG_exp_top[rownames(import_triG500_exp)[1:60],]
test_triG_exp_mrs=test_triG_exp_top[rownames(import_triG500_exp)[1:60],]

# For activity
set.seed(246)
triG_vip_rndFor500=randomForest(x=t(train_triG_vip_top),y=train_TriG,ntree = 1000,keep.forest = T,importance = T)
par(mfrow=c(1,1))
plot(triG_vip_rndFor500)
import_triG500_vip=importance(triG_vip_rndFor500,type = 1,scale = F)
import_triG500_vip=import_triG500_vip[order(import_triG500_vip,decreasing = T),,drop=F]
train_triG_vip_mrs=train_triG_vip_top[rownames(import_triG500_vip)[1:60],]
test_triG_vip_mrs=test_triG_vip_top[rownames(import_triG500_vip)[1:60],]

# How many TriG MRs are in common between the expression and activity based RFs?
sum(rownames(train_triG_exp_mrs) %in% rownames(train_triG_vip_mrs)) # 11 putative TriG MRs are identical regardless of using expression or activities

# Let's train and test a final TriG random forest
# For expression
set.seed(314)
final_triG_exp_rndFor=randomForest(x=t(train_triG_exp_mrs),y=train_TriG,xtest = t(test_triG_exp_mrs),ytest = test_TriG,ntree = 1000,keep.forest = T,importance = T)
train_lm_triG_exp=lm(train_TriG~final_triG_exp_rndFor$predicted) # Beta=1.12947, P=1.72e-62, r-squared=0.4352
test_lm_triG_exp=lm(test_TriG~final_triG_exp_rndFor$test$predicted) # Beta=1.03320, P=7.58-20, r-squared=0.3269

# For activity
set.seed(314)
final_triG_vip_rndFor=randomForest(x=t(train_triG_vip_mrs),y=train_TriG,xtest = t(test_triG_vip_mrs),ytest = test_TriG,ntree = 1000,keep.forest = T,importance = T)
train_lm_triG_vip=lm(train_TriG~final_triG_vip_rndFor$predicted) # Beta=0.9115, P=1.18e-43, r-squared=0.3235
test_lm_triG_vip=lm(test_TriG~final_triG_vip_rndFor$test$predicted) # Beta=0.8156, P=1.36-16, r-squared=0.2772

# Let's make some hi-res plots
pdf("./Triglycerides/Eurobats_adipose_TriG_MR_by_expression_importance_plot.pdf",width = 6,height = 6)
varImpPlot(final_triG_exp_rndFor,type=1,scale = F,n.var = 60,cex = 0.5)
dev.off()

pdf("./Triglycerides/Eurobats_adipose_TriG_MR_by_activity_importance_plot.pdf",width = 6,height = 6)
varImpPlot(final_triG_vip_rndFor,type=1,scale = F,n.var = 60,cex = 0.5)
dev.off()

# Write all of the MR lists to file
write.table(rownames(train_bmi_vip_mrs),"./BMI/Eurobats_adipose_time-matched_BMI_MRs_from_RF_modeling_by_activity.txt",col.names = F,row.names = F,quote = F)
write.table(rownames(train_whr_vip_mrs),"./WHR/Eurobats_adipose_pseudo-time-matched_WHR_MRs_from_RF_modeling_by_activity_in_subjects_with_steady_BMI.txt",col.names = F,row.names = F,quote = F)
write.table(rownames(train_homair_vip_mrs),"./HOMA-IR/Eurobats_adipose_time-matched_HOMA-IR_MRs_from_RF_modeling_by_activity.txt",col.names = F,row.names = F,quote = F)
write.table(rownames(train_hdl_vip_mrs),"./HDL/Eurobats_adipose_time-matched_HDL_MRs_from_RF_modeling_by_activity.txt",col.names = F,row.names = F,quote = F)
write.table(rownames(train_triG_vip_mrs),"./Triglycerides/Eurobats_adipose_time-matched_Triglycerides_MRs_from_RF_modeling_by_activity.txt",col.names = F,row.names = F,quote = F)
write.table(rownames(train_bmi_exp_mrs),"./BMI/Eurobats_adipose_time-matched_BMI_MRs_from_RF_modeling_by_expression.txt",col.names = F,row.names = F,quote = F)
write.table(rownames(train_whr_exp_mrs),"./WHR/Eurobats_adipose_pseudo-time-matched_WHR_MRs_from_RF_modeling_by_expression_in_subjects_with_steady_BMI.txt",col.names = F,row.names = F,quote = F)
write.table(rownames(train_homair_exp_mrs),"./HOMA-IR/Eurobats_adipose_time-matched_HOMA-IR_MRs_from_RF_modeling_by_expression.txt",col.names = F,row.names = F,quote = F)
write.table(rownames(train_hdl_exp_mrs),"./HDL/Eurobats_adipose_time-matched_HDL_MRs_from_RF_modeling_by_expression.txt",col.names = F,row.names = F,quote = F)
write.table(rownames(train_triG_exp_mrs),"./Triglycerides/Eurobats_adipose_time-matched_Triglycerides_MRs_from_RF_modeling_by_expression.txt",col.names = F,row.names = F,quote = F)

# Finally, let's make an MR comparison table for the Venn diagram
# For MRs by activity
all_mrs=rownames(train_bmi_vip_mrs)
all_mrs=append(all_mrs, rownames(train_whr_vip_mrs))
all_mrs=append(all_mrs, rownames(train_homair_vip_mrs))
all_mrs=append(all_mrs, rownames(train_hdl_vip_mrs))
all_mrs=append(all_mrs, rownames(train_triG_vip_mrs))
all_mrs=unique(all_mrs) # 291 unique MRs
mr_comp=as.data.frame(matrix(nrow=291,ncol = 5),row.names = all_mrs)
colnames(mr_comp)=c("BMI_MR","WHR_MR","HOMA-IR_MR","HDL_MR","TriG_MR")
mr_comp$BMI_MR=(all_mrs %in% rownames(train_bmi_vip_mrs))
mr_comp$WHR_MR=(all_mrs %in% rownames(train_whr_vip_mrs))
mr_comp$`HOMA-IR_MR`=(all_mrs %in% rownames(train_homair_vip_mrs))
mr_comp$HDL_MR=(all_mrs %in% rownames(train_hdl_vip_mrs))
mr_comp$TriG_MR=(all_mrs %in% rownames(train_triG_vip_mrs))

# For MRs by expression
all_mrs_exp=rownames(train_bmi_exp_mrs)
all_mrs_exp=append(all_mrs_exp, rownames(train_whr_exp_mrs))
all_mrs_exp=append(all_mrs_exp, rownames(train_homair_exp_mrs))
all_mrs_exp=append(all_mrs_exp, rownames(train_hdl_exp_mrs))
all_mrs_exp=append(all_mrs_exp, rownames(train_triG_exp_mrs))
all_mrs_exp=unique(all_mrs_exp) # 305 unique MRs
mr_exp_comp=as.data.frame(matrix(nrow=305,ncol = 5),row.names = all_mrs_exp)
colnames(mr_exp_comp)=c("BMI_MR","WHR_MR","HOMA-IR_MR","HDL_MR","TriG_MR")
mr_exp_comp$BMI_MR=(all_mrs_exp %in% rownames(train_bmi_exp_mrs))
mr_exp_comp$WHR_MR=(all_mrs_exp %in% rownames(train_whr_exp_mrs))
mr_exp_comp$`HOMA-IR_MR`=(all_mrs_exp %in% rownames(train_homair_exp_mrs))
mr_exp_comp$HDL_MR=(all_mrs_exp %in% rownames(train_hdl_exp_mrs))
mr_exp_comp$TriG_MR=(all_mrs_exp %in% rownames(train_triG_exp_mrs))

write.table(cbind("MR"=rownames(mr_comp),mr_comp),"MR_comparison_table.txt",sep = "\t",row.names = F,col.names = T,quote = F)
write.table(cbind("MR"=rownames(mr_exp_comp),mr_exp_comp),"MRs_BY_EXPRESSION_comparison_table.txt",sep = "\t",row.names = F,col.names = T,quote = F)

### Now it's time to validate the new MR analyses with the METSIM data

# Read in the unprocessed data from GEO for extracting phenotype data
pheno=read.table("../METSIM expression array data/GSE70353_series_matrix.txt",sep = "\t",header = F,skip = 24,fill = T)

# Format phenotype data
rownames(pheno[pheno[,1]=="ID_REF",]) # 54
pheno=pheno[1:53,-1]
pheno=pheno[-c(36:53),]
colnames(pheno)=t(pheno[2,])[,1]
pheno=pheno[-c(1:10),]
pheno=droplevels.data.frame(pheno)
rownames(pheno)=gsub(":.*","",pheno[,1])
for(i in colnames(pheno)){
  pheno[,i]=gsub(".*: ","",pheno[,i])
}
pheno=as.data.frame(t(pheno))
pheno=pheno[,-1]
for(i in 1:dim(pheno)[2]){
  pheno[,i]=as.numeric(pheno[,i])
}
# Some NAs were introduced by the coercion. Let's check where.
for(i in 1:dim(pheno)[2]){
  print(paste(colnames(pheno)[i],sum(is.na(pheno[,i])),sep=":"))
}
# There are NAs for "totfa", "matsuda", "p_crp", "s_ldlc","homair" and "P_ins0", but the only NA that would be a problem is 1 for HOMA-IR, 
# so let's drop that sample
pheno=pheno[!is.na(pheno$homair),]

# Read in the METSIM expression data
metsim=read.table("../METSIM expression array data/METSIM_probe_average_gene_expression_for_ARACNe.txt",sep = "\t",header = T, row.names = 1)

# Filter metsim to same samples as in pheno
metsim=metsim[,rownames(pheno)]
all(rownames(pheno)==colnames(metsim)) # TRUE

# Let's check that the METSIM expression data is well normalized and without batch effects

# sample-wise distributions
par(mfrow=c(5,3),mar=c(1,1,1,1))
boxplot(metsim[,1:50])
boxplot(metsim[,51:100])
boxplot(metsim[,101:150])
boxplot(metsim[,151:200])
boxplot(metsim[,201:250])
boxplot(metsim[,251:300])
boxplot(metsim[,301:350])
boxplot(metsim[,351:400])
boxplot(metsim[,401:450])
boxplot(metsim[,451:500])
boxplot(metsim[,501:550])
boxplot(metsim[,551:600])
boxplot(metsim[,601:650])
boxplot(metsim[,651:699])
boxplot(metsim[,700:769])
# All samples have very comparable distributions

# Now let's check the heteroscedacity
meanSdPlot(as.matrix(metsim),bins = 500)$gg + scale_y_continuous(limits = c(0,2.5))
# No problem

# PCA to check for batch effects
plotPCA4(as.matrix(metsim),colorVar = pheno$bmi)
plotPCA4(as.matrix(metsim),colorVar = pheno$whr)
plotPCA4(as.matrix(metsim),colorVar = pheno$homair)
plotPCA4(as.matrix(metsim),colorVar = pheno$s_hdlc)
plotPCA4(as.matrix(metsim),colorVar = pheno$s_tottg)
# The samples are distributed in a triangular shape in the PC1 vs. PC2 plane, but they seem to be a 
# single, continuous cluster. There does seem to be some correlation between each of the tested phenotypes
# and PC2 (and maybe to a lesser extent with PC1).

# Regenerate the Eurobats adipose interactome
adipo=read.table("../Adipose expression data/FINAL_logTPMs_and_activities/Eurobats_adipose_expressed_genes_logTPM.txt",header = T, row.names = 1,sep = "\t")
adipo_set=ExpressionSet(assayData=as.matrix(adipo))
adipo_regulon=aracne2regulon("../Adipose expression data/FINAL_logTPMs_and_activities/Eurobats_adipose_900bootstraps_ARACNe_network_with_LINC-PINT_no_header.txt",adipo_set,format="3col") # NOTE: The aracne network file cannot have a header!!!
aracne=read.table("../Adipose expression data/FINAL_logTPMs_and_activities/Eurobats_adipose_900bootstraps_ARACNe_network_with_LINC-PINT_no_header.txt",sep = "\t",header = F)

# Let's check for missing regulators in the METSIM data
miss_reg=unique(as.character(aracne[!(aracne[,1] %in% rownames(metsim)),1])) 
length(unique(aracne[,1])) # 4221 regulators in the Eurobats adipose interactome
# There are 123/4221 regulators missing from the METSIM expression data

# Scale (i.e. Z-transform) genes in the METSIM dataset.
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
metsim_z=list()
metsim_z=foreach(i=1:dim(metsim)[1]) %dopar%
  as.numeric(scale(as.numeric(metsim[i,])))
metsim_z=as.data.frame(t(structure(metsim_z, row.names = c(NA, -length(metsim_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(metsim_z)=colnames(metsim)
rownames(metsim_z)=rownames(metsim)
met_set=ExpressionSet(assayData=as.matrix(metsim_z))

# Run VIPER activity inference on METSIM data using the TwinsUK adipose interactome
met_viper=viper(met_set,adipo_regulon)
met_vip=as.data.frame(exprs(met_viper))

# Let's check for missing regulators in the METSIM data
miss_reg_vip=unique(as.character(aracne[!(aracne[,1] %in% rownames(met_vip)),1]))
miss_reg_targs=unique(as.character(aracne[na.omit(match(miss_reg_vip,aracne[,1])),2]))

# Z-transform the activities
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
met_zvip=list()
met_zvip=foreach(i=1:dim(met_vip)[1]) %dopar%
  as.numeric(scale(as.numeric(met_vip[i,])))
met_zvip=as.data.frame(t(structure(met_zvip, row.names = c(NA, -length(met_zvip[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(met_zvip)=colnames(met_vip)
rownames(met_zvip)=rownames(met_vip)

# Now we can finally predict the METSIM BMIs
# By expression
# This failed due to missing putative BMI MRs in the METSIM expression data
unique(as.character(rownames(final_bmi_exp_rndFor$importance)[!(rownames(final_bmi_exp_rndFor$importance) %in% rownames(metsim_z))]))
# TMEM189-UBE2V1 is missing from the METSIM expression data. To get around this, I will add the TMEM-UBE2V1
# Z-scaled activities into metsim_z as a proxy (which is reasonable assuming its activities and expression levels are highly correlated)
metsim_z_plus=rbind(metsim_z,met_zvip["TMEM189-UBE2V1",])
# Now we can predict with the BMI MR by expression RF model
met_pred_bmi_exp=predict(final_bmi_exp_rndFor,t(metsim_z_plus))
met_lm_bmi_exp=lm(pheno$bmi~met_pred_bmi_exp) # beta=0.8554, P=2.64e-69, r2=0.3314
# This is very comparable to the prediction accuracy based on activities.

# By activities
met_pred_bmi_vip=predict(final_bmi_vip_rndFor,t(met_zvip))
met_lm_bmi_vip=lm(pheno$bmi~met_pred_bmi_vip) # beta=0.6075, P=1.71E-71, r2=0.3402 


# Now we can finally predict the METSIM WHRs
# By expression
# This failed due to missing putative WHR MRs in the METSIM expression data
unique(as.character(rownames(final_whr_exp_rndFor$importance)[!(rownames(final_whr_exp_rndFor$importance) %in% rownames(metsim_z))]))
# HIST1H2BC is missing from the METSIM expression data. To get around this, I will add the HIST1H2BC
# Z-scaled activities into metsim_z_plus as a proxy (which is reasonable assuming its activities and expression levels are highly correlated)
metsim_z_plus=rbind(metsim_z_plus,met_zvip["HIST1H2BC",])
# Now we can predict with the WHR MR by expression RF model
met_pred_whr_exp=predict(final_whr_exp_rndFor,t(metsim_z_plus))
met_lm_whr_exp=lm(pheno$whr~met_pred_whr_exp) # beta=1.44846, P=5.32e-61, r2=0.2974 
# This is very comparable to the prediction accuracy based on activities.

# By activity
met_pred_whr_vip=predict(final_whr_vip_rndFor,t(met_vip))
met_lm_whr_vip=lm(pheno$whr~met_pred_whr_vip) # beta=1.4926, P=3.05E-57, r2=0.2814 

# Let's quickly compare Z-scaled predicted and actual WHRs, which should only effect the beta and the plot scales.
actual_zwhr=scale(pheno$whr)
predicted_zwhr=scale(met_pred_whr)
met_lm_zwhr=lm(actual_zwhr~predicted_zwhr) # beta=0.5313, P=3.05E-57, r2=0.2814 
# The beta is not really much better (i.e. closer to 1), but it is a little more intuitive to interpret. Basically, the WHR MR RF model systematically
# under-predicts the METSIM WHR, but the predictions are extremely well correlated with the actual WHR. This is probably largely due to the fact that
# the UK Twins study (Eurobats' source) is all female while METSIM is all male. Men and women differ in WHR distributions. There may also be a 
# calibration/measurement difference between TwinsUK and METSIM.

# Now we can finally predict the METSIM HOMA-IR
pheno$homair=log(pheno$homair) # Need to ln transform the METSIM HOMA-IR as I did for Eurobats
# By expression
# This failed due to missing putative HOMA-IR MRs in the METSIM expression data
unique(as.character(rownames(final_homair_exp_rndFor$importance)[!(rownames(final_homair_exp_rndFor$importance) %in% rownames(metsim_z))]))
# TAX1BP3 and GABRE are missing from the METSIM expression data. To get around this, I will add their
# Z-scaled activities into metsim_z_plus as a proxy (which is reasonable assuming its activities and expression levels are highly correlated)
metsim_z_plus=rbind(metsim_z_plus,met_zvip["TAX1BP3",])
metsim_z_plus=rbind(metsim_z_plus,met_zvip["GABRE",])
# Now we can predict with the HOMA-IR MR by expression RF model
met_pred_homair_exp=predict(final_homair_exp_rndFor,t(metsim_z_plus))
met_lm_homair_exp=lm(pheno$homair~met_pred_homair_exp) # beta=1.32031, P=4.08e-89, r2=0.4063 
# This actually a little better than with the HOMA-IR MR RF model based on activities.

# By activity
met_pred_homair_vip=predict(final_homair_vip_rndFor,t(met_vip))
met_lm_homair_vip=lm(pheno$homair~met_pred_homair_vip) # beta=1.5281, P=7.52E-79, r2=0.3686 

# Now we can finally predict the METSIM HDL
# By expression
# This failed due to missing putative HDL MRs in the METSIM expression data
unique(as.character(rownames(final_hdl_exp_rndFor$importance)[!(rownames(final_hdl_exp_rndFor$importance) %in% rownames(metsim_z))]))
# TAX1BP3, DDR1 and MAP2K3 are missing from the METSIM expression data. To get around this, I will add their
# Z-scaled activities into metsim_z_plus as a proxy (which is reasonable assuming its activities and expression levels are highly correlated)
metsim_z_plus=rbind(metsim_z_plus,met_zvip["DDR1",])
metsim_z_plus=rbind(metsim_z_plus,met_zvip["MAP2K3",])
# Now we can predict with the HDL MR by expression RF model
met_pred_hdl_exp=predict(final_hdl_exp_rndFor,t(metsim_z_plus))
met_lm_hdl_exp=lm(pheno$s_hdlc~met_pred_hdl_exp) # beta=0.800213, P=5.78e-20, r2=0.1023 
# This is very comparable to the prediction accuracy based on activities.

# By activity
met_pred_hdl_vip=predict(final_hdl_vip_rndFor,t(met_vip))
met_lm_hdl_vip=lm(pheno$s_hdlc~met_pred_hdl_vip) # beta=0.9372, P=9.81E-19, r2=0.0957 
plot(pheno$s_hdlc~met_pred_hdl_vip,xlab="Predicted HDL",ylab="Actual HDL",main="METSIM HDL prediction")
abline(met_lm_hdl_vip,col="red")
# Again the scales appears quite different between the predicted and actual HDL, which makes me wonder if these measurements are the same units or by the same
# test. In this case the beta is very close to 1, the P is quite significant, but the r2 leaves something to be desired. The P and r2 is a bit worse here than
# in the old analysis with missing regulators. The poor r2 seems to be largely due to the extreme actual HDLs (i.e. those above ~2.3). I compared the ranges 
# of HDL measurements for Eurobats and METSIM (see below), and they seem very comparable. Therefore, it is the RF model that seems to be rather conservative, 
# squeezing its predictions of HDLs to a tighter range than we'd expect.
range(filt_phenos$HDLcholesterol) # 0.78 to 3.5
range(pheno$s_hdlc) # 0.64 to 3.75
range(met_pred_hdl) # 1.42 to 2.26
# This may again be partially due to the gender difference between TwinsUK and METSIM.

# Now we can finally predict the METSIM Triglycerides
# By expression
met_pred_triG_exp=predict(final_triG_exp_rndFor,t(metsim_z))
met_lm_triG_exp=lm(pheno$s_tottg~met_pred_triG_exp) # beta=1.2581, P=2.40e-23, r2=0.1201 
# This is very comparable to the prediction accuracy based on activities.

# By activity
met_pred_trig_vip=predict(final_triG_vip_rndFor,t(met_vip))
met_lm_trig_vip=lm(pheno$s_tottg~met_pred_trig_vip) # beta=1.6804, P=2.45E-18, r2=0.0936 
# Again the scales appear quite different between the predicted and actual TriG, which makes me wonder if these measurements are the same units or by the same
# test. In this case the beta is very close to 1, the P is quite significant, but the r2 leaves something to be desired. The P and r2 is a bit worse here than
# in the old analysis with missing regulators. This seems to be entirely due to the extreme actual TriGs (i.e. those above ~2). I compared the ranges of TriG 
# measurements for Eurobats and METSIM (see below), and METSIM has much more extreme measurements. Regardless, the RF model seems to again be rather 
# conservative, squeezing its predictions of TriGs to a tighter range than we'd expect.
range(filt_phenos$TotalTriglycerides) # 0.30 to 4.70
range(pheno$s_tottg) # 0.39 to 9.74
range(met_pred_trig) # 0.67 to 1.61
# This may again be partially due to the gender difference between TwinsUK and METSIM.

# Now let's make some hi-res plots for all of the MR RF models
pdf("./BMI/Adipose_BMI_MR_RF_prediction_plots_for_Eurobats_and_METSIM.pdf",width = 18,height = 5)
par(mfrow=c(1,3))
plot(train_bmi$BMI~final_bmi_vip_rndFor$predicted,xlab="Predicted BMI",ylab="Actual BMI",main="Training set",xlim=c(15,50),ylim=c(15,50),pch = 20,col="purple")
abline(train_lm_bmi_vip)
plot(test_bmi$BMI~final_bmi_vip_rndFor$test$predicted,xlab="Predicted BMI",ylab="Actual BMI",main="Test set",xlim=c(15,50),ylim=c(15,50),pch = 20,col="purple")
abline(test_lm_bmi_vip)
plot(pheno$bmi~met_pred_bmi_vip,xlab="Predicted BMI",ylab="Actual BMI",main="METSIM",xlim=c(15,50),ylim=c(15,50),pch = 20,col="purple")
abline(met_lm_bmi_vip)
dev.off()

pdf("./WHR/Adipose_WHR_MR_RF_prediction_plots_for_Eurobats_and_METSIM.pdf",width = 18,height = 5)
par(mfrow=c(1,3))
plot(train_whr~final_whr_vip_rndFor$predicted,xlab="Predicted WHR",ylab="Actual WHR",main="Training set",pch = 20,col="blue")
abline(train_lm_whr_vip)
plot(test_whr~final_whr_vip_rndFor$test$predicted,xlab="Predicted WHR",ylab="Actual WHR",main="Test set",pch = 20,col="blue")
abline(test_lm_whr_vip)
plot(pheno$whr~met_pred_whr_vip,xlab="Predicted WHR",ylab="Actual WHR",main="METSIM",pch = 20,col="blue")
abline(met_lm_whr_vip)
dev.off()

pdf("./HOMA-IR/Adipose_HOMA-IR_MR_RF_prediction_plots_for_Eurobats_and_METSIM.pdf",width = 18,height = 5)
par(mfrow=c(1,3))
plot(train_HOMAIR~final_homair_vip_rndFor$predicted,xlab="Predicted ln(HOMA-IR)",ylab="Actual ln(HOMA-IR)",main="Training set",pch = 20,col="forestgreen")
abline(train_lm_homair_vip)
plot(test_HOMAIR~final_homair_vip_rndFor$test$predicted,xlab="Predicted ln(HOMA-IR)",ylab="Actual ln(HOMA-IR)",main="Test set",pch = 20,col="forestgreen")
abline(test_lm_homair_vip)
plot(pheno$homair~met_pred_homair_vip,xlab="Predicted ln(HOMA-IR)",ylab="Actual ln(HOMA-IR)",main="METSIM",pch = 20,col="forestgreen")
abline(met_lm_homair_vip)
dev.off()

pdf("./HDL/Adipose_HDL_MR_RF_prediction_plots_for_Eurobats_and_METSIM.pdf",width = 18,height = 5)
par(mfrow=c(1,3))
plot(train_HDL~final_hdl_vip_rndFor$predicted,xlab="Predicted HDL",ylab="Actual HDL",main="Training set",pch = 20,col="darkorange1")
abline(train_lm_hdl_vip)
plot(test_HDL~final_hdl_vip_rndFor$test$predicted,xlab="Predicted HDL",ylab="Actual HDL",main="Test set",pch = 20,col="darkorange1")
abline(test_lm_hdl_vip)
plot(pheno$s_hdlc~met_pred_hdl_vip,xlab="Predicted HDL",ylab="Actual HDL",main="METSIM",pch = 20,col="darkorange1")
abline(met_lm_hdl_vip)
dev.off()

pdf("./Triglycerides/Adipose_TriG_MR_RF_prediction_plots_for_Eurobats_and_METSIM.pdf",width = 18,height = 5)
par(mfrow=c(1,3))
plot(train_TriG~final_triG_vip_rndFor$predicted,xlab="Predicted Triglycerides",ylab="Actual Triglycerides",main="Training set",pch = 20,col="red")
abline(train_lm_triG_vip)
plot(test_TriG~final_triG_vip_rndFor$test$predicted,xlab="Predicted Triglycerides",ylab="Actual Triglycerides",main="Test set",pch = 20,col="red")
abline(test_lm_triG_vip)
plot(pheno$s_tottg~met_pred_trig_vip,xlab="Predicted Triglycerides",ylab="Actual Triglycerides",main="METSIM",pch = 20,col="red")
abline(met_lm_trig_vip)
dev.off()

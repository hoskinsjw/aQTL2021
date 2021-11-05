if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install(c("genefilter","preprocessCore","biomaRt","DESeq2","aracne.networks","viper"))
library(biomaRt)
library(DESeq2)
library(vsn)
library(hexbin)
library(grid)
library(gridExtra)
library(genefilter)
library(ggplot2)
library(foreach)
library(snowfall)
library(doParallel)
library(preprocessCore)
library(RNOmni)
library(viper)
library(aracne.networks)

setwd("YOUR WORKING DIRECTORY")

# Get table of Ensembl Gene IDs and Gene Symbols
ensembl=useMart("ENSEMBL_MART_ENSEMBL",host="http://grch37.ensembl.org", dataset = "hsapiens_gene_ensembl")
genes=getBM(c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position"),mart=ensembl)

# My modified write.regulon function that is way faster than the one in aracne.networks. However, the "n" argument may not 
# work. NOTE: this function requires the "aracne.networks" package to be loaded.
write.regulon3<-function(regulon,file="",sep="\t",header=TRUE,n=Inf,regulator=NULL,cpus=6,toScreen=F){
  fullTab<-data.frame()
  if(is.null(regulator)){
    cl=makeCluster(cpus)
    on.exit(stopCluster(cl))
    registerDoParallel(cl)
    tab=list()
    tab=foreach(tf=names(regulon)) %dopar%
      cbind(rep(tf,length(names(regulon[[tf]]$tfmode))),names(regulon[[tf]]$tfmode),regulon[[tf]]$tfmode,regulon[[tf]]$likelihood)
    fullTab=as.data.frame(do.call("rbind",tab)) # Very rapidly convert a list of data.frames to data.frame
  } else {
    tf<-regulator
    x<-regulon[[tf]]
    targets<-names(x$tfmode)
    moas<-x$tfmode
    likelihoods<-x$likelihood
    tab<-cbind(rep(tf,length(targets)),targets,moas,likelihoods)
    fullTab<-rbind(fullTab,tab)
  }
  colnames(fullTab)=c("Regulator","Target","MoA","likelihood")
  rownames(fullTab)=NULL
  
  if(file!=""){
    write.table(fullTab,file = file,sep = sep, col.names = header, row.names = F,quote = F)
  }
  
  if(toScreen){
    print(fullTab)
  }
  
  return(fullTab)
}

# This is a re-write of the PCAplot() function present within DEseq that allows you to look at more than PC1 and PC2,
# and use other expression matrices than DESeq objects. Also, you can give it a vector of values by which to shade
# the data points, which shows whether or not a covariate is evenly distributed in the PC space.
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

# Read in the data. The file locations are relative to your working directory, so adjust accordingly.
uk=read.table("./EUROBATS.F.annonymousIDs.rpkm",header = T,row.names = 1)

# Convert RPKM to TPM. See this paper for the reason this should be done. 
# (PMID: 22872506) http://diytranscriptomics.com/Reading/files/wagnerTPM.pdf
# The simplified formula I used below came from one of the creators of RSEM in a forum post here
# https://groups.google.com/forum/#!topic/rsem-users/W9RQrZIOzA4
# Note that the RPKM sum in the formula is per sample.
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
rpkmSum=colSums(uk)
temp=c()
tpm=foreach(i=1:dim(uk)[1]) %dopar% {
  for(j in 1:dim(uk)[2]){
    temp[j]=(uk[i,j]/rpkmSum[j])*(10^6)
  }
  temp
}
stopCluster(cl)
names(tpm)=rownames(uk)
tpm=structure(tpm, row.names = c(NA, -length(tpm[[1]])), class = "data.frame")
tpm=as.data.frame(t(tpm))
colnames(tpm)=colnames(uk)

# Calculate the average scaling factor (i.e. TPM/RPKM) across all samples
scales=colMeans(tpm)/colMeans(uk)
minScale=min(scales)
meanScale=mean(scales)
medianScale=median(scales)
maxScale=max(scales)

# Remove non-expressed genes, which I defined here as genes with fewer than 5% of samples having TPM>minScale.
exp_rows=rownames(tpm)
for(i in 1:dim(tpm)[1]){
  temp=c()
  for(j in 1:dim(tpm)[2]){
    temp[j]=tpm[i,j]>minScale
  }
  exp_rows[i]=(sum(temp)>38)
}
uk_exp=tpm[as.logical(exp_rows),]

# Let's add gene names as the first column
exp_genes=genes[match(rownames(uk_exp),genes[,1]),2]
uk_exp=cbind(exp_genes,uk_exp)
colnames(uk_exp)[1]="Gene_name"

# Get rid of genes with deprecated Ensembl IDs, which led to NA gene names
uk_noNA=subset(uk_exp,!is.na(Gene_name))

# There are genes that were given incorrect gene names based on their ENSEMBL Gene IDs, usually because they are
# paralogs, which leads to duplicate gene names. Since we are using gene names for the ARACNe analysis, this needs
# to be fixed. With the exception of ENSG00000255154 (which I rename from RPP14 to HTD2), these genes with duplicate
# names are actually nameless and uncharacterized, so I am removing them.
uk_noNA$Gene_name=as.character(uk_noNA$Gene_name)
uk_noNA[match("ENSG00000255154",rownames(uk_noNA)),1]="HTD2"
uk_noNA=uk_noNA[-which(duplicated(uk_noNA[,1])),]

# log2 transform the TPMs with a pseudocount of 1 added to avoid problems with TPM=0
logTPM=log2(uk_noNA[,2:767]+1)
rownames(logTPM)=uk_noNA$Gene_name

# Write it to file for ARACNe and downstream analyses
write.table(cbind("Gene"=rownames(logTPM),logTPM),"Eurobats_adipose_expressed_genes_logTPM.txt",sep = "\t",col.names = T,quote = F)

# Read in the un-normalized log2(TPMs) that were used to infer the ARACNe networks, and the covars file that has the samples 
# in the final order for the eQTL and aQTL analyses
adipo=read.table("Eurobats_adipose_expressed_genes_logTPM.txt",header = T, row.names = 1,sep = "\t")
adipo_set=ExpressionSet(assayData=as.matrix(adipo))
covars=as.data.frame(t(read.table("Filtered_Eurobats_adipose_covars_no_PEER.txt",sep="\t",header = T,row.names = 1)))

# The above covars file doesn't include BMI, so I also read in a version with BMI, and then add BMI to covars for PCA later
bmi=read.table("Eurobats_phenotypes_no_NA.txt",sep="\t",header = T,row.names = 1)
covars=cbind(covars[,1],bmi[na.omit(match(rownames(covars),rownames(bmi))),2],covars[,2:6])
colnames(covars)=c("Age","BMI","PC1","PC2","PC3","PC4","PC5")

### The methods section of the VIPER paper says to use expression levels that have been normalized to some
### reference samples, or to the average expression across all samples when clear reference samples are not available. 
### Therefore, below I transform the expression values to Z-scores prior to VIPER activity inference.

# Scale (i.e. Z-transform) genes in the full log2(TPM) dataset.
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
adipo_z=list()
adipo_z=foreach(i=1:dim(adipo)[1]) %dopar%
  as.numeric(scale(as.numeric(adipo[i,])))
adipo_z=as.data.frame(t(structure(adipo_z, row.names = c(NA, -length(adipo_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(adipo_z)=colnames(adipo)
rownames(adipo_z)=rownames(adipo)
z_set=ExpressionSet(assayData=as.matrix(adipo_z))

# Make new Eurobats adipose regulon from the ARACNE network
adipo_regulon=aracne2regulon("Eurobats_adipose_900bootstraps_ARACNe_network_with_LINC-PINT_no_header.txt",adipo_set,format="3col") # NOTE: The aracne network file cannot have a header!!!
#saveRDS(adipo_regulon,"Eurobats_adipose_logTPM_900boots_regulon.rds")
regTable=write.regulon3(adipo_regulon,file="Eurobats_adipose_900boots_regulon_with_LINC-PINT.txt")
regTable=read.table("Eurobats_adipose_900boots_regulon_with_LINC-PINT.txt",sep = "\t",header = T)

sum(!duplicated(regTable$Regulator)) # 4221 Regulators in the network

# Run VIPER to infer activity of regulators using the Eurobats adipose regulon. 
z_vip=viper(z_set,adipo_regulon)
zvip_results=exprs(z_vip)  # exprs() extracts the results matrix from the ExpressionSet object that is generated

# Filter viper activities and log2(TPM) to genes with inferred activity and samples that have genotype data and passed ancestry 
# filter in order according to the covars file.
filt_adipo=adipo[na.omit(match(rownames(zvip_results),rownames(adipo))),na.omit(match(rownames(covars),colnames(adipo)))]
filt_zvip=as.data.frame(zvip_results[,na.omit(match(rownames(covars),colnames(zvip_results)))])

# Let's also have a data.frame with all expressed genes, but only the samples that have genotype data and passed ancestry 
# filter for quantile normalization and inverse normal transformation prior to eQTL analyses of select loci.
geno_adipo=adipo[,na.omit(match(rownames(covars),colnames(adipo)))]

# Let's check to see how well correlated the activities are to their respective expression values.
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
expActCor=list()
startTime=Sys.time()
expActCor=foreach(i=1:dim(filt_adipo)[1]) %dopar%
  cor(t(filt_adipo)[,i],t(filt_zvip)[,i])
expActCor=unlist(expActCor) # Very rapidly convert a list to data.frame
stopCluster(cl)
Sys.time()-startTime
names(expActCor)=rownames(filt_adipo)

plot(density(expActCor),main="Activity-Expression Correlation",xlab="Pearson r")

# Quantile normalize the filtered and unfiltered log2(TPM) datasets, which both have 699 samples but differ on gene count.
qn_adipo=as.data.frame(normalize.quantiles(as.matrix(filt_adipo)))
rownames(qn_adipo)=rownames(filt_adipo)
colnames(qn_adipo)=colnames(filt_adipo)
qn_adipo_full=as.data.frame(normalize.quantiles(as.matrix(geno_adipo)))
rownames(qn_adipo_full)=rownames(geno_adipo)
colnames(qn_adipo_full)=colnames(geno_adipo)

# Inverse normal transform genes. Note that I'm doing this on the filtered and unfiltered log2(TPM) datasets, which both have 
# 699 samples but differ on gene count.
cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
int_adipo=list()
int_adipo=foreach(i=1:dim(qn_adipo)[1],.packages = 'RNOmni') %dopar%
  rankNorm(as.numeric(qn_adipo[i,]))
int_adipo=as.data.frame(t(structure(int_adipo, row.names = c(NA, -length(int_adipo[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(int_adipo)=colnames(qn_adipo)
rownames(int_adipo)=rownames(qn_adipo)

cl=makeCluster(6) # Note this was run on an 8 CPU machine, but the number of threads used could be reduced as needed
registerDoParallel(cl)
int_adipo_full=list()
int_adipo_full=foreach(i=1:dim(qn_adipo_full)[1],.packages = 'RNOmni') %dopar%
  rankNorm(as.numeric(qn_adipo_full[i,]))
int_adipo_full=as.data.frame(t(structure(int_adipo_full, row.names = c(NA, -length(int_adipo_full[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
stopCluster(cl)
colnames(int_adipo_full)=colnames(qn_adipo_full)
rownames(int_adipo_full)=rownames(qn_adipo_full)

# Let's make the map file for Matrix eQTL
# First, filter the genes data.frame from biomaRt to just those in uk_noNA based on Ensembl IDs
filt_genes=genes[na.omit(match(rownames(uk_noNA),genes$ensembl_gene_id)),]
map=filt_genes[,-1] 
rownames(map)=filt_genes$ensembl_gene_id 
colnames(map)=c("Gene","chr","s1","s2")
map[is.na(map$chr),] # all genes apparently matched
# Need the paste "chr" to all chromosomes
map$chr=sapply(map$chr,function(x) paste("chr",x,sep = ))

# Write the filtered, normalized and transformed data to file (in the correct sample order for eQTL/aQTL analyses), along with gene map files.
write.table(int_adipo,"Filtered_Eurobats_adipose_qnorm_INT_logTPMs_for_4213_regulators.txt",sep="\t",quote = F,row.names = T)
write.table(int_adipo_full,"Filtered_Eurobats_adipose_qnorm_INT_logTPMs_for_all_expressed_genes.txt",sep="\t",quote = F,row.names = T)
write.table(filt_zvip,"Filtered_Eurobats_adipose_unnormalized_activities_from_logTPM_for_4213_regulators.txt",sep="\t",quote = F,row.names = T)
write.table(map,"Hg19_gene_map_for_13776_expressed_genes_in_Eurobats_adipose.map",sep="\t",quote = F,row.names = F)


















### This script is for the co-localization analyses of the selected loci between BMI GWAS and cis and trans eQTLs/aQTLs.

#install.packages("backports")
#library(backports)
#install.packages("devtools")
#library(devtools)
#install_github("jrs95/hyprcoloc", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = F)
#browseVignettes("hyprcoloc") The install kept failing when trying to build the vignettes, so I disabled that.
#install_github("boxiangliu/locuscomparer")
library(hyprcoloc)
library(locuscomparer)
library(ggplot2)

setwd("YOUR WORKING DIRECTORY")

# Read in the LD matrices, and the GWAS and QTL data. The file locations are relative to your working directory, so adjust accordingly.
gwas=read.table("./Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt",sep = "\t",header = T)
cis_eqtl=read.table("./Eurobats_adipose_select_loci_cis-eQTLs_from_INT_logTPM.txt",sep = "\t",header = T)
cis_aqtl=read.table("./Eurobats_adipose_select_loci_cis-aQTLs_from_unnormalized_activities.txt",sep = "\t",header = T)
trans_eqtl=read.table("./Eurobats_adipose_select_loci_trans-eQTLs_for_BMI_MRs.txt",sep = "\t",header = T)
trans_aqtl=read.table("./Eurobats_adipose_select_loci_trans-aQTLs_for_BMI_MRs.txt",sep = "\t",header = T)
ld_files=read.table("./Select_loci_LD_matrix_file_list.txt",header=F)

ld=list()
index=1
for(i in ld_files$V1){
  ld[[index]]=read.table(paste("./",i,sep=""),sep = "\t",header = F)
  rownames(ld[[index]])=ld[[index]][,3]
  ld[[index]]=ld[[index]][,-c(1:5)]
  colnames(ld[[index]])=rownames(ld[[index]])
  index=index+1
}

# While preparing the manuscript and redoing the cis-QTL multiple testing analyses, I discovered a BMI locus previous
# missed (somehow) wherein we have an FDR sig cis-aQTL (for USP6) that is better than its eQTL. So I'm running co-localization 
# on this new locus (chr17p13.2) after the fact. Here's the data for the locus.
cis_eqtl17p13=read.table("./Eurobats_adipose_chr17p13.2_cis-eQTLs_from_INT_logTPM.txt",sep = "\t",header = T)
cis_aqtl17p13=read.table("./Eurobats_adipose_chr17p13.2_cis-aQTLs_from_unnormalized_activities.txt",sep = "\t",header = T)
trans_eqtl17p13=read.table("./Eurobats_adipose_chr17p13.2_trans-eQTLs_for_BMI_MRs.txt",sep = "\t",header = T)
trans_aqtl17p13=read.table("./Eurobats_adipose_chr17p13.2_trans-aQTLs_from_unnormalized_activities.txt",sep = "\t",header = T)
ld17p13=read.table("./Eurobats_chr17p13.2_LD_matrix.txt",sep = "\t",header = F)
rownames(ld17p13)=ld17p13[,3]
ld17p13=ld17p13[,-c(1:5)]
colnames(ld17p13)=rownames(ld17p13)

# Check that the SNPs included in each cis and trans QTL analysis are the same
all(unique(cis_eqtl$snps) %in% unique(cis_aqtl$snps)) # FALSE
all(unique(cis_eqtl$snps) %in% unique(trans_eqtl$snps)) # TRUE
all(unique(cis_aqtl$snps) %in% unique(trans_eqtl$snps)) # TRUE
all(unique(trans_eqtl$snps) %in% unique(trans_aqtl$snps)) # TRUE
# The cis-eQTL/aQTL SNPs do not match because their sets of tested genes are different and Matrix eQTL only tests cis-SNPs within 500kb. 
# Note also that both the cis-eQTL and cis-aQTL SNPs are proper subsets of the trans-QTL SNPs.
length(unique(cis_eqtl$snps)) # 78293
length(unique(cis_aqtl$snps)) # 74036
length(unique(trans_eqtl$snps)) # 79838
length(unique(trans_aqtl$snps)) # 79838

# Filter GWAS and QTL data to the same SNPs
filt_gwas=gwas[na.omit(match(unique(cis_aqtl$snps),gwas$SNP)),] # match to the cis-aQTL which has the fewest unique SNPs
filt_cis_eqtl=cis_eqtl[cis_eqtl$snps %in% filt_gwas$SNP,]
filt_cis_aqtl=cis_aqtl[cis_aqtl$snps %in% filt_gwas$SNP,]
filt_trans_eqtl=trans_eqtl[trans_eqtl$snps %in% filt_gwas$SNP,]
filt_trans_aqtl=trans_aqtl[trans_aqtl$snps %in% filt_gwas$SNP,]

filt_gwas17p13=gwas[na.omit(match(unique(cis_aqtl17p13$snps),gwas$SNP)),] # match to the cis-aQTL which has the fewest unique SNPs
filt_cis_eqtl17p13=cis_eqtl17p13[cis_eqtl17p13$snps %in% filt_gwas17p13$SNP,]
filt_cis_aqtl17p13=cis_aqtl17p13[cis_aqtl17p13$snps %in% filt_gwas17p13$SNP,]
filt_trans_eqtl17p13=trans_eqtl17p13[trans_eqtl17p13$snps %in% filt_gwas17p13$SNP,]
filt_trans_aqtl17p13=trans_aqtl17p13[trans_aqtl17p13$snps %in% filt_gwas17p13$SNP,]

# Calculate the SE of the Beta for the QTLs with the formula SE=beta/statistic
filt_cis_eqtl$SE=filt_cis_eqtl$beta/filt_cis_eqtl$statistic
filt_cis_aqtl$SE=filt_cis_aqtl$beta/filt_cis_aqtl$statistic
filt_trans_eqtl$SE=filt_trans_eqtl$beta/filt_trans_eqtl$statistic
filt_trans_aqtl$SE=filt_trans_aqtl$beta/filt_trans_aqtl$statistic

filt_cis_eqtl17p13$SE=filt_cis_eqtl17p13$beta/filt_cis_eqtl17p13$statistic
filt_cis_aqtl17p13$SE=filt_cis_aqtl17p13$beta/filt_cis_aqtl17p13$statistic
filt_trans_eqtl17p13$SE=filt_trans_eqtl17p13$beta/filt_trans_eqtl17p13$statistic
filt_trans_aqtl17p13$SE=filt_trans_aqtl17p13$beta/filt_trans_aqtl17p13$statistic

# Add chromosome and position to the QTLs for sorting
filt_cis_eqtl$chr=filt_gwas[match(filt_cis_eqtl$snps,filt_gwas$SNP),1]
filt_cis_aqtl$chr=filt_gwas[match(filt_cis_aqtl$snps,filt_gwas$SNP),1]
filt_trans_eqtl$chr=filt_gwas[match(filt_trans_eqtl$snps,filt_gwas$SNP),1]
filt_trans_aqtl$chr=filt_gwas[match(filt_trans_aqtl$snps,filt_gwas$SNP),1]
filt_cis_eqtl$position=filt_gwas[match(filt_cis_eqtl$snps,filt_gwas$SNP),2]
filt_cis_aqtl$position=filt_gwas[match(filt_cis_aqtl$snps,filt_gwas$SNP),2]
filt_trans_eqtl$position=filt_gwas[match(filt_trans_eqtl$snps,filt_gwas$SNP),2]
filt_trans_aqtl$position=filt_gwas[match(filt_trans_aqtl$snps,filt_gwas$SNP),2]

filt_cis_eqtl17p13$chr=filt_gwas17p13[match(filt_cis_eqtl17p13$snps,filt_gwas17p13$SNP),1]
filt_cis_aqtl17p13$chr=filt_gwas17p13[match(filt_cis_aqtl17p13$snps,filt_gwas17p13$SNP),1]
filt_trans_eqtl17p13$chr=filt_gwas17p13[match(filt_trans_eqtl17p13$snps,filt_gwas17p13$SNP),1]
filt_trans_aqtl17p13$chr=filt_gwas17p13[match(filt_trans_aqtl17p13$snps,filt_gwas17p13$SNP),1]
filt_cis_eqtl17p13$position=filt_gwas17p13[match(filt_cis_eqtl17p13$snps,filt_gwas17p13$SNP),2]
filt_cis_aqtl17p13$position=filt_gwas17p13[match(filt_cis_aqtl17p13$snps,filt_gwas17p13$SNP),2]
filt_trans_eqtl17p13$position=filt_gwas17p13[match(filt_trans_eqtl17p13$snps,filt_gwas17p13$SNP),2]
filt_trans_aqtl17p13$position=filt_gwas17p13[match(filt_trans_aqtl17p13$snps,filt_gwas17p13$SNP),2]

# Sort by chr and position
filt_gwas=filt_gwas[order(filt_gwas$CHR,filt_gwas$POS),]
filt_cis_eqtl=filt_cis_eqtl[order(filt_cis_eqtl$chr,filt_cis_eqtl$position),]
filt_cis_aqtl=filt_cis_aqtl[order(filt_cis_aqtl$chr,filt_cis_aqtl$position),]
filt_trans_eqtl=filt_trans_eqtl[order(filt_trans_eqtl$chr,filt_trans_eqtl$position),]
filt_trans_aqtl=filt_trans_aqtl[order(filt_trans_aqtl$chr,filt_trans_aqtl$position),]

filt_gwas17p13=filt_gwas17p13[order(filt_gwas17p13$CHR,filt_gwas17p13$POS),]
filt_cis_eqtl17p13=filt_cis_eqtl17p13[order(filt_cis_eqtl17p13$chr,filt_cis_eqtl17p13$position),]
filt_cis_aqtl17p13=filt_cis_aqtl17p13[order(filt_cis_aqtl17p13$chr,filt_cis_aqtl17p13$position),]
filt_trans_eqtl17p13=filt_trans_eqtl17p13[order(filt_trans_eqtl17p13$chr,filt_trans_eqtl17p13$position),]
filt_trans_aqtl17p13=filt_trans_aqtl17p13[order(filt_trans_aqtl17p13$chr,filt_trans_aqtl17p13$position),]

# Split GWAS and QTL Betas and SEs by locus, and also by genes for the QTL data. This section is pretty complex and hard to follow,
# for which I am both proud and embarassed.

# Betas
loci_betas=list()
for(i in 1:43){
  # Grab the genes per locus separately for cis-eQTL, cis-aQTL, trans-eQTL and trans-aQTL
  locus_genes=list()
  locus_genes[[1]]=as.character(unique(filt_cis_eqtl[filt_cis_eqtl$snps %in% rownames(ld[[i]]),"gene"]))
  locus_genes[[2]]=as.character(unique(filt_cis_aqtl[filt_cis_aqtl$snps %in% rownames(ld[[i]]),"gene"]))
  locus_genes[[3]]=as.character(unique(filt_trans_eqtl[filt_trans_eqtl$snps %in% rownames(ld[[i]]),"gene"]))
  locus_genes[[4]]=as.character(unique(filt_trans_aqtl[filt_trans_aqtl$snps %in% rownames(ld[[i]]),"gene"]))
  
  # Begin by grabbing all of the GWAS betas per locus
  loci_betas[[i]]=as.data.frame(filt_gwas[filt_gwas$SNP %in% rownames(ld[[i]]),"BETA"],
                                row.names = as.character(filt_gwas[filt_gwas$SNP %in% rownames(ld[[i]]),"SNP"]))
  colnames(loci_betas[[i]])="BMI"
  
  # Next we need to grab all of the QTL data per cis-e/aQTL and trans-e/aQTL
  for(j in 1:4){
    for(k in locus_genes[[j]]){
      # cis-eQTLs
      if(j==1){
        temp=filt_cis_eqtl[filt_cis_eqtl$gene %in% k,]
        loci_betas[[i]]=cbind(loci_betas[[i]],temp[match(rownames(loci_betas[[i]]),temp$snps),"beta"]) # Purposely allow for NAs from lack of match
        colnames(loci_betas[[i]])[dim(loci_betas[[i]])[2]]=paste("cis-e",k,sep="_") # This labels the new column appropriately
      }
      
      # cis-aQTLs
      if(j==2){
        temp=filt_cis_aqtl[filt_cis_aqtl$gene %in% k,]
        loci_betas[[i]]=cbind(loci_betas[[i]],temp[match(rownames(loci_betas[[i]]),temp$snps),"beta"]) # Purposely allow for NAs from lack of match
        colnames(loci_betas[[i]])[dim(loci_betas[[i]])[2]]=paste("cis-a",k,sep="_") # This labels the new column appropriately
      }
      
      # trans-eQTLs
      if(j==3){
        temp=filt_trans_eqtl[filt_trans_eqtl$gene %in% k,]
        loci_betas[[i]]=cbind(loci_betas[[i]],temp[match(rownames(loci_betas[[i]]),temp$snps),"beta"]) # Purposely allow for NAs from lack of match
        colnames(loci_betas[[i]])[dim(loci_betas[[i]])[2]]=paste("trans-e",k,sep="_") # This labels the new column appropriately
      }
      
      # trans-aQTLs
      if(j==4){
        temp=filt_trans_aqtl[filt_trans_aqtl$gene %in% k,]
        loci_betas[[i]]=cbind(loci_betas[[i]],temp[match(rownames(loci_betas[[i]]),temp$snps),"beta"]) # Purposely allow for NAs from lack of match
        colnames(loci_betas[[i]])[dim(loci_betas[[i]])[2]]=paste("trans-a",k,sep="_") # This labels the new column appropriately
      }
    }
  }
}

loci_betas17p13=list()
for(i in 1){
  # Grab the genes per locus separately for cis-eQTL, cis-aQTL, trans-eQTL and trans-aQTL
  locus_genes=list()
  locus_genes[[1]]=as.character(unique(filt_cis_eqtl17p13[filt_cis_eqtl17p13$snps %in% rownames(ld17p13),"gene"]))
  locus_genes[[2]]=as.character(unique(filt_cis_aqtl17p13[filt_cis_aqtl17p13$snps %in% rownames(ld17p13),"gene"]))
  locus_genes[[3]]=as.character(unique(filt_trans_eqtl17p13[filt_trans_eqtl17p13$snps %in% rownames(ld17p13),"gene"]))
  locus_genes[[4]]=as.character(unique(filt_trans_aqtl17p13[filt_trans_aqtl17p13$snps %in% rownames(ld17p13),"gene"]))
  
  # Begin by grabbing all of the GWAS betas per locus
  loci_betas17p13[[i]]=as.data.frame(filt_gwas17p13[filt_gwas17p13$SNP %in% rownames(ld17p13),"BETA"],
                                row.names = as.character(filt_gwas17p13[filt_gwas17p13$SNP %in% rownames(ld17p13),"SNP"]))
  colnames(loci_betas17p13[[i]])="BMI"
  
  # Next we need to grab all of the QTL data per cis-e/aQTL and trans-e/aQTL
  for(j in 1:4){
    for(k in locus_genes[[j]]){
      # cis-eQTLs
      if(j==1){
        temp=filt_cis_eqtl17p13[filt_cis_eqtl17p13$gene %in% k,]
        loci_betas17p13[[i]]=cbind(loci_betas17p13[[i]],temp[match(rownames(loci_betas17p13[[i]]),temp$snps),"beta"]) # Purposely allow for NAs from lack of match
        colnames(loci_betas17p13[[i]])[dim(loci_betas17p13[[i]])[2]]=paste("cis-e",k,sep="_") # This labels the new column appropriately
      }
      
      # cis-aQTLs
      if(j==2){
        temp=filt_cis_aqtl17p13[filt_cis_aqtl17p13$gene %in% k,]
        loci_betas17p13[[i]]=cbind(loci_betas17p13[[i]],temp[match(rownames(loci_betas17p13[[i]]),temp$snps),"beta"]) # Purposely allow for NAs from lack of match
        colnames(loci_betas17p13[[i]])[dim(loci_betas17p13[[i]])[2]]=paste("cis-a",k,sep="_") # This labels the new column appropriately
      }
      
      # trans-eQTLs
      if(j==3){
        temp=filt_trans_eqtl17p13[filt_trans_eqtl17p13$gene %in% k,]
        loci_betas17p13[[i]]=cbind(loci_betas17p13[[i]],temp[match(rownames(loci_betas17p13[[i]]),temp$snps),"beta"]) # Purposely allow for NAs from lack of match
        colnames(loci_betas17p13[[i]])[dim(loci_betas17p13[[i]])[2]]=paste("trans-e",k,sep="_") # This labels the new column appropriately
      }
      
      # trans-aQTLs
      if(j==4){
        temp=filt_trans_aqtl17p13[filt_trans_aqtl17p13$gene %in% k,]
        loci_betas17p13[[i]]=cbind(loci_betas17p13[[i]],temp[match(rownames(loci_betas17p13[[i]]),temp$snps),"beta"]) # Purposely allow for NAs from lack of match
        colnames(loci_betas17p13[[i]])[dim(loci_betas17p13[[i]])[2]]=paste("trans-a",k,sep="_") # This labels the new column appropriately
      }
    }
  }
}

# SEs
loci_ses=list()
for(i in 1:43){
  # Grab the genes per locus separately for cis-eQTL, cis-aQTL, trans-eQTL and trans-aQTL
  locus_genes=list()
  locus_genes[[1]]=as.character(unique(filt_cis_eqtl[filt_cis_eqtl$snps %in% rownames(ld[[i]]),"gene"]))
  locus_genes[[2]]=as.character(unique(filt_cis_aqtl[filt_cis_aqtl$snps %in% rownames(ld[[i]]),"gene"]))
  locus_genes[[3]]=as.character(unique(filt_trans_eqtl[filt_trans_eqtl$snps %in% rownames(ld[[i]]),"gene"]))
  locus_genes[[4]]=as.character(unique(filt_trans_aqtl[filt_trans_aqtl$snps %in% rownames(ld[[i]]),"gene"]))
  
  # Begin by grabbing all of the GWAS ses per locus
  loci_ses[[i]]=as.data.frame(filt_gwas[filt_gwas$SNP %in% rownames(ld[[i]]),"SE"],
                              row.names = as.character(filt_gwas[filt_gwas$SNP %in% rownames(ld[[i]]),"SNP"]))
  colnames(loci_ses[[i]])="BMI"
  
  # Next we need to grab all of the QTL data per cis-e/aQTL and trans-e/aQTL
  for(j in 1:4){
    for(k in locus_genes[[j]]){
      # cis-eQTLs
      if(j==1){
        temp=filt_cis_eqtl[filt_cis_eqtl$gene %in% k,]
        loci_ses[[i]]=cbind(loci_ses[[i]],temp[match(rownames(loci_ses[[i]]),temp$snps),"SE"]) # Purposely allow for NAs from lack of match
        colnames(loci_ses[[i]])[dim(loci_ses[[i]])[2]]=paste("cis-e",k,sep="_") # This labels the new column appropriately
      }
      
      # cis-aQTLs
      if(j==2){
        temp=filt_cis_aqtl[filt_cis_aqtl$gene %in% k,]
        loci_ses[[i]]=cbind(loci_ses[[i]],temp[match(rownames(loci_ses[[i]]),temp$snps),"SE"]) # Purposely allow for NAs from lack of match
        colnames(loci_ses[[i]])[dim(loci_ses[[i]])[2]]=paste("cis-a",k,sep="_") # This labels the new column appropriately
      }
      
      # trans-eQTLs
      if(j==3){
        temp=filt_trans_eqtl[filt_trans_eqtl$gene %in% k,]
        loci_ses[[i]]=cbind(loci_ses[[i]],temp[match(rownames(loci_ses[[i]]),temp$snps),"SE"]) # Purposely allow for NAs from lack of match
        colnames(loci_ses[[i]])[dim(loci_ses[[i]])[2]]=paste("trans-e",k,sep="_") # This labels the new column appropriately
      }
      
      # trans-aQTLs
      if(j==4){
        temp=filt_trans_aqtl[filt_trans_aqtl$gene %in% k,]
        loci_ses[[i]]=cbind(loci_ses[[i]],temp[match(rownames(loci_ses[[i]]),temp$snps),"SE"]) # Purposely allow for NAs from lack of match
        colnames(loci_ses[[i]])[dim(loci_ses[[i]])[2]]=paste("trans-a",k,sep="_") # This labels the new column appropriately
      }
    }
  }
}

loci_ses17p13=list()
for(i in 1){
  # Grab the genes per locus separately for cis-eQTL, cis-aQTL, trans-eQTL and trans-aQTL
  locus_genes=list()
  locus_genes[[1]]=as.character(unique(filt_cis_eqtl17p13[filt_cis_eqtl17p13$snps %in% rownames(ld17p13),"gene"]))
  locus_genes[[2]]=as.character(unique(filt_cis_aqtl17p13[filt_cis_aqtl17p13$snps %in% rownames(ld17p13),"gene"]))
  locus_genes[[3]]=as.character(unique(filt_trans_eqtl17p13[filt_trans_eqtl17p13$snps %in% rownames(ld17p13),"gene"]))
  locus_genes[[4]]=as.character(unique(filt_trans_aqtl17p13[filt_trans_aqtl17p13$snps %in% rownames(ld17p13),"gene"]))
  
  # Begin by grabbing all of the GWAS ses per locus
  loci_ses17p13[[i]]=as.data.frame(filt_gwas17p13[filt_gwas17p13$SNP %in% rownames(ld17p13),"SE"],
                                     row.names = as.character(filt_gwas17p13[filt_gwas17p13$SNP %in% rownames(ld17p13),"SNP"]))
  colnames(loci_ses17p13[[i]])="BMI"
  
  # Next we need to grab all of the QTL data per cis-e/aQTL and trans-e/aQTL
  for(j in 1:4){
    for(k in locus_genes[[j]]){
      # cis-eQTLs
      if(j==1){
        temp=filt_cis_eqtl17p13[filt_cis_eqtl17p13$gene %in% k,]
        loci_ses17p13[[i]]=cbind(loci_ses17p13[[i]],temp[match(rownames(loci_ses17p13[[i]]),temp$snps),"SE"]) # Purposely allow for NAs from lack of match
        colnames(loci_ses17p13[[i]])[dim(loci_ses17p13[[i]])[2]]=paste("cis-e",k,sep="_") # This labels the new column appropriately
      }
      
      # cis-aQTLs
      if(j==2){
        temp=filt_cis_aqtl17p13[filt_cis_aqtl17p13$gene %in% k,]
        loci_ses17p13[[i]]=cbind(loci_ses17p13[[i]],temp[match(rownames(loci_ses17p13[[i]]),temp$snps),"SE"]) # Purposely allow for NAs from lack of match
        colnames(loci_ses17p13[[i]])[dim(loci_ses17p13[[i]])[2]]=paste("cis-a",k,sep="_") # This labels the new column appropriately
      }
      
      # trans-eQTLs
      if(j==3){
        temp=filt_trans_eqtl17p13[filt_trans_eqtl17p13$gene %in% k,]
        loci_ses17p13[[i]]=cbind(loci_ses17p13[[i]],temp[match(rownames(loci_ses17p13[[i]]),temp$snps),"SE"]) # Purposely allow for NAs from lack of match
        colnames(loci_ses17p13[[i]])[dim(loci_ses17p13[[i]])[2]]=paste("trans-e",k,sep="_") # This labels the new column appropriately
      }
      
      # trans-aQTLs
      if(j==4){
        temp=filt_trans_aqtl17p13[filt_trans_aqtl17p13$gene %in% k,]
        loci_ses17p13[[i]]=cbind(loci_ses17p13[[i]],temp[match(rownames(loci_ses17p13[[i]]),temp$snps),"SE"]) # Purposely allow for NAs from lack of match
        colnames(loci_ses17p13[[i]])[dim(loci_ses17p13[[i]])[2]]=paste("trans-a",k,sep="_") # This labels the new column appropriately
      }
    }
  }
}

# The SNPs and/or genes with NAs need to be filtered out. I've decided to discard genes with >20% of locus SNPs as NAs. This is followed by filtering any 
# remaining SNPs that still contain any NAs, though this is typically a very small number after the gene filtering. This assures we remove all NAs, while
# keeping the vast majority of SNPs and probably all genes of main interest.

# Keep genes with <20% NAs for each locus
filt_betas=list()
filt_ses=list()
for(i in 1:43){
  gene_nas=c()
  for(j in colnames(loci_betas[[i]])){
    gene_nas[j]=sum(is.na(loci_betas[[i]][,j]))
  }
  filt_betas[[i]]=loci_betas[[i]][,gene_nas<(0.2*dim(loci_betas[[i]])[1])]
  filt_ses[[i]]=loci_ses[[i]][,gene_nas<(0.2*dim(loci_ses[[i]])[1])]
}

filt_betas17p13=list()
filt_ses17p13=list()
for(i in 1){
  gene_nas=c()
  for(j in colnames(loci_betas17p13[[i]])){
    gene_nas[j]=sum(is.na(loci_betas17p13[[i]][,j]))
  }
  filt_betas17p13[[i]]=loci_betas17p13[[i]][,gene_nas<(0.2*dim(loci_betas17p13[[i]])[1])]
  filt_ses17p13[[i]]=loci_ses17p13[[i]][,gene_nas<(0.2*dim(loci_ses17p13[[i]])[1])]
}

# Filter out SNPs with any remaining NAs
for(i in 1:43){
  snp_nas=c()
  for(j in rownames(filt_betas[[i]])){
    snp_nas[j]=sum(is.na(filt_betas[[i]][j,]))
  }
  filt_betas[[i]]=filt_betas[[i]][snp_nas==0,]
  filt_ses[[i]]=filt_ses[[i]][snp_nas==0,]
}

for(i in 1){
  snp_nas=c()
  for(j in rownames(filt_betas17p13[[i]])){
    snp_nas[j]=sum(is.na(filt_betas17p13[[i]][j,]))
  }
  filt_betas17p13[[i]]=filt_betas17p13[[i]][snp_nas==0,]
  filt_ses17p13[[i]]=filt_ses17p13[[i]][snp_nas==0,]
}

# Filter the LD matrices to include only the SNPs remaining and put them in a list
LD=list()
for(i in 1:length(ld)){
  LD[[i]]=ld[[i]][rownames(filt_betas[[i]]),rownames(filt_betas[[i]])]
}

LD17p13=list()
for(i in 1){
  LD17p13[[i]]=ld17p13[rownames(filt_betas17p13[[i]]),rownames(filt_betas17p13[[i]])]
}

# It seems that HyPrColoc's assumption of a single causal variant per locus can prove problematic when assessing
# multiple traits that may have different information paths (i.e. differential functional variants). For example:
testZNF436=hyprcoloc(as.matrix(filt_betas[[1]]),as.matrix(filt_ses[[1]]),
                     trait.names=colnames(filt_betas[[1]]),snp.id=rownames(filt_betas[[1]]),ld.matrix = LD[[1]],
                     trait.subset = c("BMI","cis-a_ZNF436"),snpscores = T)
test_eEPHB2=hyprcoloc(as.matrix(filt_betas[[1]]),as.matrix(filt_ses[[1]]),
                     trait.names=colnames(filt_betas[[1]]),snp.id=rownames(filt_betas[[1]]),ld.matrix = LD[[1]],
                     trait.subset = c("BMI","cis-e_EPHB2"),snpscores = T)
test_aEPHB2=hyprcoloc(as.matrix(filt_betas[[1]]),as.matrix(filt_ses[[1]]),
                    trait.names=colnames(filt_betas[[1]]),snp.id=rownames(filt_betas[[1]]),ld.matrix = LD[[1]],
                    trait.subset = c("BMI","cis-a_EPHB2"),snpscores = T)
test_bothEPHB2=hyprcoloc(as.matrix(filt_betas[[1]]),as.matrix(filt_ses[[1]]),
                      trait.names=colnames(filt_betas[[1]]),snp.id=rownames(filt_betas[[1]]),ld.matrix = LD[[1]],
                      trait.subset = c("BMI","cis-e_EPHB2","cis-a_EPHB2"),snpscores = T)
eEPHB2_aEPHB2=hyprcoloc(as.matrix(filt_betas[[1]]),as.matrix(filt_ses[[1]]),
                         trait.names=colnames(filt_betas[[1]]),snp.id=rownames(filt_betas[[1]]),ld.matrix = LD[[1]],
                         trait.subset = c("cis-e_EPHB2","cis-a_EPHB2"),snpscores = T)
testEPHB2_DAPK2=hyprcoloc(as.matrix(filt_betas[[1]]),as.matrix(filt_ses[[1]]),
                          trait.names=colnames(filt_betas[[1]]),snp.id=rownames(filt_betas[[1]]),ld.matrix = LD[[1]],
                          trait.subset = c("BMI","cis-a_ZNF436","cis-a_EPHB2","trans-a_DAPK2"),snpscores = T)
testAll=hyprcoloc(as.matrix(filt_betas[[1]]),as.matrix(filt_ses[[1]]),
                          trait.names=colnames(filt_betas[[1]]),snp.id=rownames(filt_betas[[1]]),ld.matrix = LD[[1]],
                          snpscores = T)
strsplit(testAll[[1]][1,2],", ")[[1]]
# BMI and the EPHB2 cis-aQTL have a PP of 0.8566 with 0.0834 of the PP explained by rs6692586 (the top GWAS variant). The EPHB2
# cis-eQTL does not co-localize with either, which I suspect is due to the cis-eQTL being discarded due to smaller effect sizes.
# For BMI, EPHB2 cis-aQTL and DAPK2 trans-aQTL, the three co-localize with a PP of 0.535 and 0.1363 was explained by rs2235549.
# When considering all traits, there is one large cluster with 0/16 cis-eQTLs, 1/8 cis-aQTLs, 52/100 trans-eQTLs, 85/100 trans-aQTLs 
# (WITHOUT BMI) with a 0.4108 PP of a common functional variant and essentially all of that PP explained by rs309527 (which has a BMI P ~5 orders
# of magnitude worse than the top GWAS variants). There are two other smaller clusters, one composed of 5 trans-eQTLs and 2 trans-aQTLs
# (PP=0.4689) and another cluster with cis-e_TCEA3 and trans-e_PITX2 (PP=0.5626), neither of which clustered with BMI.All together, the 
# co-localization of EPHB2 cis-aQTL with the trans-QTLs took priority over that of EPHB2 with BMI, despite the high pairwise PP. So, 
# the best variant for the co-localization varies by the traits considered. It seems that if there is evidence of co-localization among 
# the MR trans-QTLs, these co-localizations will swamp out any other co-localizations with the GWAS trait unless there is a separate cluster 
# identified. For these reasons, I think it wise to first run HyPrColoc pairwise between GWAS trait and QTLs. Then groups of traits may be 
# analyzed together in a more targeted fashion.

# Pairwise co-localization, which is essentially equivalent with using coloc.
# Make a data.frame with only the QTLs co-localized with BMI.
GWAScoloc=list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43) # initialize a list of 43 objects
posColoc=data.frame() # initialize a data.frame for the positive co-localized traits
for(i in 1:43){
  # Grab the traits list for each locus
  traits=colnames(filt_betas[[i]])

  # Iterate through each of the QTLs. Exclude BMI from the index of traits since all will be tested against BMI.
  for(j in traits[-1]){
    GWAScoloc[[i]][[j]]=hyprcoloc(as.matrix(filt_betas[[i]]),as.matrix(filt_ses[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas[[i]]),ld.matrix = LD[[i]],
                                  trait.subset = c("BMI",j))
    # Grab only the QTLs with co-localization into posColoc
    if(!is.na(GWAScoloc[[i]][[j]]$results[1,3])){
      posColoc=rbind(posColoc,cbind(GWAScoloc[[i]][[j]]$results,"locus"=i))
    }
  }
}

GWAScoloc17p13=list(1) # initialize a list of 1 object
posColoc17p13=data.frame() # initialize a data.frame for the positive co-localized traits
for(i in 1){
  # Grab the traits list for each locus
  traits=colnames(filt_betas17p13[[i]])
  
  # Iterate through each of the QTLs. Exclude BMI from the index of traits since all will be tested against BMI.
  for(j in traits[-1]){
    GWAScoloc17p13[[i]][[j]]=hyprcoloc(as.matrix(filt_betas17p13[[i]]),as.matrix(filt_ses17p13[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas17p13[[i]]),ld.matrix = LD17p13[[i]],
                                  trait.subset = c("BMI",j))
    # Grab only the QTLs with co-localization into posColoc
    if(!is.na(GWAScoloc17p13[[i]][[j]]$results[1,3])){
      posColoc17p13=rbind(posColoc17p13,cbind(GWAScoloc17p13[[i]][[j]]$results,"locus"=i))
    }
  }
}

# Let's change the locus number to legit chromosomal locus ID
loci=read.table(".Select_loci_list.txt",sep = "\t",header = T)
posColoc$locus=loci[match(posColoc$locus,rownames(loci)),1]

posColoc17p13$locus=rep("17p13.2",dim(posColoc17p13)[1])
# Only trans-e_MKX pairwise colocalizes with the BMI GWAS signal at 17p13.2. This locus was of interest because it had
# an FDR sig cis-aQTL with USP6 that was better than its eQTL. However, neither the eQTL nor the aQTL signals for USP6 
# colocalize with the BMI GWAS signal at the locus. I'll just add the trans-e_MKX colocalization results to the rest
# manually.

# Let's tally the co-localized QTLs per locus
pair_tallies=data.frame(row.names = loci$Locus)
for(i in loci$Locus){
  pair_tallies[i,1]=loci[loci$Locus==i,2]
  pair_tallies[i,2]=dim(posColoc[posColoc$locus==i,])[1]
}
colnames(pair_tallies)=c("GWAS","Colocalizations with BMI")
# HyPrColoc reports colocalizations with PP>0.25, which is what is tallied here.
# Locus                                    GWAS                      Colocalizations with BMI
# 1p36.1                                   BMI                       41
# 1p13.3                                   HDL                        1
# 1q24_loc1                                BMI                       56
# 1q24_loc2                                WHR                        0
# 1q25.2                                   BMI                        0
# 1q41                                     WHR                        0
# 2p24                           Triglycerides                        0
# 2p23.3                                   BMI                        2
# 2p21                                     T2D                        0
# 2q24.3                                   WHR                        1
# 3p25                                     HDL                        0
# 3q27                                     T2D                        0
# 4q22_loc1                                BMI                        2
# 4q22_loc2                                BMI                        0
# 5q11.2                           HDL and T2D                        4
# 5q31.2                                   BMI                       28
# 6p21_loc1                      Triglycerides                        3
# 6p21_loc2                                HDL                        2
# 6q22.3                                   WHR                        0
# 7q32         BMI, T2D, HDL and Triglycerides                        7
# 7q36                                     HDL                        0
# 8q21.1                                   BMI                        0
# 8q21.2                                   BMI                        0
# 10p13                                    BMI                        0
# 10q22.2                                  BMI                        1
# 10q24.2                                  BMI                        0
# 10q26.3                                  BMI                       15
# 11p11.2                                  HDL                        1
# 11q13                                    BMI                        0
# 12p13.33                                 BMI                        0
# 12p13.1                                  BMI                       81
# 12q13.13                                 BMI                        0
# 15q21.3                                  HDL                        0
# 15q24.1                                  BMI                        0
# 17q21.2                                  BMI                        0
# 17q25.3                                  HDL                        1
# 18q21.1                                  HDL                        0
# 18q21.3                                  BMI                        0
# 19p13.1                T2D and Triglycerides                        0
# 19q13.3_loc1                             HDL                        0
# 19q13.3_loc2                             BMI                        0
# 20q13.32                                 T2D                        0
# 22q13.3                                  T2D                        0

# Let's write the pairwise analysis to file for reference in the rest of these analyses.
write.table(posColoc,"Pairwise_HyPrColoc_between_BMI_and_each_QTL_for_select_loci.txt",sep = "\t",row.names = F,quote = F)

# Given the importance of 1p36 to this study, and the fact that there is a co-localized cis-aQTL and many co-localized
# trans-QTLs, I need to run the pairwise co-localizations between cis-a_EPHB2 and all trans-QTLs that co-localize
# with the BMI signal.
bmi_traits=gsub("BMI, ","",posColoc[posColoc$locus=="1p36.1",2])
bmi_traits=bmi_traits[-1] # Get rid of cis-a_EPHB2 for the iterations since it is the constant for this co-localization.
ephb2_coloc=data.frame()
for(i in bmi_traits){
  temp=hyprcoloc(as.matrix(filt_betas[[1]]),as.matrix(filt_ses[[1]]),
                trait.names=colnames(filt_betas[[1]]),snp.id=rownames(filt_betas[[1]]),ld.matrix = LD[[1]],
                trait.subset = c("cis-a_EPHB2",i),snpscores = T)
  ephb2_coloc=rbind(ephb2_coloc,temp$results)
}
# All of the trans-QTLs that co-localized with BMI also colocalize with cis-a_EPHB2, though the top candidate SNP varied between 
# rs4654828, rs309527 and rs6426776. The most common candidate variant was rs4654828. The vast majority have PP>0.7.

### Mult-trait analyses of loci with more than one QTL colocalizing with BMI.

### 1p36
# Start with all traits included
locus1all=hyprcoloc(as.matrix(filt_betas[[1]]),as.matrix(filt_ses[[1]]),
                  trait.names=colnames(filt_betas[[1]]),snp.id=rownames(filt_betas[[1]]),ld.matrix = LD[[1]],
                  snpscores = T)
# This is the same analysis as performed above and stored in testAll. See description of results above.

# Next let's look at the QTLs co-localized best with rs6692586
rs6692586_traits=c("BMI",gsub("BMI, ","",posColoc[posColoc$candidate_snp=="rs6692586",2]))
locus1rs6692586=hyprcoloc(as.matrix(filt_betas[[1]]),as.matrix(filt_ses[[1]]),
                    trait.names=colnames(filt_betas[[1]]),snp.id=rownames(filt_betas[[1]]),ld.matrix = LD[[1]],
                    trait.subset = rs6692586_traits,snpscores = T)
# The co-localized cluster of traits doesn't include BMI. I tried increasing the number of the QTLs included with BMI
# stepwise, and BMI remains in the cluster until the first 5 QTLs are included in the analysis. Then it disappears again
# and the top candidate SNP changes from rs6692586 to rs309527.

# Since EPHB2 seems to be the cis gene mediating much of this BMI signal, and since including all genes leads to the exclusion
# of BMI from the main cluster, let's try 3-way co-localization with BMI, EPHB2 and each other QTL.
bmi_ephb2_3way=data.frame()
for(i in bmi_traits){
  temp=hyprcoloc(as.matrix(filt_betas[[1]]),as.matrix(filt_ses[[1]]),
                 trait.names=colnames(filt_betas[[1]]),snp.id=rownames(filt_betas[[1]]),ld.matrix = LD[[1]],
                 trait.subset = c("BMI","cis-a_EPHB2",i),snpscores = T)
  if(temp[[1]]$traits!="BMI, cis-a_EPHB2"){ # exclude clusters that are only BMI and cis-a_EPHB2
    bmi_ephb2_3way=rbind(bmi_ephb2_3way,temp$results)
  }
}
# All three-way colocalized clusters with BMI, cis-a_EPHB2 and another trans-QTL had PP between 0.3863 and 0.7210 and rs2235549
# as the top candidate SNP. There are also 3 clusters of only cis-a_EPHB2 and a trans-QTL without BMI with rs4654828 as top
# candidate variant and PP>0.8.

# Let's run an analysis with only the QTLs that co-localized with BMI in the pairwise analysis.
coloc1_traits=c("BMI",gsub("BMI, ","",posColoc[posColoc$locus=="1p36.1",2]))
locus1coloc1=hyprcoloc(as.matrix(filt_betas[[1]]),as.matrix(filt_ses[[1]]),
                          trait.names=colnames(filt_betas[[1]]),snp.id=rownames(filt_betas[[1]]),ld.matrix = LD[[1]],
                          trait.subset = coloc1_traits,snpscores = T)
locus1coloc_cluster1=strsplit(locus1coloc1[[1]][1,2],", ")[[1]]
# As with locus1all, the only cluster does not include BMI, but instead has all 41 QTLs including cis-a_EPHB2. The cluster 
# of QTLs has PP=0.8512 with rs4654828 as the top candidate variant explaining essentially all of the PP.

# Let's run an analysis with only the QTLs that were FDR significant in the original eQTL/aQTL analyses.
fdr1_traits=c("BMI","trans-a_ABCC9","trans-a_ADH1A","trans-a_ADRB1","trans-a_AKAP12","trans-a_ARHGEF26","trans-a_BMP3","trans-a_BTG2",
              "trans-a_CECR2","trans-a_CHAF1B","trans-a_CSNK2A2","trans-a_DAPK2","trans-a_DDR1","trans-a_DMRT2","trans-a_EIF4EBP2",
              "trans-a_ESR1","trans-a_GHR","trans-a_GLRA1","trans-a_GNA14","trans-a_GNAI1","trans-a_GPHN","trans-a_HLF","trans-a_HOMEZ",
              "trans-a_LASP1","trans-a_MAP3K5","trans-a_MET","trans-a_MKX","trans-a_NCAM1","trans-a_NTRK3","trans-a_P2RX6",
              "trans-a_PDLIM1","trans-a_PRKAG2","trans-a_RASL10B","trans-a_RASSF4","trans-a_RGS3","trans-a_RORB","trans-a_RRBP1",
              "trans-a_RTN1","trans-a_SIRT1","trans-a_STMN2","trans-a_TBX4","trans-a_ZNF334","cis-e_EPHB2","cis-e_KDM1A","cis-e_PNRC2",
              "cis-e_TCEA3","cis-a_EPHB2","cis-a_ZNF436")
locus1fdr=hyprcoloc(as.matrix(filt_betas[[1]]),as.matrix(filt_ses[[1]]),
                    trait.names=colnames(filt_betas[[1]]),snp.id=rownames(filt_betas[[1]]),ld.matrix = LD[[1]],
                    trait.subset = fdr1_traits,snpscores = T)
# As with the other multi-trait analyses of this locus, inclusion of many trans-QTLs ends up excluding BMI while all of the
# FDR significant trans-QTLs and cis-a_EPHB2 cluster together with a PP=0.6782 with rs4654828 as the top candidate SNP (GWAS
# P=3.3E-10). This is several orders of magnitude weaker than the top GWAS SNP of the locus (rs6692586 P=1.1E-16) and has an
# r2 with the top SNP of about 0.34. However, it is notable that while the BMI GWAS P is quite different between these two SNPs
# the EPHB2 aQTL P is very similar. This is not the case for the trans-aQTLs where the P values differ by 1-2 orders of
# magnitude and most are not FDR significant for rs6692586.

### I think the take home message is that this 1p36 BMI locus is complex and likely features multiple functional variants
### whose info each flows a little differently through the gene regulatory network. However, it is clear that the EPHB2
### cis-aQTL co-localizes pretty well with the BMI locus.

### 1q24_loc1 (centered on BMI signal), which is the 3rd locus in the list.
# Start with all traits included
locus3all=hyprcoloc(as.matrix(filt_betas[[3]]),as.matrix(filt_ses[[3]]),
                    trait.names=colnames(filt_betas[[3]]),snp.id=rownames(filt_betas[[3]]),ld.matrix = LD[[3]],
                    snpscores = T)
locus3all_cluster1=strsplit(locus3all[[1]][1,2],", ")[[1]]
# There is one cluster with a PP=0.4547 and rs4916229 as top candidate variant that includes BMI, cis-e_PRRC2C, 18 trans-eQTLs,
# and 58 trans-aQTLs. This is more than the pairwise colocalizations (there were 56 QTLs that colocalized with BMI).

# Let's run an analysis with only the QTLs that co-localized with BMI in the pairwise analysis.
coloc3_traits=c("BMI",gsub("BMI, ","",posColoc[posColoc$locus=="1q24_loc1",2]))
locus3coloc=hyprcoloc(as.matrix(filt_betas[[3]]),as.matrix(filt_ses[[3]]),
                       trait.names=colnames(filt_betas[[3]]),snp.id=rownames(filt_betas[[3]]),ld.matrix = LD[[3]],
                       trait.subset = coloc3_traits,snpscores = T)
locus3coloc_cluster1=strsplit(locus3coloc[[1]][1,2],", ")[[1]]
# By restricting the multi-trait analysis to the pairwise colocalized QTLs, the PP goes up to 0.7195 with rs4916229 again as
# the top candidate variant explaining essentially all of the PP. That SNP is in extremely high LD with the top BMI SNP (r2=0.987).

# Let's run an analysis with only the QTLs that were FDR significant in the original eQTL/aQTL analyses.
fdr3_traits=c("BMI","trans-a_EPHB2","trans-a_P2RX6","trans-a_RASSF4","trans-a_STMN2")
locus3fdr=hyprcoloc(as.matrix(filt_betas[[3]]),as.matrix(filt_ses[[3]]),
                      trait.names=colnames(filt_betas[[3]]),snp.id=rownames(filt_betas[[3]]),ld.matrix = LD[[3]],
                      trait.subset = fdr3_traits,snpscores = T)
# The colocalization of all of these FDR significant trans-aQTLs is rather remarkable with a PP=0.9701.

# Let's run the same analysis with the FDR sig tran-aQTLs, but also add the cis-e_PRRC2C as it is the only cis-QTL to colocalize
# with this BMI signal.
locus3fdr_prrc2c=hyprcoloc(as.matrix(filt_betas[[3]]),as.matrix(filt_ses[[3]]),
                    trait.names=colnames(filt_betas[[3]]),snp.id=rownames(filt_betas[[3]]),ld.matrix = LD[[3]],
                    trait.subset = c(fdr3_traits,"cis-e_PRRC2C"),snpscores = T)
# Adding cis-e_PRRC2C to the FDR sig analysis does not appreciably diminish the colocalization (PP=0.9579 vs. 0.9701), which
# suggests that PRRC2C cis-eQTL could very well mediate the BMI association via these trans-aQTLs.

### 2p23.3, which is the 8th locus in the list.
# Start with all traits included
locus8all=hyprcoloc(as.matrix(filt_betas[[8]]),as.matrix(filt_ses[[8]]),
                    trait.names=colnames(filt_betas[[8]]),snp.id=rownames(filt_betas[[8]]),ld.matrix = LD[[8]],
                    snpscores = T)
# There were 3 clusters that each excluded BMI. The first was the largest with 2 cis-eQTLs, 16 trans-eQTLs and 52 trans-aQTLs. The 
# second cluster contained only cis-eQTLs, 3 trans-eQTLs and 1 trans-aQTL. The third cluster contained 2 cis-eQTLs, 1 cis-aQTL, 12
# trans-eQTLs and 7 trans-aQTLs. Only the third cluster had a top candidate variant that was a significant BMI GWAS variant. It was
# also the cluster that contained the strongest pairwaise colocalization with BMI, namely trans-e_PITX2. Cluster 1 contained the
# only other pairwise colocalized QTL, namely trans-e_BMP3.

# Let's run an analysis with only the QTLs that co-localized with BMI in the pairwise analysis.
coloc8_traits=c("BMI",gsub("BMI, ","",posColoc[posColoc$locus=="2p23.3",2]))
locus8coloc=hyprcoloc(as.matrix(filt_betas[[8]]),as.matrix(filt_ses[[8]]),
                      trait.names=colnames(filt_betas[[8]]),snp.id=rownames(filt_betas[[8]]),ld.matrix = LD[[8]],
                      trait.subset = coloc8_traits,snpscores = T)
# By restricting the multi-trait analysis to the pairwise colocalized QTLs, there is one cluster that includes BMI as well as both
# pairwise colocalized QTLs with a PP=0.5444 and rs12468863 as top candidate variant.

# Let's run an analysis with only the QTLs that were FDR significant in the original eQTL/aQTL analyses.
fdr8_traits=c("BMI","cis-e_MAPRE3","cis-e_SNX17","cis-e_IFT172","cis-e_SUPT7L","cis-e_ASXL2","cis-e_CAD","cis-e_PPP1CB","cis-a_SNX17",
              "trans-a_ADH1A")
locus8fdr=hyprcoloc(as.matrix(filt_betas[[8]]),as.matrix(filt_ses[[8]]),
                    trait.names=colnames(filt_betas[[8]]),snp.id=rownames(filt_betas[[8]]),ld.matrix = LD[[3]],
                    trait.subset = fdr8_traits,snpscores = T)
# There is no cluster identified. However, it is worth noting that 4 of the FDR significant cis genes are absent from the analysis due
# to the NA filtering above.

### 4q22_loc1, which is the 13th locus in the list.
# Start with all traits included
locus13all=hyprcoloc(as.matrix(filt_betas[[13]]),as.matrix(filt_ses[[13]]),
                    trait.names=colnames(filt_betas[[13]]),snp.id=rownames(filt_betas[[13]]),ld.matrix = LD[[13]],
                    snpscores = T)
# No clusters.

# Let's run an analysis with only the QTLs that co-localized with BMI in the pairwise analysis.
coloc13_traits=c("BMI",gsub("BMI, ","",posColoc[posColoc$locus=="4q22_loc1",2]))
locus13coloc=hyprcoloc(as.matrix(filt_betas[[13]]),as.matrix(filt_ses[[13]]),
                      trait.names=colnames(filt_betas[[13]]),snp.id=rownames(filt_betas[[13]]),ld.matrix = LD[[13]],
                      trait.subset = coloc13_traits,snpscores = T)
# When restricted to only the 2 QTLs (trans-e_PITX2 and trans-a_GNA14) that pairwise colocalized with BMI, they both cluster with BMI 
# with a PP=0.4033 and rs1903579 as the top candidate variant.

# The only FDR significant QTL for this BMI locus is cis-a_SCNA, but it did not pairwise colocalize with BMI nor cluster in the
# multi-trait analysis.

### 5q31.2, which is the 16th locus in the list.
# Start with all traits included
locus16all=hyprcoloc(as.matrix(filt_betas[[16]]),as.matrix(filt_ses[[16]]),
                    trait.names=colnames(filt_betas[[16]]),snp.id=rownames(filt_betas[[16]]),ld.matrix = LD[[16]],
                    snpscores = T)
locus16all_cluster1=strsplit(locus16all[[1]][1,2],", ")[[1]]
locus16all_cluster2=strsplit(locus16all[[1]][2,2],", ")[[1]]
# There were 3 clusters found. The first was the largest with 1 cis-eQTL, 29 trans-eQTLs and 67 trans-aQTLs, a PP=0.3612 that is entirely
# explained by rs217272 that is not even suggestively a BMI signal. The second cluster included BMI, 6 cis-eQTLs, 10 trans-eQTLs and 1
# trans-aQTL with a PP=0.3971 that is mostly explained by rs11750814 that is within 1 order of magnitude of the top BMI variant. The third
# cluster included 1 cis-aQTL, 8 trans-eQTLs and 2 trans-aQTLs with a PP=0.4685 that is mostly explained by rs529526 that is definitely not
# associated with BMI.

# Let's run an analysis with only the QTLs that co-localized with BMI in the pairwise analysis.
coloc16_traits=c("BMI",gsub("BMI, ","",posColoc[posColoc$locus=="5q31.2",2]))
locus16coloc=hyprcoloc(as.matrix(filt_betas[[16]]),as.matrix(filt_ses[[16]]),
                      trait.names=colnames(filt_betas[[16]]),snp.id=rownames(filt_betas[[16]]),ld.matrix = LD[[16]],
                      trait.subset = coloc16_traits,snpscores = T)
locus16coloc_cluster1=strsplit(locus16coloc[[1]][1,2],", ")[[1]]
# By restricting the multi-trait analysis to the pairwise colocalized QTLs, there is one cluster that includes BMI as well as all 28
# pairwise colocalized QTLs with a PP=0.4483 and rs11750814 explaining essentially all of the PP. Many of the SNPs in cluster 2 of the
# above all-QTL analysis are included in this cluster. This cluster has 3 cis-eQTLs, 10 trans-eQTLs and 15 trans-aQTLs.

# The only FDR significant QTL for this BMI locus is trans-a_ARRB1, which did successfully colocalize in the pairwise analysis with a 
# PP=0.7635 that was best explained by rs11750814 as well.

### 7q32, which is the 20th locus in the list.
# Start with all traits included
locus20all=hyprcoloc(as.matrix(filt_betas[[20]]),as.matrix(filt_ses[[20]]),
                    trait.names=colnames(filt_betas[[20]]),snp.id=rownames(filt_betas[[20]]),ld.matrix = LD[[20]],
                    snpscores = T)
# There were 2 clusters identified. The first included BMI, 3 cis-eQTLs, 1 cis-aQTL, 5 trans-eQTLs and 4 trans-aQTLs with a PP=0.4303 that
# was mostly explained by rs972283, which is the top BMI GWAS SNP. The second cluster is simply a pairwise colocalization between
# cis-e_MEST and cis-e_COPG2 with a PP=0.9975 mostly explained by rs1558917 that is not significantly associated with BMI.

# Let's run an analysis with only the QTLs that co-localized with BMI in the pairwise analysis.
coloc20_traits=c("BMI",gsub("BMI, ","",posColoc[posColoc$locus=="7q32",2]))
locus20coloc=hyprcoloc(as.matrix(filt_betas[[20]]),as.matrix(filt_ses[[20]]),
                       trait.names=colnames(filt_betas[[20]]),snp.id=rownames(filt_betas[[20]]),ld.matrix = LD[[20]],
                       trait.subset = coloc20_traits,snpscores = T)
# By restricting the multi-trait analysis to the pairwise colocalized QTLs, there is one cluster that includes BMI as well as all 7 
# pairwise colocalized QTLs (cis-e_LINC-PINT, cis-e_KLF14, cis-e_AC016831.7, cis-a_KLF14, trans-e_MKX, trans-e_TBX4 and trans-a_FGFRL1).
# All of these QTLs were also included in the cluster with BMI above in the all QTL analysis, but the PP is much improved to 0.83 that
# is mostly explained by rs972283 as well.

# Let's run an analysis with only the QTLs that were FDR significant in the original eQTL/aQTL analyses.
fdr20_traits=c("BMI","cis-e_LINC-PINT","cis-e_KLF14","cis-a_KLF14","trans-e_TBX4")
locus20fdr=hyprcoloc(as.matrix(filt_betas[[20]]),as.matrix(filt_ses[[20]]),
                    trait.names=colnames(filt_betas[[20]]),snp.id=rownames(filt_betas[[20]]),ld.matrix = LD[[20]],
                    trait.subset = fdr20_traits,snpscores = T)
# Further restricting the above analysis of pairwise colocalized QTLs to those that were FDR significant (which is a proper subset)
# resulted in a slight improvement in PP to 0.8941.

### 10q26.3, which is the 27th locus in the list.
# Start with all traits included
locus27all=hyprcoloc(as.matrix(filt_betas[[27]]),as.matrix(filt_ses[[27]]),
                     trait.names=colnames(filt_betas[[27]]),snp.id=rownames(filt_betas[[27]]),ld.matrix = LD[[27]],
                     snpscores = T)
locus27all_cluster1=strsplit(locus27all[[1]][1,2],", ")[[1]]
# There were 3 clusters identified. The first included BMI, 13 trans-eQTLs and 43 trans-aQTLs with a PP=0.4353 that was mostly explained 
# by rs871883, which is about 1 order of magnitude higer P than the top BMI GWAS SNP. The second cluster contains cis-e_DPYSL4 and 2 trans-aQTLs
# with a PP=0.4496 best explained by rs2814163. The third cluster includes 4 trans-eQTLs and 10 trans-aQTLs with a PP=0.4564 almost
# entirely explained by rs4897737.

# Let's run an analysis with only the QTLs that co-localized with BMI in the pairwise analysis.
coloc27_traits=c("BMI",gsub("BMI, ","",posColoc[posColoc$locus=="10q26.3",2]))
locus27coloc=hyprcoloc(as.matrix(filt_betas[[27]]),as.matrix(filt_ses[[27]]),
                       trait.names=colnames(filt_betas[[27]]),snp.id=rownames(filt_betas[[27]]),ld.matrix = LD[[27]],
                       trait.subset = coloc27_traits,snpscores = T)
# By restricting the multi-trait analysis to the pairwise colocalized QTLs, there is one cluster that includes BMI as well as all 15 
# pairwise colocalized QTLs (4 trans-eQTLs and 11 trans-aQTLs). All of these QTLs were also included in the cluster with BMI above in the 
# all QTL analysis, but the PP is much improved to 0.8606 that is mostly explained by rs902627 as well.

# Let's run an analysis with only the QTLs that were FDR significant in the original eQTL/aQTL analyses.
fdr27_traits=c("BMI","trans-a_ARHGEF37","trans-a_POFUT1","trans-a_TENM4")
locus27fdr=hyprcoloc(as.matrix(filt_betas[[27]]),as.matrix(filt_ses[[27]]),
                     trait.names=colnames(filt_betas[[27]]),snp.id=rownames(filt_betas[[27]]),ld.matrix = LD[[27]],
                     trait.subset = fdr27_traits,snpscores = T)
# Further restricting the above analysis of pairwise colocalized QTLs to those 3 trans-aQTLs that were FDR significant (which is a proper subset)
# resulted in a slightly diminished PP=0.817. The cluster is still best explained by rs902627.

### 12p13.1, which is the 31st locus in the list.
# Start with all traits included
locus31all=hyprcoloc(as.matrix(filt_betas[[31]]),as.matrix(filt_ses[[31]]),
                     trait.names=colnames(filt_betas[[31]]),snp.id=rownames(filt_betas[[31]]),ld.matrix = LD[[31]],
                     snpscores = T)
locus31all_cluster1=strsplit(locus31all[[1]][1,2],", ")[[1]]
# There was a single cluster including BMI, 23 trans-eQTLs and 74 trans-aQTLs with PP=0.4272 entirely explained by rs12422552, which is the
# only genome-wide significant BMI SNP for this locus.

# Let's run an analysis with only the QTLs that co-localized with BMI in the pairwise analysis.
coloc31_traits=c("BMI",gsub("BMI, ","",posColoc[posColoc$locus=="12p13.1",2]))
locus31coloc=hyprcoloc(as.matrix(filt_betas[[31]]),as.matrix(filt_ses[[31]]),
                       trait.names=colnames(filt_betas[[31]]),snp.id=rownames(filt_betas[[31]]),ld.matrix = LD[[31]],
                       trait.subset = coloc31_traits,snpscores = T)
locus31coloc_cluster1=strsplit(locus31coloc[[1]][1,2],", ")[[1]]
all(locus31coloc_cluster1 %in% locus31all_cluster1) # TRUE
# By restricting the multi-trait analysis to the pairwise colocalized QTLs, there is one cluster that includes BMI as well as all 81 
# pairwise colocalized QTLs (15 trans-eQTLs and 66 trans-aQTLs). All of these QTLs were also included in the cluster with BMI above in the 
# all QTL analysis, but the PP is much improved to 0.6096 that is entirely explained by rs12422552 as well.

# Let's run an analysis with only the QTLs that were FDR significant in the original eQTL/aQTL analyses.
fdr31_traits=c("BMI","trans-a_ANG","trans-a_CSNK2A2","trans-a_ID2","trans-a_PIM1","trans-a_PTPRJ","trans-a_TENM4","trans-a_TNFRSF10C","trans-a_ZFAT")
locus31fdr=hyprcoloc(as.matrix(filt_betas[[31]]),as.matrix(filt_ses[[31]]),
                     trait.names=colnames(filt_betas[[31]]),snp.id=rownames(filt_betas[[31]]),ld.matrix = LD[[31]],
                     trait.subset = fdr31_traits,snpscores = T)
# Further restricting the above analysis of pairwise colocalized QTLs to the 8 trans-aQTLs that were FDR significant (which is a proper subset)
# resulted in a very impressive PP=0.9787. Of course, rs12422552 is still fully explaining this PP.

# Write results to text files
write.table(ephb2_coloc,"Pairwise_HyPrColoc_between_EPHB2_aQTL_and_each_other_QTL_for_1p36.txt",sep = "\t",row.names = F,quote = F)
write.table(locus1all[[1]],"HyPrColoc_of_BMI_and_all_QTLs_for_1p36.txt",sep = "\t",row.names = F,quote = F)
write.table(bmi_ephb2_3way,"3-way_HyPrColoc_of_BMI_EPHB2_aQTL_and_each_other_QTL_for_1p36.txt",sep = "\t",row.names = F,quote = F)
write.table(locus1coloc1[[1]],"HyPrColoc_of_BMI_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_1p36.txt",sep = "\t",row.names = F,quote = F)
write.table(locus1fdr[[1]],"HyPrColoc_of_BMI_and_all_FDR_sig_QTLs_for_1p36.txt",sep = "\t",row.names = F,quote = F)
write.table(locus3all[[1]],"HyPrColoc_of_BMI_and_all_QTLs_for_1q24_loc1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus3coloc[[1]],"HyPrColoc_of_BMI_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_1q24_loc1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus3fdr[[1]],"HyPrColoc_of_BMI_and_all_FDR_sig_QTLs_for_1q24_loc1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus3fdr_prrc2c[[1]],"HyPrColoc_of_BMI_and_all_FDR_sig_QTLs_plus_PRCC2C_cis-eQTL_for_1q24_loc1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus8all[[1]],"HyPrColoc_of_BMI_and_all_QTLs_for_2p23.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus8coloc[[1]],"HyPrColoc_of_BMI_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_2p23.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus13all[[1]],"HyPrColoc_of_BMI_and_all_QTLs_for_4q22_loc1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus13coloc[[1]],"HyPrColoc_of_BMI_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_4q22_loc1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus16all[[1]],"HyPrColoc_of_BMI_and_all_QTLs_for_5q31.2.txt",sep = "\t",row.names = F,quote = F)
write.table(locus16coloc[[1]],"HyPrColoc_of_BMI_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_5q31.2.txt",sep = "\t",row.names = F,quote = F)
write.table(locus20all[[1]],"HyPrColoc_of_BMI_and_all_QTLs_for_7q32.txt",sep = "\t",row.names = F,quote = F)
write.table(locus20coloc[[1]],"HyPrColoc_of_BMI_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_7q32.txt",sep = "\t",row.names = F,quote = F)
write.table(locus20fdr[[1]],"HyPrColoc_of_BMI_and_all_FDR_sig_QTLs_for_7q32.txt",sep = "\t",row.names = F,quote = F)
write.table(locus27all[[1]],"HyPrColoc_of_BMI_and_all_QTLs_for_10q26.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus27coloc[[1]],"HyPrColoc_of_BMI_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_10q26.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus27fdr[[1]],"HyPrColoc_of_BMI_and_all_FDR_sig_QTLs_for_10q26.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus31all[[1]],"HyPrColoc_of_BMI_and_all_QTLs_for_12p13.1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus31coloc[[1]],"HyPrColoc_of_BMI_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_12p13.1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus31fdr[[1]],"HyPrColoc_of_BMI_and_all_FDR_sig_QTLs_for_12p13.1.txt",sep = "\t",row.names = F,quote = F)

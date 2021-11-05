### This script is for the co-localization analyses of the selected loci between WHRadjBMI GWAS and cis and trans eQTLs/aQTLs.

#install.packages("devtools")
#library(devtools)
#install_github("jrs95/hyprcoloc", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = F)
#browseVignettes("hyprcoloc") The install kept failing when trying to build the vignettes, so I disabled that.
#devtools::install_github("boxiangliu/locuscomparer")
library(hyprcoloc)
library(locuscomparer)

setwd("YOUR WORKING DIRECTORY")

# Read in the LD matrices, and the GWAS and QTL data. The file locations are relative to your working directory, so adjust accordingly.
gwas=read.table("GIANT_2015_WHRadjBMI_COMBINED_EUR_select_loci_with_locations.txt",sep = "\t",header = T)
cis_eqtl=read.table("./Eurobats_adipose_select_loci_cis-eQTLs_from_INT_logTPM.txt",sep = "\t",header = T)
cis_aqtl=read.table("./Eurobats_adipose_select_loci_cis-aQTLs_from_unnormalized_activities.txt",sep = "\t",header = T)
trans_eqtl=read.table("./Eurobats_adipose_select_loci_trans-eQTLs_for_WHR_MRs.txt",sep = "\t",header = T)
trans_aqtl=read.table("./Eurobats_adipose_select_loci_trans-aQTLs_for_WHR_MRs.txt",sep = "\t",header = T)
snp_pos=read.table("./Eurobats_QCd_select_loci_locations_for_colocalizations.txt",sep = "\t",header = T)
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

# Check that the SNPs included in each cis and trans QTL analysis are the same
all(unique(cis_eqtl$snps) %in% unique(cis_aqtl$snps)) # FALSE
all(unique(cis_eqtl$snps) %in% unique(trans_eqtl$snps)) # TRUE
all(unique(cis_aqtl$snps) %in% unique(trans_eqtl$snps)) # TRUE
all(unique(trans_eqtl$snps) %in% unique(trans_aqtl$snps)) # TRUE
# The cis-eQTL/aQTL SNPs do not match because their sets of tested genes are different and Matrix eQTL only tests cis-SNPs with 500kb. 
# Note also that both the cis-eQTL and cis-aQTL SNPs are proper subsets of the trans-QTL SNPs.
length(unique(cis_eqtl$snps)) # 78293
length(unique(cis_aqtl$snps)) # 74036
length(unique(trans_eqtl$snps)) # 79838
length(unique(trans_aqtl$snps)) # 79838

# Before filtering, let's standardize the gwas column names a bit
colnames(gwas)=c("SNP","CHR","POS","Allele1","Allele2","FreqAlleleHapMapCEU","BETA","SE","p","N")

# Let's get rid of "chr" before the chr numbers in gwas
gwas$CHR=as.numeric(gsub("chr","",gwas$CHR))

# Filter GWAS and QTL data to the same SNPs
filt_gwas=gwas[na.omit(match(unique(cis_eqtl$snps),gwas$SNP)),]
filt_cis_eqtl=cis_eqtl[cis_eqtl$snps %in% filt_gwas$SNP,]
filt_cis_aqtl=cis_aqtl[cis_aqtl$snps %in% filt_gwas$SNP,]
filt_trans_eqtl=trans_eqtl[trans_eqtl$snps %in% filt_gwas$SNP,]
filt_trans_aqtl=trans_aqtl[trans_aqtl$snps %in% filt_gwas$SNP,]
levels(filt_trans_aqtl$gene)[!(levels(filt_trans_aqtl$gene) %in% levels(filt_trans_eqtl$gene))]

# Calculate the SE of the Beta for the QTLs with the formula SE=beta/statistic
filt_cis_eqtl$SE=filt_cis_eqtl$beta/filt_cis_eqtl$statistic
filt_cis_aqtl$SE=filt_cis_aqtl$beta/filt_cis_aqtl$statistic
filt_trans_eqtl$SE=filt_trans_eqtl$beta/filt_trans_eqtl$statistic
filt_trans_aqtl$SE=filt_trans_aqtl$beta/filt_trans_aqtl$statistic

# Add chromosome and position to the QTLs for sorting
filt_cis_eqtl$chr=filt_gwas[match(filt_cis_eqtl$snps,filt_gwas$SNP),"CHR"]
filt_cis_aqtl$chr=filt_gwas[match(filt_cis_aqtl$snps,filt_gwas$SNP),"CHR"]
filt_trans_eqtl$chr=filt_gwas[match(filt_trans_eqtl$snps,filt_gwas$SNP),"CHR"]
filt_trans_aqtl$chr=filt_gwas[match(filt_trans_aqtl$snps,filt_gwas$SNP),"CHR"]
filt_cis_eqtl$position=filt_gwas[match(filt_cis_eqtl$snps,filt_gwas$SNP),"POS"]
filt_cis_aqtl$position=filt_gwas[match(filt_cis_aqtl$snps,filt_gwas$SNP),"POS"]
filt_trans_eqtl$position=filt_gwas[match(filt_trans_eqtl$snps,filt_gwas$SNP),"POS"]
filt_trans_aqtl$position=filt_gwas[match(filt_trans_aqtl$snps,filt_gwas$SNP),"POS"]

# Sort by chr and position
filt_gwas=filt_gwas[order(filt_gwas$CHR,filt_gwas$POS),]
filt_cis_eqtl=filt_cis_eqtl[order(filt_cis_eqtl$chr,filt_cis_eqtl$position),]
filt_cis_aqtl=filt_cis_aqtl[order(filt_cis_aqtl$chr,filt_cis_aqtl$position),]
filt_trans_eqtl=filt_trans_eqtl[order(filt_trans_eqtl$chr,filt_trans_eqtl$position),]
filt_trans_aqtl=filt_trans_aqtl[order(filt_trans_aqtl$chr,filt_trans_aqtl$position),]

# Split GWAS and QTL Betas and SEs by locus, and also by genes for the QTL data. This section is pretty complex and hard to follow,
# for which I am both proud and embarrassed.

# After running into some errors down below and then investigating the source, I discovered that the 1q24_loc1 and 1q24_loc2 loci 
# are not correctly split by gaps between SNPs because they are too close together. Therefore, I need to rework the locus splitting
# process based on the LD matrices.

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
  colnames(loci_betas[[i]])="WHR"
  
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
  colnames(loci_ses[[i]])="WHR"
  
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

# Filter out SNPs with any remaining NAs
for(i in 1:43){
  snp_nas=c()
  for(j in rownames(filt_betas[[i]])){
    snp_nas[j]=sum(is.na(filt_betas[[i]][j,]))
  }
  filt_betas[[i]]=filt_betas[[i]][snp_nas==0,]
  filt_ses[[i]]=filt_ses[[i]][snp_nas==0,]
}

# Filter the LD matrices to include only the SNPs remaining
LD=list()
for(i in 1:length(ld)){
  LD[[i]]=ld[[i]][rownames(filt_betas[[i]]),rownames(filt_betas[[i]])]
}

# First run HyPrColoc pairwise between GWAS trait and QTLs. Then groups of traits may be analyzed together in a more targeted fashion.

# Pairwise co-localization, which is essentially equivalent with using coloc.
# Make a data.frame with only the QTLs co-localized with WHRadjBMI
GWAScoloc=list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43) # initialize a list of 43 objects
posColoc=data.frame() # initialize a data.frame for the positive co-localized traits
for(i in 1:43){
  # Grab the traits list for each locus
  traits=colnames(filt_betas[[i]])
  
  # Iterate through each of the QTLs. Exclude BMI from the index of traits since all will be tested against BMI.
  for(j in traits[-1]){
    GWAScoloc[[i]][[j]]=hyprcoloc(as.matrix(filt_betas[[i]]),as.matrix(filt_ses[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas[[i]]),ld.matrix = LD[[i]],
                                  trait.subset = c("WHR",j))
    # Grab only the QTLs with co-localization into posColoc
    if(!is.na(GWAScoloc[[i]][[j]]$results[1,3])){
      posColoc=rbind(posColoc,cbind(GWAScoloc[[i]][[j]]$results,"locus"=i))
    }
  }
}

# Let's change the locus number to legit chromosomal locus ID
loci=read.table("../Select_loci_list.txt",sep = "\t",header = T)
posColoc$locus=loci[match(posColoc$locus,rownames(loci)),1]

# Let's tally the co-localized QTLs per locus
pair_tallies=data.frame(row.names = loci$Locus)
for(i in loci$Locus){
  pair_tallies[i,1]=loci[loci$Locus==i,2]
  pair_tallies[i,2]=dim(posColoc[posColoc$locus==i,])[1]
}
colnames(pair_tallies)=c("GWAS","Colocalizations with WHR")
# HyPrColoc reports colocalizations with PP>0.25, which is what is tallied here.
# Locus                                   GWAS                        Colocalizations with WHR
# 1p36.1                                   BMI                        0
# 1p13.3                                   HDL                        0
# 1q24_loc1                                BMI                        0
# 1q24_loc2                                WHR                        9
# 1q25.2                                   BMI                        0
# 1q41                                     WHR                        0
# 2p24                           Triglycerides                        0
# 2p23.3                                   BMI                        0
# 2p21                                     T2D                        0
# 2q24.3                                   WHR                        0
# 3p25                                     HDL                        0
# 3q27                                     T2D                        0
# 4q22_loc1                                BMI                        0
# 4q22_loc2                                BMI                        0
# 5q11.2                           HDL and BMI                        0
# 5q31.2                                   BMI                        0
# 6p21_loc1                      Triglycerides                        1
# 6p21_loc2                                HDL                        0
# 6q22.3                                   WHR                        2
# 7q32         BMI, T2D, HDL and Triglycerides                        0
# 7q36                                     HDL                        0
# 8q21.1                                   BMI                        0
# 8q21.2                                   BMI                        0
# 10p13                                    BMI                        0
# 10q22.2                                  BMI                        0
# 10q24.2                                  BMI                        0
# 10q26.3                                  BMI                        0
# 11p11.2                                  HDL                        0
# 11q13                                    BMI                        1
# 12p13.33                                 BMI                        0
# 12p13.1                                  BMI                        0
# 12q13.13                                 BMI                        1
# 15q21.3                                  HDL                        0
# 15q24.1                                  BMI                        0
# 17q21.2                                  BMI                        0
# 17q25.3                                  HDL                        0
# 18q21.1                                  HDL                        0
# 18q21.3                                  BMI                        0
# 19p13.1                T2D and Triglycerides                        0
# 19q13.3_loc1                             HDL                        0
# 19q13.3_loc2                             BMI                        0
# 20q13.32                                 T2D                        0
# 22q13.3                                  T2D                        0

# For WHRadjBMI, 1q24_loc2, 1q41, 2q24.3 and 6q22.3 were of interest, especially 1q24_loc2 due to the presence of 44 FDR significant trans-aGenes (1 trans-eGenes).
# In the pairwise analysis above, despite 44 significant WHR MR tran-aGenes for 1q24_loc2, only 9 trans-QTLs show co-localization
# and 2 of those are trans-eQTLs (DAPK2 and HOXA3). Let's see if this changes at all in the multi-trait co-localization below. 

# Let's write the pairwise analysis to file for reference in the rest of these analyses.
write.table(posColoc,"Pairwise_HyPrColoc_between_BMIadjWHR_and_each_QTL_for_select_loci.txt",sep = "\t",row.names = F,quote = F)

### Mult-trait analyses

# 1q24_loc2 (which is the 4th locus in the list)
# Start with all traits included
locus4all=hyprcoloc(as.matrix(filt_betas[[4]]),as.matrix(filt_ses[[4]]),
                    trait.names=colnames(filt_betas[[4]]),snp.id=rownames(filt_betas[[4]]),ld.matrix = LD[[4]],
                    snpscores = T)
locus4all_cluster1=strsplit(locus4all[[1]][1,2],", ")[[1]]
# In the pairwise analysis above, those 9 trans-QTLs that co-localized all had rs714515 as the top candidate SNP. This seems to be
# the top WHRadjBMI GWAS SNP in the region, but only has FDR significant trans-QTLs with HOXA3 and PMEPA1. However, in this
# multi-trait analysis, there is a cluster that includes WHRadjBMI, 20 trans-eQTLs and 68 trans-aQTLs. The PP (0.3968) for this cluster is
# not terribly impressive, but the top candidate SNP rs9425301 explains essentially all of the PP and has FDR significant trans-aQTLs
# with 41 of the WHR MRs and a WHRadjBMI GWAS P=3.8E-13. The r2 between rs714515 and rs9425301 is 0.697. One other cluster contained
# 8 trans-eQTLs, and the top candidate SNP (rs17368942) explains about half of the PP and has a GWAS P=9.3E-9 (5 orders worse than top SNP).

# Let's run an analysis with only the QTLs that co-localized with WHRadjBMI in the pairwise analysis.
coloc4_traits=c("WHR",gsub("WHR, ","",posColoc[posColoc$locus=="1q24_loc2",2]))
locus4coloc=hyprcoloc(as.matrix(filt_betas[[4]]),as.matrix(filt_ses[[4]]),
                          trait.names=colnames(filt_betas[[4]]),snp.id=rownames(filt_betas[[4]]),ld.matrix = LD[[4]],
                          trait.subset = coloc4_traits,snpscores = T)
locus4coloc_cluster1=strsplit(locus4coloc[[1]][1,2],", ")[[1]]
# By restricting the multi-trait co-localization to those QTLs that co-localized in the pairwise analysis, we get one
# cluster that includes 9/9 trans-QTLs that co-localized before. This cluster has a PP=0.7896 and the top
# candidate SNP is rs9425301, which is consistent with the all QTL multi-trait analysis.

# Let's see what happens when run HyPrColoc on WHRadjBMI with all of the FDR significant trans-aQTLs.
locus4fdr=hyprcoloc(as.matrix(filt_betas[[4]]),as.matrix(filt_ses[[4]]),
                    trait.names=colnames(filt_betas[[4]]),snp.id=rownames(filt_betas[[4]]),ld.matrix = LD[[4]],
                    trait.subset = c("WHR","trans-e_TNFAIP6","trans-a_ADRB1","trans-a_ANG","trans-a_ANXA5","trans-a_ARHGEF26","trans-a_BMP3",
                                     "trans-a_CCND2","trans-a_CECR2","trans-a_CIDEA","trans-a_CPE","trans-a_CSNK2A2","trans-a_DAB2IP",
                                     "trans-a_DAPK2","trans-a_DOK5","trans-a_EDARADD","trans-a_GIPR","trans-a_GLIS1","trans-a_GLRA1",
                                     "trans-a_GNA14","trans-a_GPR39","trans-a_HEXIM1","trans-a_HMOX1","trans-a_HOXA3","trans-a_ID2",
                                     "trans-a_KCNIP2","trans-a_MKX","trans-a_MSC","trans-a_MSN","trans-a_P2RX6","trans-a_PMEPA1",
                                     "trans-a_PTPRJ","trans-a_RASL10B","trans-a_RASSF4","trans-a_RGS3","trans-a_RND3","trans-a_RTN1",
                                     "trans-a_TBX4","trans-a_TNFAIP6","trans-a_TYRO3","trans-a_WWC1","trans-a_YWHAH","trans-a_ZNF334"),
                    snpscores = T)
locus4fdr_cluster1=strsplit(locus4fdr[[1]][1,2],", ")[[1]]
sum(locus4fdr_cluster1 %in% coloc4_traits)-1 # All 7 pairwise trans-aQTLs overlap between the pairwise and the fdr cluster1
sum(locus4fdr_cluster1 %in% locus4all_cluster1)-1 # All 45 of the fdr HyPrColoc cluster1 QTLs are also in the all QTL cluster1
# In this case the one FDR sig trans-eQTL and 41/41 trans-aQTLs for this locus co-localize in a single cluster with WHRadjBMI with a PP=0.7436
# which is fully explained by rs9425301, as for the all QTL HyPrColoc analysis above.

# 6q22.3 (which is the 4th locus in the list)
# Start with all traits included
locus19all=hyprcoloc(as.matrix(filt_betas[[19]]),as.matrix(filt_ses[[19]]),
                    trait.names=colnames(filt_betas[[19]]),snp.id=rownames(filt_betas[[19]]),ld.matrix = LD[[19]],
                    snpscores = T)
# This BMIadjWHR locus was of interest due to its having 1 FDR sig trans-eQTL with DMRTA1 and 1 trans-aQTL with TACSTD2. However, the only
# pairwise colocalized QTLs was cis-e_RSPO3 and cis-e_RSPO1, which is interesting given they are clearly in the same gene family. However,
# in the all QTL multi-trait analysis, only trans-e_RSPO1 clusters with BMIadjWHR with PP=0.8526 for which 0.2802 is explained by
# rs2745353 that is in perfect LD with the top GWAS SNP. There is another cluster that includes 2 cis-eQTLs, 6 trans-eQTLs and 5 trans-aQTLs 
# with PP=0.4524 that is mostly explained by rs17054787 that is not a significant BMIadjWHR GWAS variant.

# Let's run an analysis with only the QTLs that co-localized with WHRadjBMI in the pairwise analysis.
coloc19_traits=c("WHR",gsub("WHR, ","",posColoc[posColoc$locus=="6q22.3",2]))
locus19coloc=hyprcoloc(as.matrix(filt_betas[[19]]),as.matrix(filt_ses[[19]]),
                      trait.names=colnames(filt_betas[[19]]),snp.id=rownames(filt_betas[[19]]),ld.matrix = LD[[19]],
                      trait.subset = coloc19_traits,snpscores = T)
# Both pairwise QTLs cluster with BMIadjWHR with PP=0.8642 that is again best explained by rs2745353.

# Let's see what happens when run HyPrColoc on WHRadjBMI with all of the FDR significant trans-aQTLs.
locus19fdr=hyprcoloc(as.matrix(filt_betas[[19]]),as.matrix(filt_ses[[19]]),
                    trait.names=colnames(filt_betas[[19]]),snp.id=rownames(filt_betas[[19]]),ld.matrix = LD[[19]],
                    trait.subset = c("WHR","trans-e_DMRTA1","trans-a_TACSTD2"),snpscores = T)
# Neither of the FDR significant trans-QTLs cluster with BMIadjWHR. which is consistent with both the all QTL multi-trait analysis
# and the pairwise analysis. I checked whether RSPO3 was included in the ARACNe and VIPER analysis, but unfortunately it was not. If
# it has similar functions as RSPO1, which was considered not only a regulator but a WHR MR, then why wouldn't it also be counted as
# an expression regulator? I looked them up in NCBI Gene, and they both code for ligands that bind leucine-rich repeat-containing
# G-protein coupled receptors and regulate Wnt signaling. Interestingly, serum RSPO1 was recently identified as a new surrogate
# marker for obesity and insulin resistance, which would be consistent with the notion of RSPO1 as a WHR MR. I checked the Pearson
# correlation between RSPO1 and RSPO3 expression and it was only 0.20, so not real impressive. Consequently, I don't think RSPO3
# and RSPO1 could be interchangeable as WHR MRs. 

# Write results to text files
write.table(posColoc,"Pairwise_HyPrColoc_between_WHRadjBMI_and_each_QTL_for_select_loci.txt",sep = "\t",row.names = F,quote = F)
write.table(locus4all[[1]],"HyPrColoc_of_WHRadjBMI_and_all_QTLs_for_1q24_loc2.txt",sep = "\t",row.names = F,quote = F)
write.table(locus4coloc[[1]],"HyPrColoc_of_WHRadjBMI_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_1q24_loc2.txt",sep = "\t",row.names = F,quote = F)
write.table(locus4fdr[[1]],"HyPrColoc_of_WHRadjBMI_and_all_FDR_sig_trans-aQTLs_for_1q24_loc2.txt",sep = "\t",row.names = F,quote = F)
write.table(locus19all[[1]],"HyPrColoc_of_WHRadjBMI_and_all_QTLs_for_6q22.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus19coloc[[1]],"HyPrColoc_of_WHRadjBMI_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_6q22.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus19fdr[[1]],"HyPrColoc_of_WHRadjBMI_and_all_FDR_sig_trans-aQTLs_for_6q22.3.txt",sep = "\t",row.names = F,quote = F)


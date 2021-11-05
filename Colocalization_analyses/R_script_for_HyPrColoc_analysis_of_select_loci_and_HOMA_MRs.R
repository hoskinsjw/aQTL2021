### This script is for the co-localization analyses of the selected loci between T2D GWAS and cis and trans eQTLs/aQTLs.

#install.packages("devtools")
#library(devtools)
#install_github("jrs95/hyprcoloc", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = F)
#browseVignettes("hyprcoloc") The install kept failing when trying to build the vignettes, so I disabled that.
#devtools::install_github("boxiangliu/locuscomparer")
library(hyprcoloc)
library(locuscomparer)

setwd("YOUR WORKING DIRECTORY")

# Read in the LD matrices, and the GWAS and QTL data. The file locations are relative to your working directory, so adjust accordingly.
gwas=read.table("Mahajan.NatGenet2018b.T2Dbmiadj.European.with.rsIDs.txt",sep = "\t",header = T)
cis_eqtl=read.table("./Eurobats_adipose_select_loci_cis-eQTLs_from_INT_logTPM.txt",sep = "\t",header = T)
cis_aqtl=read.table("./Eurobats_adipose_select_loci_cis-aQTLs_from_unnormalized_activities.txt",sep = "\t",header = T)
trans_eqtl=read.table("./Eurobats_adipose_select_loci_trans-eQTLs_for_HOMA-IR_MRs.txt",sep = "\t",header = T)
trans_aqtl=read.table("./Eurobats_adipose_select_loci_trans-aQTLs_for_HOMA-IR_MRs.txt",sep = "\t",header = T)
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
colnames(gwas)=c("SNP","Chr","Pos","EA","NEA","EA_freq","Beta","SE","PValue","Neff")

# Filter GWAS and QTL data to the same SNPs
filt_gwas=gwas[na.omit(match(unique(cis_aqtl$snps),gwas$SNP)),] # match to the cis-aQTL which has the fewest unique SNPs
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
filt_cis_eqtl$chr=filt_gwas[match(filt_cis_eqtl$snps,filt_gwas$rsID),"Chr"]
filt_cis_aqtl$chr=filt_gwas[match(filt_cis_aqtl$snps,filt_gwas$rsID),"Chr"]
filt_trans_eqtl$chr=filt_gwas[match(filt_trans_eqtl$snps,filt_gwas$rsID),"Chr"]
filt_trans_aqtl$chr=filt_gwas[match(filt_trans_aqtl$snps,filt_gwas$rsID),"Chr"]
filt_cis_eqtl$position=filt_gwas[match(filt_cis_eqtl$snps,filt_gwas$rsID),"Pos"]
filt_cis_aqtl$position=filt_gwas[match(filt_cis_aqtl$snps,filt_gwas$rsID),"Pos"]
filt_trans_eqtl$position=filt_gwas[match(filt_trans_eqtl$snps,filt_gwas$rsID),"Pos"]
filt_trans_aqtl$position=filt_gwas[match(filt_trans_aqtl$snps,filt_gwas$rsID),"Pos"]

# Sort by chr and position
filt_gwas=filt_gwas[order(filt_gwas$Chr,filt_gwas$Pos),]
filt_cis_eqtl=filt_cis_eqtl[order(filt_cis_eqtl$chr,filt_cis_eqtl$position),]
filt_cis_aqtl=filt_cis_aqtl[order(filt_cis_aqtl$chr,filt_cis_aqtl$position),]
filt_trans_eqtl=filt_trans_eqtl[order(filt_trans_eqtl$chr,filt_trans_eqtl$position),]
filt_trans_aqtl=filt_trans_aqtl[order(filt_trans_aqtl$chr,filt_trans_aqtl$position),]

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
  loci_betas[[i]]=as.data.frame(filt_gwas[filt_gwas$SNP %in% rownames(ld[[i]]),"Beta"],
                                row.names = as.character(filt_gwas[filt_gwas$SNP %in% rownames(ld[[i]]),"SNP"]))
  colnames(loci_betas[[i]])="T2D"
  
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
  colnames(loci_ses[[i]])="T2D"
  
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

# Filter the LD matrices to include only the SNPs remaining and put them in a list
LD=list()
for(i in 1:length(ld)){
  LD[[i]]=ld[[i]][rownames(filt_betas[[i]]),rownames(filt_betas[[i]])]
}

# First run HyPrColoc pairwise between GWAS trait and QTLs. Then groups of traits may be analyzed together in a more targeted fashion.

# Pairwise co-localization, which is essentially equivalent with using coloc.
# Make a data.frame with only the QTLs co-localized with T2D
GWAScoloc=list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43) # initialize a list of 43 objects
posColoc=data.frame() # initialize a data.frame for the positive co-localized traits
for(i in 1:43){
  # Grab the traits list for each locus
  traits=colnames(filt_betas[[i]])
  
  # Iterate through each of the QTLs. Exclude BMI from the index of traits since all will be tested against BMI.
  for(j in traits[-1]){
    GWAScoloc[[i]][[j]]=hyprcoloc(as.matrix(filt_betas[[i]]),as.matrix(filt_ses[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas[[i]]),ld.matrix = LD[[i]],
                                  trait.subset = c("T2D",j))
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
colnames(pair_tallies)=c("GWAS","Colocalizations with T2D")
# HyPrColoc reports colocalizations with PP>0.25, which is what is tallied here.
#Locus                                   GWAS            Colocalizations with T2D
#1p36.1                                   BMI                        0
#1p13.3                                   HDL                        0
#1q24_loc1                                BMI                        0
#1q24_loc2                                WHR                        0
#1q25.2                                   BMI                        0
#1q41                                     WHR                        0
#2p24                           Triglycerides                        0
#2p23.3                                   BMI                        0
#2p21                                     T2D                        0
#2q24.3                                   WHR                        0
#3p25                                     HDL                        0
#3q27                                     T2D                        8
#4q22_loc1                                BMI                        0
#4q22_loc2                                BMI                        0
#5q11.2                           HDL and T2D                        0
#5q31.2                                   BMI                        0
#6p21_loc1                      Triglycerides                        0
#6p21_loc2                                HDL                        0
#6q22.3                                   WHR                        0
#7q32         BMI, T2D, HDL and Triglycerides                       12
#7q36                                     HDL                        0
#8q21.1                                   BMI                        0
#8q21.2                                   BMI                        0
#10p13                                    BMI                        0
#10q22.2                                  BMI                        0
#10q24.2                                  BMI                        0
#10q26.3                                  BMI                        0
#11p11.2                                  HDL                        0
#11q13                                    BMI                        1
#12p13.33                                 BMI                        0
#12p13.1                                  BMI                        0
#12q13.13                                 BMI                        0
#15q21.3                                  HDL                        0
#15q24.1                                  BMI                        0
#17q21.2                                  BMI                        0
#17q25.3                                  HDL                        0
#18q21.1                                  HDL                        0
#18q21.3                                  BMI                        0
#19p13.1                T2D and Triglycerides                        0
#19q13.3_loc1                             HDL                        0
#19q13.3_loc2                             BMI                        1
#20q13.32                                 T2D                        2
#22q13.3                                  T2D                        0

# Let's write the pairwise analysis to file for reference in the rest of these analyses.
write.table(posColoc,"Pairwise_HyPrColoc_between_BMIadjT2D_and_each_QTL_for_select_loci.txt",sep = "\t",row.names = F,quote = F)

### Mult-trait analyses

# 2p21, which is the 9th locus in the list
# Start with all cis and HOMA-IR traits included
locus9all=hyprcoloc(as.matrix(filt_betas[[9]]),as.matrix(filt_ses[[9]]),
                       trait.names=colnames(filt_betas[[9]]),snp.id=rownames(filt_betas[[9]]),ld.matrix = LD[[9]],
                       snpscores = T)
locus9all_cluster1=strsplit(locus9all[[1]][1,2],", ")[[1]]
# There is only one cluster that does not include T2D, but does include cis-a_LRPPRC, 40 trans-eQTLs and 85 trans-aQTLs with a PP=0.43 
# that is entirely explained by rs71420034 that has a BMIadjT2D GWAS P=0.066, so clearly not relevant to that GWAS signal. To be sure there
# is nothing relevant to this cluster, I also checked the T2D GWAS where the P=0.23, so it is also irrelevant to BMI-associated T2D. This
# locus was selected based on 2 suggestive cis-aQTLs (MTA3 and ZFP36L2) that were better than their corresponding eQTLs.
# They didn't pairwise colocalize with BMIadjT2D, but let's see if they cluster with the GWAS in the multi-trait analysis.

# Let's run an analysis with only the QTLs that were suggestive in the original eQTL/aQTL analyses.
fdr9_traits=c("T2D","cis-a_MTA3","cis-a_ZFP36L2")
locus9fdr=hyprcoloc(as.matrix(filt_betas[[9]]),as.matrix(filt_ses[[9]]),
                    trait.names=colnames(filt_betas[[9]]),snp.id=rownames(filt_betas[[9]]),ld.matrix = LD[[9]],
                    trait.subset = fdr9_traits,snpscores = T)
# No cluster. Note that MTA3 had been filtered out due to too many NAs (26%), suggesting it is rather far from the the BMIadjT2D
# signal. However, given that the ZFP36L2 aQTL signal looks similar to that of MTA3 (based the significant SNPs) and the fact that
# the BMIadjT2D GWAS P for those aQTL significant SNPs was 19 orders of magnitude worse than the top GWAS SNP, I think its safe to
# drop this locus from consideration.

# 3q27, which is the 12th locus in the list
# Start with all cis and HOMA-IR traits included
locus12all=hyprcoloc(as.matrix(filt_betas[[12]]),as.matrix(filt_ses[[12]]),
                    trait.names=colnames(filt_betas[[12]]),snp.id=rownames(filt_betas[[12]]),ld.matrix = LD[[12]],
                    snpscores = T)
locus12all_cluster1=strsplit(locus12all[[1]][1,2],", ")[[1]]
# There is 1 cluster that excludes BMIadjT2D, but includes 2 cis-eQTLs (C3orf70 and MAP3K13), 1 cis-aQTL (ETV5), 14 HOMA-IR MR trans-eQTLs 
# and 58 HOMA-IR MR trans-aQTLs with PP=0.4038 that is entirely explained by rs71320320, which is the top aQTL for ETV5 and has a 
# BMIadjT2D P=1.9E-28.  All 7 HOMA-IR MR trans-QTLs that co-localized with BMIadjT2D in the pairwise analysis fall in cluster 1 even though 
# for 5 of them their pairwise top candidate SNP was rs34782298, which has a BMIadjT2D P=1.6E-50. Based on LDLink, rs34782298, and the 
# other two pairwise top candidate SNPs (rs4376068 and rs9858406) are very highly correlated (r2>0.95), but rs71320320 only has an r2 of 
# ~0.49 with them, consistent with it possibly representing an independent signal. This will have to be followed up by other means.

# Let's run an analysis with only the QTLs that co-localized with BMIadjT2D in the pairwise analysis.
coloc12_traits=c("T2D",gsub("T2D, ","",posColoc[posColoc$locus=="3q27",2]))
locus12coloc=hyprcoloc(as.matrix(filt_betas[[12]]),as.matrix(filt_ses[[12]]),
                          trait.names=colnames(filt_betas[[12]]),snp.id=rownames(filt_betas[[12]]),ld.matrix = LD[[12]],
                          trait.subset = coloc12_traits,snpscores = T)
locus12coloc_cluster1=strsplit(locus12coloc[[1]][1,2],", ")[[1]]
# By restricting the multi-trait analysis to the pairwise co-localizing QTLs all 8 pairwise QTLs cluster together with a PP=0.609 best 
# explained by rs34782298, which is consistent with the pairwise analysis.

# Let's run an analysis with only the QTLs that were FDR significant in the original eQTL/aQTL analyses.
fdr12_traits=c("T2D","cis-e_ETV5","cis-e_MAP3K13","cis-e_SENP2","cis-a_ETV5")
locus12fdr=hyprcoloc(as.matrix(filt_betas[[12]]),as.matrix(filt_ses[[12]]),
                    trait.names=colnames(filt_betas[[12]]),snp.id=rownames(filt_betas[[12]]),ld.matrix = LD[[12]],
                    trait.subset = fdr12_traits,snpscores = T)
# No clustering.

# 7q32, which is the 20th locus in the list
# Start with all cis and HOMA-IR traits included
locus20all=hyprcoloc(as.matrix(filt_betas[[20]]),as.matrix(filt_ses[[20]]),
                       trait.names=colnames(filt_betas[[20]]),snp.id=rownames(filt_betas[[20]]),ld.matrix = LD[[20]],
                       snpscores = T)
locus20all_cluster1=strsplit(locus20all[[1]][1,2],", ")[[1]]
coloc20_traits=c("T2D",gsub("T2D, ","",posColoc[posColoc$locus=="7q32",2]))
coloc20_traits[(coloc20_traits %in% locus20all_cluster1)] # 12
# BMIadjT2D did not cluster with any QTLs in this multi-trait analysis, including those that colocalized in the pairwise analysis. 
# There were 2 clusters (1 and 69) of QTLs. Cluster 1 included cis-e_KLF14, cis-e_AC016831.7, cis-e_LINC-PINT, cis-a_KLF14 and 16
# trans-eQTLs and 9 trans-aQTLs with a PP=0.3215 best explained by rs6467315 that has a BMIadjT2D GWAS P=6.1E-14, which is 2 orders of
# magnitude worse than the GWAS tag SNP rs61462211. Cluster 69 only has 2 cis-eQTLs with a PP=0.9951 best explained by rs1558917
# that has a BMIadjT2D GWAS P=4.4E-4.

# Let's run an analysis with only the QTLs that co-localized with BMIadjT2D in the pairwise analysis.
locus20coloc=hyprcoloc(as.matrix(filt_betas[[20]]),as.matrix(filt_ses[[20]]),
                      trait.names=colnames(filt_betas[[20]]),snp.id=rownames(filt_betas[[20]]),ld.matrix = LD[[20]],
                      trait.subset = coloc20_traits,snpscores = T)
locus20coloc_cluster1=strsplit(locus20coloc[[1]][1,2],", ")[[1]]
# All 12 QTLs that colocalized with BMIadjT2D in the pairwise analysis cluster together in this analysis with BMIadjT2D.
# The cluster has PP=0.5449 and is best explained by rs6467317 that has a BMIadjT2D GWAS P=1.9E-15, which is only 1 order of
# magnitude worse than the top GWAS SNP.

# Let's run an analysis with only the QTLs that were FDR significant in the original eQTL/aQTL analyses.
fdr20_traits=c("T2D","cis-e_KLF14","cis-e_LINC-PINT","cis-a_KLF14","trans-e_GNB1","trans-e_TBX4")
locus20fdr=hyprcoloc(as.matrix(filt_betas[[20]]),as.matrix(filt_ses[[20]]),
                     trait.names=colnames(filt_betas[[20]]),snp.id=rownames(filt_betas[[20]]),ld.matrix = LD[[20]],
                     trait.subset = fdr20_traits,snpscores = T)
# BMIadjT2D signal clusters with all 5 FDR sig QTLs with a PP=0.5178 that is best explained by rs35722851 that has a GWAS P=4.1E-15
# that is 1 order of magnitude worse than the best GWAS SNP.

# 19p13.1, which is the 39th locus in the list
# Start with all cis and HOMA-IR traits included
locus39all=hyprcoloc(as.matrix(filt_betas[[39]]),as.matrix(filt_ses[[39]]),
                     trait.names=colnames(filt_betas[[39]]),snp.id=rownames(filt_betas[[39]]),ld.matrix = LD[[39]],
                     snpscores = T)
# This yielded 2 cluster (1 and 11), but neither included BMIadjT2D. This locus was selected because it had a cis-aQTL for 
# ZNF101 that was stronger than the eQTL. However, it doesn't show up in either the pairwise or the all trait analyses, so I guess
# there's nothing more to do with this.

# 20q13.32, which is the 42nd locus in the list
# Start with all cis and HOMA-IR traits included
locus42all=hyprcoloc(as.matrix(filt_betas[[42]]),as.matrix(filt_ses[[42]]),
                     trait.names=colnames(filt_betas[[42]]),snp.id=rownames(filt_betas[[42]]),ld.matrix = LD[[42]],
                     snpscores = T)
# There were 2 clusters found (49 and 104). The first did not include BMIadjT2D, but did have 2 trans-eQTLs and 9 trans-aQTLs with a
# PP=0.4161 that was best explained by rs7268267 that has a GWAS P=0.021. The second cluster was simply a pairwise colocalization between
# BMIadjT2D and trans-e_RASSF4, which has the same statistics as in the pairwise analysis.

# Let's run an analysis with only the QTLs that co-localized with BMIadjT2D in the pairwise analysis.
coloc42_traits=c("T2D",gsub("T2D, ","",posColoc[posColoc$locus=="20q13.32",2]))
locus42coloc=hyprcoloc(as.matrix(filt_betas[[42]]),as.matrix(filt_ses[[42]]),
                       trait.names=colnames(filt_betas[[42]]),snp.id=rownames(filt_betas[[42]]),ld.matrix = LD[[42]],
                       trait.subset = coloc42_traits,snpscores = T)
# Both pairwise colocalized trans-eQTLs cluster with BMIadjT2D with a PP=0.5143 that is best explained by rs736266, which is the top
# GWAS SNP. Note that this is dramatically worse than the pairwise colocalization between BMIadjT2D and trans-e_RASSF4 (PP=0.9489), but
# is a bit better than the pairwise colocalization between BMIadjT2D and trans-e_ARHGEF37 (PP=0.4531).

# 22q13.3, which is the 43rd locus in the list
# Start with all cis and HOMA-IR traits included
locus43all=hyprcoloc(as.matrix(filt_betas[[43]]),as.matrix(filt_ses[[43]]),
                     trait.names=colnames(filt_betas[[43]]),snp.id=rownames(filt_betas[[43]]),ld.matrix = LD[[43]],
                     snpscores = T)
# This yielded 4 cluster (1, 2, 3 and 103), but none included BMIadjT2D. This locus was selected because it had a cis-aQTL for 
# RABL2B that was stronger than the eQTL. However, it doesn't show up in either the pairwise or the all trait analyses, so I guess
# there's nothing more to do with this.


# Write results to text files
write.table(locus12all[[1]],"HyPrColoc_of_BMIadjT2D_and_all_QTLs_for_3q27.txt",sep = "\t",row.names = F,quote = F)
write.table(locus12coloc[[1]],"HyPrColoc_of_BMIadjT2D_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_3q27.txt",sep = "\t",row.names = F,quote = F)
write.table(locus20all[[1]],"HyPrColoc_of_BMIadjT2D_and_all_QTLs_for_7q32.txt",sep = "\t",row.names = F,quote = F)
write.table(locus20coloc[[1]],"HyPrColoc_of_BMIadjT2D_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_7q32.txt",sep = "\t",row.names = F,quote = F)
write.table(locus20fdr[[1]],"HyPrColoc_of_BMIadjT2D_and_all_FDR_sig_QTLs_for_7q32.txt",sep = "\t",row.names = F,quote = F)
write.table(locus42all[[1]],"HyPrColoc_of_BMIadjT2D_and_all_QTLs_for_20q13.32.txt",sep = "\t",row.names = F,quote = F)
write.table(locus42coloc[[1]],"HyPrColoc_of_BMIadjT2D_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_20q13.32.txt",sep = "\t",row.names = F,quote = F)


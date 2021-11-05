### This script is for the co-localization analyses of the selected loci between Triglycerides (TriG) GWAS and cis and trans eQTLs/aQTLs.

#install.packages("devtools")
#library(devtools)
#install_github("jrs95/hyprcoloc", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = F)
#browseVignettes("hyprcoloc") The install kept failing when trying to build the vignettes, so I disabled that.
#devtools::install_github("boxiangliu/locuscomparer")
library(hyprcoloc)
library(locuscomparer)

setwd("YOUR WORKING DIRECTORY")

# Read in the LD matrices, and the GWAS and QTL data. The file locations are relative to your working directory, so adjust accordingly.
gwas=read.table("jointGwasMc_TG.txt",sep = "\t",header = T)
cis_eqtl=read.table("./Eurobats_adipose_select_loci_cis-eQTLs_from_INT_logTPM.txt",sep = "\t",header = T)
cis_aqtl=read.table("./Eurobats_adipose_select_loci_cis-aQTLs_from_unnormalized_activities.txt",sep = "\t",header = T)
trans_eqtl=read.table("./Eurobats_adipose_select_loci_trans-eQTLs_for_TriG_MRs.txt",sep = "\t",header = T)
trans_aqtl=read.table("./Eurobats_adipose_select_loci_trans-aQTLs_for_TriG_MRs.txt",sep = "\t",header = T)
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

# The HDL GWAS data has coordinates for hg18 and hg 19, but I need to have CHR and POS columns (based on hg19) instead.
colnames(gwas)=c("CHR","POS","SNP","A1","A2","BETA","SE","N","P","Freq.A1.1000G.EUR")
gwas$CHR=gsub("chr","",gwas$CHR)
gwas$CHR=as.numeric(gsub(":.*","",gwas$CHR))
# This introduced NAs, but only for 5 SNPs without rsIDs (labeled only as ".")
gwas=gwas[!is.na(gwas$CHR),]
gwas$POS=as.numeric(gsub("chr.*:","",gwas$POS))

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
filt_cis_eqtl$chr=filt_gwas[match(filt_cis_eqtl$snps,filt_gwas$SNP),1]
filt_cis_aqtl$chr=filt_gwas[match(filt_cis_aqtl$snps,filt_gwas$SNP),1]
filt_trans_eqtl$chr=filt_gwas[match(filt_trans_eqtl$snps,filt_gwas$SNP),1]
filt_trans_aqtl$chr=filt_gwas[match(filt_trans_aqtl$snps,filt_gwas$SNP),1]
filt_cis_eqtl$position=filt_gwas[match(filt_cis_eqtl$snps,filt_gwas$SNP),2]
filt_cis_aqtl$position=filt_gwas[match(filt_cis_aqtl$snps,filt_gwas$SNP),2]
filt_trans_eqtl$position=filt_gwas[match(filt_trans_eqtl$snps,filt_gwas$SNP),2]
filt_trans_aqtl$position=filt_gwas[match(filt_trans_aqtl$snps,filt_gwas$SNP),2]

# Sort by chr and position
filt_gwas=filt_gwas[order(filt_gwas$CHR,filt_gwas$POS),]
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
  loci_betas[[i]]=as.data.frame(filt_gwas[filt_gwas$SNP %in% rownames(ld[[i]]),"BETA"],
                                row.names = as.character(filt_gwas[filt_gwas$SNP %in% rownames(ld[[i]]),"SNP"]))
  colnames(loci_betas[[i]])="TriG"
  
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
  colnames(loci_ses[[i]])="TriG"
  
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

# Pairwise co-localization, which is essentially equivalent with using coloc.
# Make a data.frame with only the QTLs co-localized with TriG.
GWAScoloc=list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43) # initialize a list of 43 objects
posColoc=data.frame() # initialize a data.frame for the positive co-localized traits
for(i in 1:43){
  # Grab the traits list for each locus
  traits=colnames(filt_betas[[i]])
  
  # Iterate through each of the QTLs. Exclude BMI from the index of traits since all will be tested against BMI.
  for(j in traits[-1]){
    GWAScoloc[[i]][[j]]=hyprcoloc(as.matrix(filt_betas[[i]]),as.matrix(filt_ses[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas[[i]]),ld.matrix = LD[[i]],
                                  trait.subset = c("TriG",j))
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
colnames(pair_tallies)=c("GWAS","Colocalizations with TriG")
# HyPrColoc reports colocalizations with PP>0.25, which is what is tallied here.
#Locus                                   GWAS            Colocalizations with TriG
#1p36.1                                   BMI                         0
#1p13.3                                   HDL                         0
#1q24_loc1                                BMI                         0
#1q24_loc2                                WHR                         0
#1q25.2                                   BMI                         0
#1q41                                     WHR                         0
#2p24                           Triglycerides                         0
#2p23.3                                   BMI                         0
#2p21                                     T2D                         0
#2q24.3                                   WHR                         0
#3p25                                     HDL                         0
#3q27                                     T2D                         0
#4q22_loc1                                BMI                         0
#4q22_loc2                                BMI                         0
#5q11.2                           HDL and T2D                         0
#5q31.2                                   BMI                         0
#6p21_loc1                      Triglycerides                         2
#6p21_loc2                                HDL                         0
#6q22.3                                   WHR                         1
#7q32         BMI, T2D, HDL and Triglycerides                        11
#7q36                                     HDL                         0
#8q21.1                                   BMI                         0
#8q21.2                                   BMI                         0
#10p13                                    BMI                         0
#10q22.2                                  BMI                         0
#10q24.2                                  BMI                         0
#10q26.3                                  BMI                         0
#11p11.2                                  HDL                         0
#11q13                                    BMI                         0
#12p13.33                                 BMI                         0
#12p13.1                                  BMI                         0
#12q13.13                                 BMI                         0
#15q21.3                                  HDL                         1
#15q24.1                                  BMI                         0
#17q21.2                                  BMI                         0
#17q25.3                                  HDL                         2
#18q21.1                                  HDL                         0
#18q21.3                                  BMI                         0
#19p13.1                T2D and Triglycerides                         0
#19q13.3_loc1                             HDL                         1
#19q13.3_loc2                             BMI                         3
#20q13.32                                 T2D                         0
#22q13.3                                  T2D                         0

### Mult-trait analyses

# 2p24, which is the 7th locus in the list
# Start with all traits included
locus7all=hyprcoloc(as.matrix(filt_betas[[7]]),as.matrix(filt_ses[[7]]),
                     trait.names=colnames(filt_betas[[7]]),snp.id=rownames(filt_betas[[7]]),ld.matrix = LD[[7]],
                     snpscores = T)
locus7all_cluster1=strsplit(locus7all[[1]][1,2],", ")[[1]]
locus7all_cluster20=strsplit(locus7all[[1]][20,2],", ")[[1]]
# Two clusters were found (1 and 20), but neither included TriG. The first cluster included 8 trans-eQTL and 12 trans-aQTLs with 
# a PP=0.4320 that was entirely explained by rs10166792 that has a TriG P=0.101. Cluster 20 included cis-e_APOB, 2 trans-eQTLs and 
# 13 trans-aQTLs with a PP=0.3923 that was best explained by rs56090741 with a TriG P=0.0065.

# 6p21_loc1,which was the 17th locus in the list
# Start with all traits included
locus17all=hyprcoloc(as.matrix(filt_betas[[17]]),as.matrix(filt_ses[[17]]),
                    trait.names=colnames(filt_betas[[17]]),snp.id=rownames(filt_betas[[17]]),ld.matrix = LD[[17]],
                    snpscores = T)
locus17all_cluster1=strsplit(locus17all[[1]][1,2],", ")[[1]]
locus17all_cluster2=strsplit(locus17all[[1]][2,2],", ")[[1]]
locus17all_cluster8=strsplit(locus17all[[1]][8,2],", ")[[1]]
# Four clusters were identified (1, 2, 3 and 8) and TriG went to cluster 3 with only cis-e_HCG27 (as in the pairwise analysis) with
# a PP=0.9566 that is entirely explained by rs1265098 with a TriG P=6.85e-13. Cluster 1 includes 5 cis-eQTLs, 6 trans-eQTLs and 
# 9 trans-aQTLs with a PP=0.4323 that is almost entirely explained by rs9266231 with a TriG P=1.92e-5. Cluster 2 includes 3 cis-eQTLs,
# cis-a_POU5F1, 6 trans-eQTLs and 5 trans-aQTLs with a PP=0.3355 that is entirely explained by rs3094006 with a TriG P=2.97e-5. Cluster
# 8 includes 5 cis-eQTLs, 2 cis-aQTLs, 1 trans-eQTL and 1 trans-aQTL with a PP=0.3720 that is entirely explained by rs6929796 with a
# TriG P=0.012.

# Let's run an analysis with only the QTLs that co-localized with TriG in the pairwise analysis.
coloc17_traits=c("TriG",gsub("TriG, ","",posColoc[posColoc$locus=="6p21_loc1",2]))
locus17coloc=hyprcoloc(as.matrix(filt_betas[[17]]),as.matrix(filt_ses[[17]]),
                      trait.names=colnames(filt_betas[[17]]),snp.id=rownames(filt_betas[[17]]),ld.matrix = LD[[17]],
                      trait.subset = coloc17_traits,snpscores = T)
# Just as in the all trait analysis above, the TriG only clustered with cis-e_HCG27. The other pairwise colocalized QTL was cis-eFLOT1
# with a PP=0.8132 that was almost entirely explained by rs3130557 with a TriG P=3.2e-14. This is odd given that the PP is pretty good
# for both pairwise colocalizations and both top candidate variants are strong TriG GWAS SNPs, though both are 7-8 orders of
# magnitude worse than the top TriG SNP at the locus. I checked the LD between rs1265098, rs3130557 and rs2247056 and the r^2 is not
# very good (0.3 between rs2247056 and 3130557 and less than 0.1 for the other 2 pairs). This is consistent with the fact that each
# colocalization here is entirely explained by a single SNP rather than spreading the explanation among correlated SNPs. But this also
# suggests there may be multiple independent TriG GWAS signals at this locus. This locus was of interest due to a sig cis-aQTL for
# TNF where the best aQTL SNP was also rs1265098, so I'm not sure why that wouldn't also colocalize with TriG or cluster with TriG and
# cis-e_HCG27 unless there is another stronger cis-a_TNF signal in this locus that doesn't match up with any TriG signal. I suppose
# some locus compare plots would be needed to get to the bottom of this. Note, TNF was not filtered out, so it is there in filt_betas.

# I doubt this will help, but let's try adding TNF to the pairwise colocalized QTLs.
coloc17_traitsTNF=append(coloc17_traits,"cis-a_TNF")
locus17colocTNF=hyprcoloc(as.matrix(filt_betas[[17]]),as.matrix(filt_ses[[17]]),
                       trait.names=colnames(filt_betas[[17]]),snp.id=rownames(filt_betas[[17]]),ld.matrix = LD[[17]],
                       trait.subset = coloc17_traitsTNF,snpscores = T)
# As expected, no clustering with cis-a_TNF.

# 7q32, which is the 20th locus in the list
# Start with all traits included
locus20all=hyprcoloc(as.matrix(filt_betas[[20]]),as.matrix(filt_ses[[20]]),
                    trait.names=colnames(filt_betas[[20]]),snp.id=rownames(filt_betas[[20]]),ld.matrix = LD[[20]],
                    snpscores = T)
locus20all_cluster1=strsplit(locus20all[[1]][1,2],", ")[[1]]
coloc20_traits=c("TriG",gsub("TriG, ","",posColoc[posColoc$locus=="7q32",2]))
coloc20_traits[(coloc20_traits %in% locus20all_cluster1)] # 12/12
# All 11 QTLs that pairwise colocalized with TriG in the pairwise analysis are included in cluster 1 that has TriG, 3 cis eQTLs, 
# cis-a_KLF14, 10 trans-eQTLs and 7 trans-aQTLs and a PP=0.3420 that is best explained by rs10954284 with a TriG P=1.83e-6, which is
# ~3 orders of magnitude worse than the best TriG SNP at the locus. There were 2 other clusters for this locus, which I've seen before 
# in the other HyPrColoc analyses with other traits, but they again seem irrelevant to the TriG GWAS with GWAS P=0.12 and P=0.13.

# Given that the all trait analysis generated a cluster with all of the pairwise colocalized traits, it would be redundant to run
# the analysis on just the pairwise colocalized traits.

# Finally for this locus, let's check just the FDR sig QTLs
fdr20_traits=c("TriG","cis-e_KLF14","cis-e_LINC-PINT","cis-a_KLF14","trans-e_AGT","trans-e_RABIF","trans-e_TBX4")
locus20fdr=hyprcoloc(as.matrix(filt_betas[[20]]),as.matrix(filt_ses[[20]]),
                     trait.names=colnames(filt_betas[[20]]),snp.id=rownames(filt_betas[[20]]),ld.matrix = LD[[20]],
                     trait.subset = fdr20_traits,snpscores = T)
# All of the FDR sig QTLs clustered with TriG, which is not surprising since each pairwise colocalized with TriG and were in the all
# trait analysis cluster with TriG. However, by restricting the analysis to the FDR sig QTLs the PP goes up (0.6639) and the cluster
# is now best explained by rs13241538 with a TriG P=1.77e-8, which is within 1 order of magnitude of the top TriG GWAS SNP in the locus.

# 19p13.1, which is the 39th locus in the list
# Start with all traits included
locus39all=hyprcoloc(as.matrix(filt_betas[[39]]),as.matrix(filt_ses[[39]]),
                     trait.names=colnames(filt_betas[[39]]),snp.id=rownames(filt_betas[[39]]),ld.matrix = LD[[39]],
                     snpscores = T)
locus39all_cluster1=strsplit(locus39all[[1]][1,2],", ")[[1]]
# TriG did not cluster with anything, but there was one cluster that included 9 trans-eQTLs and 35 trans-aQTLs with a PP=0.4519 that
# was best explained by rs56007657 with a TriG P=0.5759. So basically nothing exciting going on here.

# Write results to text files
write.table(posColoc,"Pairwise_HyPrColoc_between_TriG_and_each_QTL_for_select_loci.txt",sep = "\t",row.names = F,quote = F)
write.table(locus17all[[1]],"HyPrColoc_of_TriG_and_all_QTLs_for_6p21_loc1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus17coloc[[1]],"HyPrColoc_of_TriG_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_6p21_loc1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus20all[[1]],"HyPrColoc_of_TriG_and_all_QTLs_for_7q32.txt",sep = "\t",row.names = F,quote = F)
write.table(locus20fdr[[1]],"HyPrColoc_of_TriG_and_all_FDR_sig_QTLs_for_7q32.txt",sep = "\t",row.names = F,quote = F)


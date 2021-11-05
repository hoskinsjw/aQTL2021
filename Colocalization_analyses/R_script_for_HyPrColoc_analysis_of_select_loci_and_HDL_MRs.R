### This script is for the co-localization analyses of the selected loci between HDL GWAS and cis and trans eQTLs/aQTLs.

#install.packages("devtools")
#library(devtools)
#install_github("jrs95/hyprcoloc", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = F)
#browseVignettes("hyprcoloc") The install kept failing when trying to build the vignettes, so I disabled that.
#devtools::install_github("boxiangliu/locuscomparer")
library(hyprcoloc)
library(locuscomparer)

setwd("YOUR WORKING DIRECTORY")

# Read in the LD matrices, and the GWAS and QTL data. The file locations are relative to your working directory, so adjust accordingly.
gwas=read.table("jointGwasMc_HDL.txt",sep = "\t",header = T)
cis_eqtl=read.table("./Eurobats_adipose_select_loci_cis-eQTLs_from_INT_logTPM.txt",sep = "\t",header = T)
cis_aqtl=read.table("./Eurobats_adipose_select_loci_cis-aQTLs_from_unnormalized_activities.txt",sep = "\t",header = T)
trans_eqtl=read.table("./Eurobats_adipose_select_loci_trans-eQTLs_for_HDL_MRs.txt",sep = "\t",header = T)
trans_aqtl=read.table("./Eurobats_adipose_select_loci_trans-aQTLs_for_HDL_MRs.txt",sep = "\t",header = T)
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
  colnames(loci_betas[[i]])="HDL"
  
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
  colnames(loci_ses[[i]])="HDL"
  
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

# Locus 6p21_loc2 (the 18th) was selected because there were aQTLs for TEAD3 and FANCE that were stronger than the eQTLs and overlapped with 
# the HDL locus. However, the above filtering scheme loses both due to NAs. Rather than change the scheme for every locus, I
# will just do a separate special filter for this locus so I can compare the colocalizations with and without TEAD3 and FANCE.
gene_nas=c()
for(i in colnames(loci_betas[[18]])){
  gene_nas[i]=sum(is.na(loci_betas[[18]][,i]))
}
betas18=loci_betas[[18]][,gene_nas<(0.36*dim(loci_betas[[18]])[1])]
ses18=loci_ses[[18]][,gene_nas<(0.36*dim(loci_ses[[18]])[1])]

snp_nas=c()
for(i in rownames(betas18)){
  snp_nas[i]=sum(is.na(betas18[i,]))
}
# Turns out there is one QTL (cis-e_LEMD2) that has NAs with SNPs at the far downstream, while the bulk has NAs with SNPs at the far upstream.
# Therefore, if we drop that QTL we can hold onto 122 more SNPs.
betas18=betas18[snp_nas<=1,colnames(betas18)!="cis-e_LEMD2"]
ses18=ses18[snp_nas<=1,colnames(ses18)!="cis-e_LEMD2"]

LD18=ld[[18]][rownames(betas18),rownames(betas18)]


# Pairwise co-localization, which is essentially equivalent with using coloc.
# Make a data.frame with only the QTLs co-localized with HDL.
GWAScoloc=list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43) # initialize a list of 43 objects
posColoc=data.frame() # initialize a data.frame for the positive co-localized traits
for(i in 1:43){
  # Grab the traits list for each locus
  traits=colnames(filt_betas[[i]])
  
  # Iterate through each of the QTLs. Exclude BMI from the index of traits since all will be tested against BMI.
  for(j in traits[-1]){
    GWAScoloc[[i]][[j]]=hyprcoloc(as.matrix(filt_betas[[i]]),as.matrix(filt_ses[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas[[i]]),ld.matrix = LD[[i]],
                                  trait.subset = c("HDL",j))
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
colnames(pair_tallies)=c("GWAS","Colocalizations with HDL")
# HyPrColoc reports colocalizations with PP>0.25, which is what is tallied here.
#Locus                                   GWAS            Colocalizations with HDL
#1p36.1                                   BMI                        0
#1p13.3                                   HDL                        0
#1q24_loc1                                BMI                        0
#1q24_loc2                                WHR                        0
#1q25.2                                   BMI                        0
#1q41                                     WHR                        0
#2p24                           Triglycerides                        0
#2p23.3                                   BMI                        0
#2p21                                     T2D                        0
#2q24.3                                   WHR                        6
#3p25                                     HDL                        4
#3q27                                     T2D                        0
#4q22_loc1                                BMI                        0
#4q22_loc2                                BMI                        0
#5q11.2                           HDL and T2D                       10
#5q31.2                                   BMI                        0
#6p21_loc1                      Triglycerides                        2
#6p21_loc2                                HDL                        2
#6q22.3                                   WHR                        1
#7q32         BMI, T2D, HDL and Triglycerides                       14
#7q36                                     HDL                        2
#8q21.1                                   BMI                        0
#8q21.2                                   BMI                        0
#10p13                                    BMI                        0
#10q22.2                                  BMI                        0
#10q24.2                                  BMI                        0
#10q26.3                                  BMI                        0
#11p11.2                                  HDL                        0
#11q13                                    BMI                        0
#12p13.33                                 BMI                        0
#12p13.1                                  BMI                        0
#12q13.13                                 BMI                        0
#15q21.3                                  HDL                        1
#15q24.1                                  BMI                        0
#17q21.2                                  BMI                        0
#17q25.3                                  HDL                        4
#18q21.1                                  HDL                        1
#18q21.3                                  BMI                        0
#19p13.1                T2D and Triglycerides                        0
#19q13.3_loc1                             HDL                        0
#19q13.3_loc2                             BMI                        2
#20q13.32                                 T2D                        0
#22q13.3                                  T2D                        0

# Out of the 11 HDL loci, 8 had at least 1 pairwise colocalization. 

# Let's write the pairwise analysis to file for reference in the rest of these analyses.
write.table(posColoc,"Pairwise_HyPrColoc_between_HDL_and_each_QTL_for_select_loci.txt",sep = "\t",row.names = F,quote = F)

# Let's also run the pairwise analysis for 6p21_loc2 with TEAD3 QTLs included
pair18=list()
posColoc18=data.frame() # initialize a data.frame for the positive co-localized traits
# Iterate through each of the QTLs. Exclude BMI from the index of traits since all will be tested against BMI.
for(j in colnames(betas18)[-1]){
  pair18[[j]]=hyprcoloc(as.matrix(betas18),as.matrix(ses18),
                                trait.names=colnames(betas18),snp.id=rownames(betas18),ld.matrix = LD18,
                                trait.subset = c("HDL",j))
  # Grab only the QTLs with co-localization into posColoc
  if(!is.na(pair18[[j]]$results[1,3])){
    posColoc18=rbind(posColoc18,cbind(pair18[[j]]$results,"locus"="6p21_loc2"))
  }
}
# The same 2 QTLs (cis-e_C6orf106 and trans-e_NTRK3) pairwise co-localized given this alternative filtering.

# Let's write these pairwise results using the alternate filtering scheme for 6p21_loc2
write.table(posColoc18,"Pairwise_HyPrColoc_between_HDL_and_each_QTL_for_6p21_loc2_using_alternate_filtering_scheme.txt",sep = "\t",row.names = F,quote = F)

### Mult-trait analyses

# 3p25, which is the 11th locus in the list
# Start with all traits included
locus11all=hyprcoloc(as.matrix(filt_betas[[11]]),as.matrix(filt_ses[[11]]),
                  trait.names=colnames(filt_betas[[11]]),snp.id=rownames(filt_betas[[11]]),ld.matrix = LD[[11]],
                  snpscores = T)
locus11all_cluster1=strsplit(locus11all[[1]][1,2],", ")[[1]]
locus11all_cluster2=strsplit(locus11all[[1]][2,2],", ")[[1]]
locus11all_cluster4=strsplit(locus11all[[1]][4,2],", ")[[1]]
coloc11_traits=c("HDL",gsub("HDL, ","",posColoc[posColoc$locus=="3p25",2]))
coloc11_traits[(coloc11_traits %in% locus11all_cluster1)] # 2 trans-eQTLs
coloc11_traits[(coloc11_traits %in% locus11all_cluster2)] # 2 trans-eQTLs
coloc11_traits[(coloc11_traits %in% locus11all_cluster4)] # 2 trans-eQTLs
# Three clusters were found (1, 2 and 4), but none included HDL. Instead, the first included cis-a_VGLL4, 15 trans-eQTLs and 45
# trans-aQTLs with a PP=0.3831 that is best explained by rs11711175 that has an HDL GWAS P=0.23. The second includes cis-a_HRH1,
# 11 trans-eQTLs and 21 trans-aQTLs with a PP=0.4711 entirely explained by rs11719708 that has an HDL GWAS P=0.9387. The third
# includes 5 trans-eQTLs and 11 trans-aQTLs with a PP=0.4650 that is entirely explained by rs9682960 that has an HDL GWAS
# P=0.9016. Two of the 4 trans-eQTLs (BMP3 and RGS22) that pairwise colocalized with HDL went to cluster 1, which is particularly
# odd for RGS22 since its pairwise PP with HDL was 0.9278. If the pairwise colocalization was so strong, how does it end up in a
# cluster explained by a SNP that's completely irrelevant to the HDL signal?! By the way, RGS22 is the only FDR sig HDL MR
# trans-QTL for this locus, and it is only sig as an eQTL.

# Let's run an analysis with only the QTLs that co-localized with HDL in the pairwise analysis.
locus11coloc=hyprcoloc(as.matrix(filt_betas[[11]]),as.matrix(filt_ses[[11]]),
                       trait.names=colnames(filt_betas[[11]]),snp.id=rownames(filt_betas[[11]]),ld.matrix = LD[[11]],
                       trait.subset = coloc11_traits,snpscores = T)
# All 4 pairwise colocalized trans-eQTLs clustered together with HDL with a PP=0.489 almost entirely explained by rs2606736, which
# is the top GWAS SNP. 

# 5q11.2, which is the 15th locus in the list
# Start with all traits included
locus15all=hyprcoloc(as.matrix(filt_betas[[15]]),as.matrix(filt_ses[[15]]),
                     trait.names=colnames(filt_betas[[15]]),snp.id=rownames(filt_betas[[15]]),ld.matrix = LD[[15]],
                     snpscores = T)
locus15all_cluster1=strsplit(locus15all[[1]][1,2],", ")[[1]]
locus15all_cluster2=strsplit(locus15all[[1]][2,2],", ")[[1]]
locus15all_cluster6=strsplit(locus15all[[1]][6,2],", ")[[1]]
coloc15_traits=c("HDL",gsub("HDL, ","",posColoc[posColoc$locus=="5q11.2",2]))
coloc15_traits[(coloc15_traits %in% locus15all_cluster1)] # 2 trans-eQTLs
coloc15_traits[(coloc15_traits %in% locus15all_cluster2)] # 2 trans-eQTLs
coloc15_traits[(coloc15_traits %in% locus15all_cluster6)] # 2 trans-eQTLs
# Three clusters were found (1,2 and 6) and HDL landed in cluster 6 with only trans-e_JDP2 (also a pairwise colocalization) and a 
# PP=0.4541 that was best explained by rs6450176 (the top GWAS SNP). The first cluster was the largest with 29 trans-eQTLs and
# 81 trans-aQTLs with a PP=0.3714 that was entirely explained by rs1645273 with a GWAS P=0.4791. The second cluster contained
# cis-e_ARL15, cis-a_FST (which was a sig aQTL), 8 trans-eQTLs and 7 trans-aQTLs with a PP=0.3764 that was best explained
# by rs255766 with a GWAS P=0.0397. Four of the pairwise colocalized QTLs went into cluster 1, three went into cluster 2 and just
# the one went to cluster 6 with HDL. This locus was of interest because the HDL signal had a sig cis-aQTL with FST that was
# stronger than its eQTL.

# Let's run an analysis with only the QTLs that co-localized with HDL in the pairwise analysis.
locus15coloc=hyprcoloc(as.matrix(filt_betas[[15]]),as.matrix(filt_ses[[15]]),
                       trait.names=colnames(filt_betas[[15]]),snp.id=rownames(filt_betas[[15]]),ld.matrix = LD[[15]],
                       trait.subset = coloc15_traits,snpscores = T)
# As with the all trait analysis, only trans-e_JDP2 clustered with HDL while the rest clustered with each other in another cluster
# that had a PP=0.4086 that was best explained by rs255766 with a GWAS P=0.0397.

# 6p21_loc2, which is the 18th locus in the list. This locus had an alternate filtering scheme in addition to the standard in order
# to include TEAD3 that had an aQTL that was FDR sig and better than its eQTL. So I am testing both filtering schemes.
# Start with all traits included
locus18all=hyprcoloc(as.matrix(filt_betas[[18]]),as.matrix(filt_ses[[18]]),
                    trait.names=colnames(filt_betas[[18]]),snp.id=rownames(filt_betas[[18]]),ld.matrix = LD[[18]],
                    snpscores = T)
locus18all_cluster2=strsplit(locus18all[[1]][2,2],", ")[[1]]

locus18all_alt=hyprcoloc(as.matrix(betas18),as.matrix(ses18),
                     trait.names=colnames(betas18),snp.id=rownames(betas18),ld.matrix = LD18,
                     snpscores = T)
locus18all_alt_cluster2=strsplit(locus18all_alt[[1]][2,2],", ")[[1]]
# For both filtering schemes there was only one cluster identified that included HDL, cis-e_C6orf106, 3 trans-eQTLs and 1 trans-aQTL.
# The PP for the two different filtered datasets were very similar (0.40 and 0.41, respectively) and in both cases the top candidate
# variant was rs2744937 that has a GWAS P=1.8E-12, which is 1 order of magnitude worse than the top GWAS SNP. Note that this cluster
# included the 2 pairwise colocalized QTLs.

# Let's run an analysis with only the QTLs that co-localized with HDL in the pairwise analysis.
coloc18_traits=c("HDL",gsub("HDL, ","",posColoc[posColoc$locus=="6p21_loc2",2]))
locus18coloc=hyprcoloc(as.matrix(filt_betas[[18]]),as.matrix(filt_ses[[18]]),
                       trait.names=colnames(filt_betas[[18]]),snp.id=rownames(filt_betas[[18]]),ld.matrix = LD[[18]],
                       trait.subset = coloc18_traits,snpscores = T)
# Interestingly, despite relatively high pairwise colocalization with HDL and their inclusion in the one cluster for the all QTL
# analysis above, the cluster of cis-e_C6orf106 and trans-e_NTRK3 with HDL is much weaker than the pairwise or all QTL colocalizations.
# The cluster's PP=0.3041 that is best explained by rs2814944 that has a GWAS P=3.95E-13, which is very close to the top GWAS SNP. There
# was no difference in pairwise clustering between the two filtering schemes, so there was no need to test this in the alternate scheme.

# Since this locus was of interest due to the sig cis-aQTLs with TEAD3 and FANCE, let's see what happens when we 
# include those with the pairwise colocalized traits.
coloc18_traitsTF=c(coloc18_traits,"cis-a_TEAD3","cis-a_FANCE")
locus18colocTF=hyprcoloc(as.matrix(betas18),as.matrix(ses18),
                         trait.names=colnames(betas18),snp.id=rownames(betas18),ld.matrix = LD18,
                         trait.subset = coloc18_traitsTF,snpscores = T)
# Neither cis-a_TEAD3 nor cis-a_FANCE were included in the cluster with the pairwise colocalized QTLs.

# Finally for this locus, let's check just the FDR sig QTLs
fdr18_traits=c("HDL","cis-e_ZNF76","cis-a_TEAD3","cis-a_FANCE")
locus18fdr=hyprcoloc(as.matrix(betas18),as.matrix(ses18),
                         trait.names=colnames(betas18),snp.id=rownames(betas18),ld.matrix = LD18,
                         trait.subset = fdr18_traits,snpscores = T)
# No cluster.

# 7q32, which is the 20th locus in the list
# Start with all traits included
locus20all=hyprcoloc(as.matrix(filt_betas[[20]]),as.matrix(filt_ses[[20]]),
                    trait.names=colnames(filt_betas[[20]]),snp.id=rownames(filt_betas[[20]]),ld.matrix = LD[[20]],
                    snpscores = T)
locus20all_cluster1=strsplit(locus20all[[1]][1,2],", ")[[1]]
coloc20_traits=c("HDL",gsub("HDL, ","",posColoc[posColoc$locus=="7q32",2]))
coloc20_traits[(coloc20_traits %in% locus20all_cluster1)] # 14/14 pairwise colocalized QTLs
# All 14 QTLs that pairwise colocalized with HDL in the pairwise analysis are included in this cluster that has HDL, cis-e_LINC-PINT, 
# cis-e_KLF14, cis-e_AC016831.7, cis-a_KLF14, 14 trans-eQTLs and 8 trans-aQTLs with a PP=0.3803 that is best explained by rs10954284 
# that is the 4th best HDL GWAS SNP (P=5.80E-17). In the pairwise analysis, all but 2 of the colocalized QTLs were best explained 
# by rs3996352, which is the 2nd best HDL GWAS SNP (P=3.59E-17). In the pairwise analysis cis-e_KLF14 and trans-e_NR2F1 were best 
# explained by the 3rd best HDL GWAS SNP rs4731702. There was 1 other cluster for this locus, which I've seen before in the other 
# HyPrColoc analyses with other traits, but it again seem irrelevant to the HDL GWAS.

# Let's run an analysis with only the QTLs that co-localized with HDL in the pairwise analysis.
locus20coloc=hyprcoloc(as.matrix(filt_betas[[20]]),as.matrix(filt_ses[[20]]),
                      trait.names=colnames(filt_betas[[20]]),snp.id=rownames(filt_betas[[20]]),ld.matrix = LD[[20]],
                      trait.subset = coloc20_traits,snpscores = T)
locus20coloc_cluster1=strsplit(locus20coloc[[1]][1,2],", ")[[1]]
# All of the pairwise colocalizing QTLs cluster together when restricting the multi-trait analysis to the pairwise colocalizing traits.
# The PP=0.6622, which is a quite a bit better than in the all trait analysis, but is now best explained by rs3996352 (2nd best GWAS SNP), 
# which is consistent with the pairwise analysis.

# Finally for this locus, let's check just the FDR sig QTLs
fdr20_traits=c("HDL","cis-e_LINC-PINT","cis-e_KLF14","cis-a_KLF14","trans-e_ESR2","trans-e_NR2F1","trans-e_TBX4")
locus20fdr=hyprcoloc(as.matrix(filt_betas[[20]]),as.matrix(filt_ses[[20]]),
                     trait.names=colnames(filt_betas[[20]]),snp.id=rownames(filt_betas[[20]]),ld.matrix = LD[[20]],
                     trait.subset = fdr20_traits,snpscores = T)
# This yielded an impressive cluster that included all FDR sig traits and HDL with PP=0.8346 that is again best explained by 
# rs3996352 (2nd best GWAS SNP).

# 7q36, which is the 21st locus in the list
# Start with all traits included
locus21all=hyprcoloc(as.matrix(filt_betas[[21]]),as.matrix(filt_ses[[21]]),
                     trait.names=colnames(filt_betas[[21]]),snp.id=rownames(filt_betas[[21]]),ld.matrix = LD[[21]],
                     snpscores = T)
# There is one cluster of cis-e_TMEM176A and cis-e_TMEM167B that has a really high PP=0.995, but HDL does not cluster.

# Let's run an analysis with only the QTLs that co-localized with HDL in the pairwise analysis.
coloc21_traits=c("HDL",gsub("HDL, ","",posColoc[posColoc$locus=="7q36",2]))
locus21coloc=hyprcoloc(as.matrix(filt_betas[[21]]),as.matrix(filt_ses[[21]]),
                       trait.names=colnames(filt_betas[[21]]),snp.id=rownames(filt_betas[[21]]),ld.matrix = LD[[21]],
                       trait.subset = coloc21_traits,snpscores = T)
# Both pairwise colocalized QTLs (cis-e_GIMAP8 and trans-e_IL6R) do cluster with HDL with a PP=0.5746 that is almost entirely
# explained by rs13225097 with a GWAS P=4.32E-8, which is less than an order of magnitude worse than the top GWAS SNP.

# Finally for this locus, let's check just the FDR sig QTLs
fdr21_traits=c("HDL","cis-e_FASTK","cis-a_CDK5","cis-a_FASTK")
locus21fdr=hyprcoloc(as.matrix(filt_betas[[21]]),as.matrix(filt_ses[[21]]),
                     trait.names=colnames(filt_betas[[21]]),snp.id=rownames(filt_betas[[21]]),ld.matrix = LD[[21]],
                     trait.subset = fdr21_traits,snpscores = T)
# None of the FDR sig QTLs colocalize with HDL.

# 15q21.3, which is the 33rd locus in the list.
# Start with all traits included
locus33all=hyprcoloc(as.matrix(filt_betas[[33]]),as.matrix(filt_ses[[33]]),
                     trait.names=colnames(filt_betas[[33]]),snp.id=rownames(filt_betas[[33]]),ld.matrix = LD[[33]],
                     snpscores = T)
locus33all_cluster2=strsplit(locus33all[[1]][2,2],", ")[[1]]
locus33all_cluster12=strsplit(locus33all[[1]][12,2],", ")[[1]]
coloc33_traits=c("HDL",gsub("HDL, ","",posColoc[posColoc$locus=="15q21.3",2]))
coloc33_traits[(coloc33_traits %in% locus33all_cluster2)] # 0
coloc33_traits[(coloc33_traits %in% locus33all_cluster12)] # 0
# Two cluster were found, but neither included HDL. Oddly, neither included the one QTL (trans-a_IFT27) that colocalized with 
# HDL in the pairwise analysis with a decent PP=0.5025. If neither HDL nor trans-a_IFT27 clustered with anything else, why
# wouldn't they cluster with each other?

# 17q25.3, which is the 36th locus in the list
# Start with all traits included
locus36all=hyprcoloc(as.matrix(filt_betas[[36]]),as.matrix(filt_ses[[36]]),
                     trait.names=colnames(filt_betas[[36]]),snp.id=rownames(filt_betas[[36]]),ld.matrix = LD[[36]],
                     snpscores = T)
locus36all_cluster1=strsplit(locus36all[[1]][1,2],", ")[[1]]
locus36all_cluster2=strsplit(locus36all[[1]][2,2],", ")[[1]]
locus36all_cluster33=strsplit(locus36all[[1]][33,2],", ")[[1]]
coloc36_traits=c("HDL",gsub("HDL, ","",posColoc[posColoc$locus=="17q25.3",2]))
coloc36_traits[(coloc36_traits %in% locus36all_cluster1)] # 1/4
coloc36_traits[(coloc36_traits %in% locus36all_cluster2)] # 0/4
coloc36_traits[(coloc36_traits %in% locus36all_cluster33)] # 3/4
# Three clusters were found. The first contains 1 cis-eQTL, 7 trans-eQTLs and 4 trans-aQTLs (including trans-a_NR2F1 that
# pairwise colocalized with HDL) with a PP=0.4530 that is almost entirely exlained by rs744120 that has a GWAS P=0.7506. The
# second included 3 cis-eQTLs, 5 trans-eQTLs and 18 trans-aQTLs with a PP=0.4137 that was best explained by rs7215260 that
# has a GWAS P=0.336. The third cluster (33) does include HDL along wth cis-e_PGS1, 3 trans-eQTLs and 1 trans-aQTL (3 of these
# paiwise colocalized with HDL) with a PP=0.4861 that is almost entirely explained by rs4969178 that is the top HDL GWAS SNP for 
# the locus. This cluster 33 includes the only FDR sig QTL for this locus (trans-e_CAPN5).

# Let's run an analysis with only the QTLs that co-localized with HDL in the pairwise analysis.
locus36coloc=hyprcoloc(as.matrix(filt_betas[[36]]),as.matrix(filt_ses[[36]]),
                       trait.names=colnames(filt_betas[[36]]),snp.id=rownames(filt_betas[[36]]),ld.matrix = LD[[36]],
                       trait.subset = coloc36_traits,snpscores = T)
# All 4 pairwise colocalized QTLs also cluster together with HDL with a PP=0.6629 that is almost entirely explained by rs4969178
# that is the top HDL GWAS SNP for the locus. This includes the only FDR sig QTL for this locus (trans-e_CAPN5).

# 18q21.1, which is the 37th locus in the list
# Start with all traits included
locus37all=hyprcoloc(as.matrix(filt_betas[[37]]),as.matrix(filt_ses[[37]]),
                     trait.names=colnames(filt_betas[[37]]),snp.id=rownames(filt_betas[[37]]),ld.matrix = LD[[37]],
                     snpscores = T)
# Three clusters were found. The first contains 4 trans-eQTL and 8 trans-aQTL (no HDL) with a PP=0.4543 that is almost entirely
# explained by rs17713794 that has a GWAS P=0.0013. The second is the pairwise colocalization of HDL and trans-e_ZNF69 with a
# PP=0.4421 that is best explained by rs4939883 that is the top HDL GWAS SNP for the locus. The third cluster (31) includes
# 1 trans-eQTL and 8 trans-aQTLs with a PP=0.4586 that is best explained by rs12606416 that has a GWAS P=3.64E-5 (which is
# 51 orders of magnitude worse than the top HDL SNP for this locus). This locus was of interest due to a cis-aQTL for MBD1 that
# was stronger than its eQTL, but I'm seeing no evidence of that signal colocalizing with the GWAS or any other signals.


# Write results to text files
write.table(locus11all[[1]],"HyPrColoc_of_HDL_and_all_QTLs_for_3q25.txt",sep = "\t",row.names = F,quote = F)
write.table(locus11coloc[[1]],"HyPrColoc_of_HDL_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_3q25.txt",sep = "\t",row.names = F,quote = F)
write.table(locus15all[[1]],"HyPrColoc_of_HDL_and_all_QTLs_for_5q11.2.txt",sep = "\t",row.names = F,quote = F)
write.table(locus15coloc[[1]],"HyPrColoc_of_HDL_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_5q11.2.txt",sep = "\t",row.names = F,quote = F)
write.table(locus18all[[1]],"HyPrColoc_of_HDL_and_all_QTLs_for_6p21_loc2.txt",sep = "\t",row.names = F,quote = F)
write.table(locus18coloc[[1]],"HyPrColoc_of_HDL_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_6p21_loc2.txt",sep = "\t",row.names = F,quote = F)
write.table(locus20all[[1]],"HyPrColoc_of_HDL_and_all_QTLs_for_7q32.txt",sep = "\t",row.names = F,quote = F)
write.table(locus20coloc[[1]],"HyPrColoc_of_HDL_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_7q32.txt",sep = "\t",row.names = F,quote = F)
write.table(locus20fdr[[1]],"HyPrColoc_of_HDL_and_all_FDR_sig_QTLs_for_7q32.txt",sep = "\t",row.names = F,quote = F)
write.table(locus21all[[1]],"HyPrColoc_of_HDL_and_all_QTLs_for_7q36.txt",sep = "\t",row.names = F,quote = F)
write.table(locus21coloc[[1]],"HyPrColoc_of_HDL_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_7q36.txt",sep = "\t",row.names = F,quote = F)
write.table(locus36all[[1]],"HyPrColoc_of_HDL_and_all_QTLs_for_17q25.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus36coloc[[1]],"HyPrColoc_of_HDL_and_all_QTLs_that_colocalized_in_pairwise_analysis_for_17q25.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus37all[[1]],"HyPrColoc_of_HDL_and_all_QTLs_for_18q21.1.txt",sep = "\t",row.names = F,quote = F)


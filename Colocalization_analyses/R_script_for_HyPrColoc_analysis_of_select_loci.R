### This script is for the pairwise co-localization analyses of the selected loci between signals from all 5 GWAS
### and to check for multi-signal colocalizations between all GWAS and cis-eQTLs/aQTLs for the selected loci. 
### The full pairwise co-localizations between individual GWAS studies and all cis-eQTLs/aQTLs and matching MR
### trans-eQTLs/aQTLs are run with a separate script per GWAS all for the same set of selected loci.

install.packages("devtools")
library(devtools)
install_github("jrs95/hyprcoloc", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = F)
devtools::install_github("boxiangliu/locuscomparer")
library(hyprcoloc)
library(locuscomparer)

setwd("YOUR WORKING DIRECTORY")

# Read in the LD matrices, and the GWAS and QTL data. The file locations are relative to your working directory, so adjust accordingly.
bmi=read.table("./Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt",sep = "\t",header = T)
whr=read.table("./GIANT_2015_WHRadjBMI_COMBINED_EUR.txt",sep = "\t",header = T)
t2d=read.table("./Mahajan.NatGenet2018b.T2Dbmiadj.European.with.rsIDs.txt",sep = "\t",header = T)
hdl=read.table("./jointGwasMc_HDL.txt",sep = "\t",header = T)
triG=read.table("./jointGwasMc_TG.txt",sep = "\t",header = T)
cis_eqtl=read.table("./Eurobats_adipose_select_loci_cis-eQTLs_from_INT_logTPM.txt",sep = "\t",header = T)
cis_aqtl=read.table("./Eurobats_adipose_select_loci_cis-aQTLs_from_unnormalized_activities.txt",sep = "\t",header = T)
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

# Add in the BMI chr17p13.2 locus that was missing from the original list.
chr17p13_eqtl=read.table("./Eurobats_adipose_chr17p13.2_cis-eQTLs_from_INT_logTPM.txt",sep = "\t",header = T)
chr17p13_aqtl=read.table("./Eurobats_adipose_chr17p13.2_cis-aQTLs_from_unnormalized_activities.txt",sep = "\t",header = T)
ld17p13=read.table("./Eurobats_chr17p13.2_LD_matrix.txt",sep = "\t",header = F)
rownames(ld17p13)=ld17p13[,3]
ld17p13=ld17p13[,-c(1:5)]
colnames(ld17p13)=rownames(ld17p13)

# Before filtering, let's standardize the gwas column names a bit
colnames(whr)=c("SNP","Allele1","Allele2","FreqAlleleHapMapCEU","BETA","SE","p","N")

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

# Filter GWAS and QTL data to the same SNPs
filt_bmi=bmi[na.omit(match(unique(cis_eqtl$snps),bmi$SNP)),]
filt_bmi=filt_bmi[na.omit(match(whr$SNP,filt_bmi$SNP)),]
filt_bmi=filt_bmi[na.omit(match(t2d$rsID,filt_bmi$SNP)),]
filt_bmi=filt_bmi[na.omit(match(hdl$SNP,filt_bmi$SNP)),]
filt_bmi=filt_bmi[na.omit(match(triG$SNP,filt_bmi$SNP)),]
filt_whr=whr[na.omit(match(filt_bmi$SNP,whr$SNP)),]
filt_t2d=t2d[na.omit(match(filt_bmi$SNP,t2d$rsID)),]
filt_hdl=hdl[na.omit(match(filt_bmi$SNP,hdl$SNP)),]
filt_triG=triG[na.omit(match(filt_bmi$SNP,triG$SNP)),]
filt_cis_eqtl=cis_eqtl[cis_eqtl$snps %in% filt_bmi$SNP,]
filt_cis_aqtl=cis_aqtl[cis_aqtl$snps %in% filt_bmi$SNP,]

filt_bmi17p13=bmi[na.omit(match(unique(chr17p13_eqtl$snps),bmi$SNP)),]
filt_bmi17p13=filt_bmi17p13[na.omit(match(unique(whr$SNP),filt_bmi17p13$SNP)),]
filt_bmi17p13=filt_bmi17p13[na.omit(match(unique(t2d$rsID),filt_bmi17p13$SNP)),]
filt_bmi17p13=filt_bmi17p13[na.omit(match(unique(hdl$SNP),filt_bmi17p13$SNP)),]
filt_bmi17p13=filt_bmi17p13[na.omit(match(unique(triG$SNP),filt_bmi17p13$SNP)),]
filt_whr17p13=whr[na.omit(match(filt_bmi17p13$SNP,whr$SNP)),]
filt_t2d17p13=t2d[na.omit(match(filt_bmi17p13$SNP,t2d$rsID)),]
filt_hdl17p13=hdl[na.omit(match(filt_bmi17p13$SNP,hdl$SNP)),]
filt_triG17p13=triG[na.omit(match(filt_bmi17p13$SNP,triG$SNP)),]
filt_cis_eqtl17p13=chr17p13_eqtl[chr17p13_eqtl$snps %in% filt_bmi17p13$SNP,]
filt_cis_aqtl17p13=chr17p13_aqtl[chr17p13_aqtl$snps %in% filt_bmi17p13$SNP,]

# The WHRadjBMI GWAS results do not include SNP chr or position, so I need to grab them from the info file for the select loci.
filt_whr$CHR=filt_bmi[na.omit(match(filt_whr$SNP,filt_bmi$SNP)),"CHR"]
filt_whr$POS=filt_bmi[na.omit(match(filt_whr$SNP,filt_bmi$SNP)),"POS"]

filt_whr17p13$CHR=filt_bmi17p13[na.omit(match(filt_whr17p13$SNP,filt_bmi17p13$SNP)),"CHR"]
filt_whr17p13$POS=filt_bmi17p13[na.omit(match(filt_whr17p13$SNP,filt_bmi17p13$SNP)),"POS"]

# Calculate the SE of the Beta for the QTLs with the formula SE=beta/statistic
filt_cis_eqtl$SE=filt_cis_eqtl$beta/filt_cis_eqtl$statistic
filt_cis_aqtl$SE=filt_cis_aqtl$beta/filt_cis_aqtl$statistic

filt_cis_eqtl17p13$SE=filt_cis_eqtl17p13$beta/filt_cis_eqtl17p13$statistic
filt_cis_aqtl17p13$SE=filt_cis_aqtl17p13$beta/filt_cis_aqtl17p13$statistic

# Add chromosome and position to the QTLs for sorting
filt_cis_eqtl$chr=filt_bmi[match(filt_cis_eqtl$snps,filt_bmi$SNP),"CHR"]
filt_cis_aqtl$chr=filt_bmi[match(filt_cis_aqtl$snps,filt_bmi$SNP),"CHR"]
filt_cis_eqtl$position=filt_bmi[match(filt_cis_eqtl$snps,filt_bmi$SNP),"POS"]
filt_cis_aqtl$position=filt_bmi[match(filt_cis_aqtl$snps,filt_bmi$SNP),"POS"]

filt_cis_eqtl17p13$chr=filt_bmi17p13[match(filt_cis_eqtl17p13$snps,filt_bmi17p13$SNP),"CHR"]
filt_cis_aqtl17p13$chr=filt_bmi17p13[match(filt_cis_aqtl17p13$snps,filt_bmi17p13$SNP),"CHR"]
filt_cis_eqtl17p13$position=filt_bmi17p13[match(filt_cis_eqtl17p13$snps,filt_bmi17p13$SNP),"POS"]
filt_cis_aqtl17p13$position=filt_bmi17p13[match(filt_cis_aqtl17p13$snps,filt_bmi17p13$SNP),"POS"]

# Sort by chr and position
filt_bmi=filt_bmi[order(filt_bmi$CHR,filt_bmi$POS),]
filt_whr=filt_whr[order(filt_whr$CHR,filt_whr$POS),]
filt_t2d=filt_t2d[order(filt_t2d$Chr,filt_t2d$Pos),]
filt_hdl=filt_hdl[order(filt_hdl$CHR,filt_hdl$POS),]
filt_triG=filt_triG[order(filt_triG$CHR,filt_triG$POS),]
filt_cis_eqtl=filt_cis_eqtl[order(filt_cis_eqtl$chr,filt_cis_eqtl$position),]
filt_cis_aqtl=filt_cis_aqtl[order(filt_cis_aqtl$chr,filt_cis_aqtl$position),]
all(as.character(filt_bmi$SNP)==as.character(filt_whr$SNP)) # TRUE
all(as.character(filt_bmi$SNP)==as.character(filt_t2d$rsID)) # TRUE
all(as.character(filt_bmi$SNP)==as.character(filt_hdl$SNP)) # TRUE
all(as.character(filt_bmi$SNP)==as.character(filt_triG$SNP)) # TRUE

filt_bmi17p13=filt_bmi17p13[order(filt_bmi17p13$CHR,filt_bmi17p13$POS),]
filt_whr17p13=filt_whr17p13[order(filt_whr17p13$CHR,filt_whr17p13$POS),]
filt_t2d17p13=filt_t2d17p13[order(filt_t2d17p13$Chr,filt_t2d17p13$Pos),]
filt_hdl17p13=filt_hdl17p13[order(filt_hdl17p13$CHR,filt_hdl17p13$POS),]
filt_triG17p13=filt_triG17p13[order(filt_triG17p13$CHR,filt_triG17p13$POS),]
filt_cis_eqtl17p13=filt_cis_eqtl17p13[order(filt_cis_eqtl17p13$chr,filt_cis_eqtl17p13$position),]
filt_cis_aqtl17p13=filt_cis_aqtl17p13[order(filt_cis_aqtl17p13$chr,filt_cis_aqtl17p13$position),]
all(as.character(filt_bmi17p13$SNP)==as.character(filt_whr17p13$SNP)) # TRUE
all(as.character(filt_bmi17p13$SNP)==as.character(filt_t2d17p13$rsID)) # TRUE
all(as.character(filt_bmi17p13$SNP)==as.character(filt_hdl17p13$SNP)) # TRUE
all(as.character(filt_bmi17p13$SNP)==as.character(filt_triG17p13$SNP)) # TRUE

# Split GWAS and QTL Betas and SEs by locus, and also by genes for the QTL data. This section is pretty complex and hard to follow,
# for which I am both proud and embarrassed.

# Betas
loci_betas=list()
for(i in 1:43){
  # Grab the genes per locus separately for cis-eQTL, cis-aQTL, trans-eQTL and trans-aQTL
  locus_genes=list()
  locus_genes[[1]]=as.character(unique(filt_cis_eqtl[filt_cis_eqtl$snps %in% rownames(ld[[i]]),"gene"]))
  locus_genes[[2]]=as.character(unique(filt_cis_aqtl[filt_cis_aqtl$snps %in% rownames(ld[[i]]),"gene"]))
  
  # Begin by grabbing all of the GWAS betas per locus
  loci_betas[[i]]=data.frame("BMI"=filt_bmi[filt_bmi$SNP %in% rownames(ld[[i]]),"BETA"],
                             "WHR"=filt_whr[filt_whr$SNP %in% rownames(ld[[i]]),"BETA"],
                             "T2D"=filt_t2d[filt_t2d$rsID %in% rownames(ld[[i]]),"Beta"],
                             "HDL"=filt_hdl[filt_hdl$SNP %in% rownames(ld[[i]]),"BETA"],
                             "TriG"=filt_triG[filt_triG$SNP %in% rownames(ld[[i]]),"BETA"],
                             row.names = as.character(filt_bmi[filt_bmi$SNP %in% rownames(ld[[i]]),"SNP"]))

  # Next we need to grab all of the QTL data per cis-e/aQTL
  for(j in 1:2){
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
    }
  }
}

loci_betas17p13=list()
for(i in 1){
  # Grab the genes per locus separately for cis-eQTL, cis-aQTL, trans-eQTL and trans-aQTL
  locus_genes=list()
  locus_genes[[1]]=as.character(unique(filt_cis_eqtl17p13[filt_cis_eqtl17p13$snps %in% rownames(ld17p13),"gene"]))
  locus_genes[[2]]=as.character(unique(filt_cis_aqtl17p13[filt_cis_aqtl17p13$snps %in% rownames(ld17p13),"gene"]))
  
  # Begin by grabbing all of the GWAS betas per locus
  loci_betas17p13[[i]]=data.frame("BMI"=filt_bmi17p13[filt_bmi17p13$SNP %in% rownames(ld17p13),"BETA"],
                             "WHR"=filt_whr17p13[filt_whr17p13$SNP %in% rownames(ld17p13),"BETA"],
                             "T2D"=filt_t2d17p13[filt_t2d17p13$rsID %in% rownames(ld17p13),"Beta"],
                             "HDL"=filt_hdl17p13[filt_hdl17p13$SNP %in% rownames(ld17p13),"BETA"],
                             "TriG"=filt_triG17p13[filt_triG17p13$SNP %in% rownames(ld17p13),"BETA"],
                             row.names = as.character(filt_bmi17p13[filt_bmi17p13$SNP %in% rownames(ld17p13),"SNP"]))
  
  # Next we need to grab all of the QTL data per cis-e/aQTL
  for(j in 1:2){
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
  
  # Begin by grabbing all of the GWAS betas per locus
  loci_ses[[i]]=data.frame("BMI"=filt_bmi[filt_bmi$SNP %in% rownames(ld[[i]]),"SE"],
                             "WHR"=filt_whr[filt_whr$SNP %in% rownames(ld[[i]]),"SE"],
                             "T2D"=filt_t2d[filt_t2d$rsID %in% rownames(ld[[i]]),"SE"],
                             "HDL"=filt_hdl[filt_hdl$SNP %in% rownames(ld[[i]]),"SE"],
                             "TriG"=filt_triG[filt_triG$SNP %in% rownames(ld[[i]]),"SE"],
                             row.names = as.character(filt_bmi[filt_bmi$SNP %in% rownames(ld[[i]]),"SNP"]))
  
  # Next we need to grab all of the QTL data per cis-e/aQTL
  for(j in 1:2){
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
    }
  }
}

loci_ses17p13=list()
for(i in 1){
  # Grab the genes per locus separately for cis-eQTL, cis-aQTL, trans-eQTL and trans-aQTL
  locus_genes=list()
  locus_genes[[1]]=as.character(unique(filt_cis_eqtl17p13[filt_cis_eqtl17p13$snps %in% rownames(ld17p13),"gene"]))
  locus_genes[[2]]=as.character(unique(filt_cis_aqtl17p13[filt_cis_aqtl17p13$snps %in% rownames(ld17p13),"gene"]))
  
  # Begin by grabbing all of the GWAS ses per locus
  loci_ses17p13[[i]]=data.frame("BMI"=filt_bmi17p13[filt_bmi17p13$SNP %in% rownames(ld17p13),"SE"],
                                  "WHR"=filt_whr17p13[filt_whr17p13$SNP %in% rownames(ld17p13),"SE"],
                                  "T2D"=filt_t2d17p13[filt_t2d17p13$rsID %in% rownames(ld17p13),"SE"],
                                  "HDL"=filt_hdl17p13[filt_hdl17p13$SNP %in% rownames(ld17p13),"SE"],
                                  "TriG"=filt_triG17p13[filt_triG17p13$SNP %in% rownames(ld17p13),"SE"],
                                  row.names = as.character(filt_bmi17p13[filt_bmi17p13$SNP %in% rownames(ld17p13),"SNP"]))
  
  # Next we need to grab all of the QTL data per cis-e/aQTL
  for(j in 1:2){
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
    }
  }
}

# The SNPs and/or genes with NAs need to be filtered out.
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

# Pairwise co-localization, which is essentially equivalent to using coloc.
# The cis-QTL pairwise analyses have already been done by other scripts (though with somewhat different SNPs included),
# so here let's just do the pairwise colocalization among the GWAS.
GWAScoloc=list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43) # initialize a list of 43 objects
posColoc=data.frame() # initialize a data.frame for the positive co-localized traits
for(i in 1:43){
  # Grab the traits list for each locus
  traits=colnames(filt_betas[[i]])

  # Iterate through keeping each GWAS constant, without duplicating pairs
  # BMI constant
  for(j in traits[2:5]){
    GWAScoloc[[i]][[j]]=hyprcoloc(as.matrix(filt_betas[[i]]),as.matrix(filt_ses[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas[[i]]),ld.matrix = LD[[i]],
                                  trait.subset = c("BMI",j))
    # Grab only the QTLs with co-localization into posColoc
    if(!is.na(GWAScoloc[[i]][[j]]$results[1,3])){
      posColoc=rbind(posColoc,cbind(GWAScoloc[[i]][[j]]$results,"locus"=i))
    }
  }

  # WHR constant
  for(j in traits[3:5]){
    GWAScoloc[[i]][[j]]=hyprcoloc(as.matrix(filt_betas[[i]]),as.matrix(filt_ses[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas[[i]]),ld.matrix = LD[[i]],
                                  trait.subset = c("WHR",j))
    # Grab only the QTLs with co-localization into posColoc
    if(!is.na(GWAScoloc[[i]][[j]]$results[1,3])){
      posColoc=rbind(posColoc,cbind(GWAScoloc[[i]][[j]]$results,"locus"=i))
    }
  }

  # T2D constant
  for(j in traits[4:5]){
    GWAScoloc[[i]][[j]]=hyprcoloc(as.matrix(filt_betas[[i]]),as.matrix(filt_ses[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas[[i]]),ld.matrix = LD[[i]],
                                  trait.subset = c("T2D",j))
    # Grab only the QTLs with co-localization into posColoc
    if(!is.na(GWAScoloc[[i]][[j]]$results[1,3])){
      posColoc=rbind(posColoc,cbind(GWAScoloc[[i]][[j]]$results,"locus"=i))
    }
  }

  # HDL constant
  for(j in traits[5]){
    GWAScoloc[[i]][[j]]=hyprcoloc(as.matrix(filt_betas[[i]]),as.matrix(filt_ses[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas[[i]]),ld.matrix = LD[[i]],
                                  trait.subset = c("HDL",j))
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
  
  # Iterate through keeping each GWAS constant, without duplicating pairs
  # BMI constant
  for(j in traits[2:5]){
    GWAScoloc17p13[[i]][[j]]=hyprcoloc(as.matrix(filt_betas17p13[[i]]),as.matrix(filt_ses17p13[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas17p13[[i]]),ld.matrix = LD17p13[[i]],
                                  trait.subset = c("BMI",j))
    # Grab only the QTLs with co-localization into posColoc
    if(!is.na(GWAScoloc17p13[[i]][[j]]$results[1,3])){
      posColoc17p13=rbind(posColoc17p13,cbind(GWAScoloc17p13[[i]][[j]]$results,"locus"=i))
    }
  }
  
  # WHR constant
  for(j in traits[3:5]){
    GWAScoloc17p13[[i]][[j]]=hyprcoloc(as.matrix(filt_betas17p13[[i]]),as.matrix(filt_ses17p13[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas17p13[[i]]),ld.matrix = LD17p13[[i]],
                                  trait.subset = c("WHR",j))
    # Grab only the QTLs with co-localization into posColoc
    if(!is.na(GWAScoloc17p13[[i]][[j]]$results[1,3])){
      posColoc17p13=rbind(posColoc17p13,cbind(GWAScoloc17p13[[i]][[j]]$results,"locus"=i))
    }
  }
  
  # T2D constant
  for(j in traits[4:5]){
    GWAScoloc17p13[[i]][[j]]=hyprcoloc(as.matrix(filt_betas17p13[[i]]),as.matrix(filt_ses17p13[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas17p13[[i]]),ld.matrix = LD17p13[[i]],
                                  trait.subset = c("T2D",j))
    # Grab only the QTLs with co-localization into posColoc
    if(!is.na(GWAScoloc17p13[[i]][[j]]$results[1,3])){
      posColoc17p13=rbind(posColoc17p13,cbind(GWAScoloc17p13[[i]][[j]]$results,"locus"=i))
    }
  }
  
  # HDL constant
  for(j in traits[5]){
    GWAScoloc17p13[[i]][[j]]=hyprcoloc(as.matrix(filt_betas17p13[[i]]),as.matrix(filt_ses17p13[[i]]),
                                  trait.names=traits,snp.id=rownames(filt_betas17p13[[i]]),ld.matrix = LD17p13[[i]],
                                  trait.subset = c("HDL",j))
    # Grab only the QTLs with co-localization into posColoc
    if(!is.na(GWAScoloc17p13[[i]][[j]]$results[1,3])){
      posColoc17p13=rbind(posColoc17p13,cbind(GWAScoloc17p13[[i]][[j]]$results,"locus"=i))
    }
  }
}

# Let's change the locus number to legit chromosomal locus ID
loci=read.table("Select_loci_list.txt",sep = "\t",header = T)
posColoc$locus=loci[match(posColoc$locus,rownames(loci)),1]

posColoc17p13$locus=rep("17p13.2",dim(posColoc17p13)[1])

# Let's tally the co-localized QTLs per locus
pair_tallies=data.frame(row.names = loci$Locus)
for(i in loci$Locus){
  pair_tallies[i,1]=loci[loci$Locus==i,2]
  pair_tallies[i,2]=dim(posColoc[posColoc$locus==i,])[1]
}
colnames(pair_tallies)=c("Sig GWAS","Pairwise GWAS colocalizations")
# HyPrColoc reports colocalizations with PP>0.25, which is what is tallied here.
#Locus                                 Sig GWAS            Pairwise GWAS colocalizations
#1p36.1                                   BMI                             0
#1p13.3                                   HDL                             0
#1q24_loc1                                BMI                             0
#1q24_loc2                                WHR                             4
#1q25.2                                   BMI                             0
#1q41                                     WHR                             4
#2p24                           Triglycerides                             1
#2p23.3                                   BMI                             1
#2p21                                     T2D                             0
#2q24.3                                   WHR                             6
#3p25                                     HDL                             1
#3q27                                     T2D                             1
#4q22_loc1                                BMI                             0
#4q22_loc2                                BMI                             0
#5q11.2                           HDL and T2D                             4
#5q31.2                                   BMI                             0
#6p21_loc1                      Triglycerides                             1
#6p21_loc2                                HDL                             2
#6q22.3                                   WHR                             3
#7q32         BMI, T2D, HDL and Triglycerides                             6
#7q36                                     HDL                             0
#8q21.1                                   BMI                             0
#8q21.2                                   BMI                             0
#10p13                                    BMI                             0
#10q22.2                                  BMI                             0
#10q24.2                                  BMI                             0
#10q26.3                                  BMI                             0
#11p11.2                                  HDL                             1
#11q13                                    BMI                             1
#12p13.33                                 BMI                             0
#12p13.1                                  BMI                             0
#12q13.13                                 BMI                             0
#15q21.3                                  HDL                             0
#15q24.1                                  BMI                             1
#17q21.2                                  BMI                             0
#17q25.3                                  HDL                             1
#18q21.1                                  HDL                             0
#18q21.3                                  BMI                             1
#19p13.1                T2D and Triglycerides                             1
#19q13.3_loc1                             HDL                             3
#19q13.3_loc2                             BMI                             5
#20q13.32                                 T2D                             0
#22q13.3                                  T2D                             0

# There was no pairwise colocalization between BMI and any other GWAS at 17p13.2, which is not surprising given that none
# of the other GWAS have significant signals there. Therefore, no need to look at this locus further in this script.

write.table(posColoc,"Pairwise_HyPrColoc_between_all_GWAS_for_select_loci.txt",sep = "\t",row.names = F,quote = F)

### Mult-trait analyses

# 1q24_loc2, which is the 4th locus in the list
# Start with all traits included
locus4all=hyprcoloc(as.matrix(filt_betas[[4]]),as.matrix(filt_ses[[4]]),
                    trait.names=colnames(filt_betas[[4]]),snp.id=rownames(filt_betas[[4]]),ld.matrix = LD[[4]],
                    snpscores = T)
# There is just one cluster including WHR, T2D, HDL and TriG with a PP=0.7304 that is best explained by rs714515. BMI and all
# cis-QTLs were excluded.

# 1q41, which is the 6th locus in the list
# Start with all traits included
locus6all=hyprcoloc(as.matrix(filt_betas[[6]]),as.matrix(filt_ses[[6]]),
                    trait.names=colnames(filt_betas[[6]]),snp.id=rownames(filt_betas[[6]]),ld.matrix = LD[[6]],
                    snpscores = T)
# The GWAS are split into two clusters. The first is just WHR and T2D with a PP=0.9579 that is mostly explained by rs2820443.
# The other cluster includes BMI, HDL and TriG with a PP=0.6253 that is mostly explained by rs11118308. All cis-QTLs were excluded.

# 2p24, which is the 7th locus in the list
# Start with all traits included
locus7all=hyprcoloc(as.matrix(filt_betas[[7]]),as.matrix(filt_ses[[7]]),
                    trait.names=colnames(filt_betas[[7]]),snp.id=rownames(filt_betas[[7]]),ld.matrix = LD[[7]],
                    snpscores = T)
# Only HDL and TriG cluster together with PP=0.9978 that is almost entirely explained by rs676210. All other traits were dropped.

# 2p23.3, which is the 8th locus in the list
# Start with all traits included
locus8all=hyprcoloc(as.matrix(filt_betas[[8]]),as.matrix(filt_ses[[8]]),
                    trait.names=colnames(filt_betas[[8]]),snp.id=rownames(filt_betas[[8]]),ld.matrix = LD[[8]],
                    snpscores = T)
# Two clusters were found (1 and 2). The first includes 5 cis-eQTLs (but no GWAS) with PP=0.4660 best explained by rs13020526.
# The other cluster was just the pairwise colocalization of BMI and TriG with PP=0.4411 best explained by rs3739081.

# 2q24.3, which is the 10th locus in the list
# Start with all traits included
locus10all=hyprcoloc(as.matrix(filt_betas[[10]]),as.matrix(filt_ses[[10]]),
                    trait.names=colnames(filt_betas[[10]]),snp.id=rownames(filt_betas[[10]]),ld.matrix = LD[[10]],
                    snpscores = T)
# Only one cluster including BMI, WHR, T2D and TriG with PP=0.846 that is mostly explained by rs1128249.

# 3p25, which is the 11th locus in the list
# Start with all traits included
locus11all=hyprcoloc(as.matrix(filt_betas[[11]]),as.matrix(filt_ses[[11]]),
                     trait.names=colnames(filt_betas[[11]]),snp.id=rownames(filt_betas[[11]]),ld.matrix = LD[[11]],
                     snpscores = T)
# Only one cluster that is the pairwise colocalization of HDL and TriG with PP=0.9705 that is almost entirely explained by rs2606736.

# 3q27, which is the 12th locus in the list
# Start with all traits included
locus12all=hyprcoloc(as.matrix(filt_betas[[12]]),as.matrix(filt_ses[[12]]),
                    trait.names=colnames(filt_betas[[12]]),snp.id=rownames(filt_betas[[12]]),ld.matrix = LD[[12]],
                    snpscores = T)
# There is only one cluster with WHR, T2D and cis-e_IGF2BP2 with a PP=0.4895 that is best explained by rs7646518.

# 5q11.2, which is the 15th locus in the list
# Start with all traits included
locus15all=hyprcoloc(as.matrix(filt_betas[[15]]),as.matrix(filt_ses[[15]]),
                    trait.names=colnames(filt_betas[[15]]),snp.id=rownames(filt_betas[[15]]),ld.matrix = LD[[15]],
                    snpscores = T)
# There is only one cluster with BMI, HDL, TriG and cis-a_FST with PP=0.3075 that is best explained by rs6450176.

# The cis-a_FST colocalizations do not look good and are probably messing up the cluster. Let's just try the GWAS.
locus15gwas=hyprcoloc(as.matrix(filt_betas[[15]]),as.matrix(filt_ses[[15]]),
                    trait.names=colnames(filt_betas[[15]]),snp.id=rownames(filt_betas[[15]]),ld.matrix = LD[[15]],
                    trait.subset = c("BMI","WHR","T2D","HDL","TriG"),snpscores = T)
# There is a cluster with BMI, HDL and TriG with PP=0.5257 that is best explained by rs6450176. This is better than when cis-a_FST
# was included. T2D still does not cluster despite pairwise colocalization with BMI with PP=0.9931, but that was best explained
# by rs4865796. Based on the pairwise colocalizations between the GWAS at this locus, it seems more likely that BMI and T2D are
# explained by the same variant (rs4865796) while HDL and TriG are best explained by a different variant (rs6450176).

# 6p21_loc1, which is the 17th locus in the list
# Start with all traits included
locus17all=hyprcoloc(as.matrix(filt_betas[[17]]),as.matrix(filt_ses[[17]]),
                     trait.names=colnames(filt_betas[[17]]),snp.id=rownames(filt_betas[[17]]),ld.matrix = LD[[17]],
                     snpscores = T)
# Five clusters were found (1,2,3,4,7), but only cluster 2 included any GWAS (WHR, HDL and TriG) along with 4 cis-eQTLs with 
# PP=0.3918 that is best explained by rs3130557.

# Let's check whether the GWAS cluster better without the cis-QTLs.
locus17gwas=hyprcoloc(as.matrix(filt_betas[[17]]),as.matrix(filt_ses[[17]]),
                     trait.names=colnames(filt_betas[[17]]),snp.id=rownames(filt_betas[[17]]),ld.matrix = LD[[17]],
                     trait.subset = c("BMI","WHR","T2D","HDL","TriG"),snpscores = T)
# Interestingly, the cis-QTLs are necessary for the GWAS to cluster at all. Although it I'm not sure why this doesn't at least
# recapitulate the pairwise colocalization already observed between WHR and TriG (PP=0.6358).

# 6p21_loc2, which is the 18th locus in the list
# Start with all traits included
locus18all=hyprcoloc(as.matrix(filt_betas[[18]]),as.matrix(filt_ses[[18]]),
                     trait.names=colnames(filt_betas[[18]]),snp.id=rownames(filt_betas[[18]]),ld.matrix = LD[[18]],
                     snpscores = T)
# There is only one cluster with BMI, HDL, cis-e_SNRPC and cis-e_UHRF1BP1 with PP=0.6441 that is entirely explained by rs2744957.

# 6q22.3, which is the 19th locus in the list
# Start with all traits included
locus19all=hyprcoloc(as.matrix(filt_betas[[19]]),as.matrix(filt_ses[[19]]),
                     trait.names=colnames(filt_betas[[19]]),snp.id=rownames(filt_betas[[19]]),ld.matrix = LD[[19]],
                     snpscores = T)
# There is only one cluster with WHR, HDL, TriG and cis-e_RSPO3 with PP=0.8424 that is best explained by rs1936805.

# 7q32, which is the 20th locus in the list
# Start with all traits included
locus20all=hyprcoloc(as.matrix(filt_betas[[20]]),as.matrix(filt_ses[[20]]),
                    trait.names=colnames(filt_betas[[20]]),snp.id=rownames(filt_betas[[20]]),ld.matrix = LD[[20]],
                    snpscores = T)
# There are 3 clusters (1,2,3). Cluster 2 and 3 were seen in other analyses, but seem irrelevant to any of our GWAS signals.
# However, cluster 1 includes BMI, T2D, HDL, TriG, cis-e_LINC-PINT, cis-e_KLF14 and cis-a_KLF14 with PP=0.8304 (WOW!) that is 
# best explained by rs4731702.

# 11p11.2, which is the 28th locus in the list
# Start with all traits included
locus28all=hyprcoloc(as.matrix(filt_betas[[28]]),as.matrix(filt_ses[[28]]),
                    trait.names=colnames(filt_betas[[28]]),snp.id=rownames(filt_betas[[28]]),ld.matrix = LD[[28]],
                    snpscores = T)
# There are 3 clusters (1,2,7). The first cluster includes just BMI and cis-a_MADD with a PP=0.5162 that is entirely explained
# by rs7124681 that has a BMI P=3.2e-58, though this locus was selected because an HDL signal at the locus overlapped a
# cis-aQTL for PHF21A that was stronger than its eQTL. The rs7124681-MADD aQTL P=3.9e-4. Apparently I didn't notice this locus
# for BMI because MADD has a stronger eQTL overlapping the BMI locus, but HyPrColoc analyses have shown no evidence of the MADD
# eQTL colocalizing, only the aQTL. The second cluster was the pairwise colocalization of HDL and TriG GWAS signals with PP=0.9612 
# that is best explained by rs10501321 with an HDL P=3.54e-38 (which is extremely similar to the top HDL GWAS SNP at the locus) 
# and a TriG P=1.41e-8 (which is the top TriG GWAS SNP at the locus). Cluster 7 included cis-e_ARHGAP1, cis-e_LRP4, cis-e_FNBP4 
# and cis-e_PACSIN3 with PP=0.5704.

# 11q13, which is the 29th locus in the list
# Start with all traits included
locus29all=hyprcoloc(as.matrix(filt_betas[[29]]),as.matrix(filt_ses[[29]]),
                     trait.names=colnames(filt_betas[[29]]),snp.id=rownames(filt_betas[[29]]),ld.matrix = LD[[29]],
                     snpscores = T)
# Two cluster were found (1,2). The first includes WHR, T2D, HDL, TriG and cis-e_VEGFB with a PP=0.6219 that is entirely explained
# by rs3751120. The rs3751120-VEGFB eQTL P=1.24e-15. This locus was of interest due to overlap between BMI GWAS signal and 3 FDR sig
# BMI MR trans-aQTLs, but none pairwise colocalized in the BMI analysis. The BMI signal did not cluster in this analysis. The other
# cluster only contained cis-e_PRDX5 and cis-e_TRMT112 with a PP=0.9761 that is best explained by rs9787810.

# 15q24.1, which is the 34th locus in the list
# Start with all traits included
locus34all=hyprcoloc(as.matrix(filt_betas[[34]]),as.matrix(filt_ses[[34]]),
                     trait.names=colnames(filt_betas[[34]]),snp.id=rownames(filt_betas[[34]]),ld.matrix = LD[[34]],
                     snpscores = T)
# Only 1 cluster was found, which only included BMI and TriG with a PP=0.9632 that is best explained by rs8030477 that has a BMI
# P=4.4e-25 (which is very similar to the top BMI GWAS SNP at the locus) and a TriG P=5.75e-5.

# 17q25.3, which is the 36th locus in the list
# Start with all traits included
locus36all=hyprcoloc(as.matrix(filt_betas[[36]]),as.matrix(filt_ses[[36]]),
                     trait.names=colnames(filt_betas[[36]]),snp.id=rownames(filt_betas[[36]]),ld.matrix = LD[[36]],
                     snpscores = T)
# Three clusters were found (1,2,3). Cluster 1 included BMI, T2D, cis-e_USP36, cis-e_TIMP2 and cis-e_CYTH1 with a PP=0.4352 that
# is best explained by rs1044486. Cluster 2 included HDL, TriG and cis-e_PGS1 with a PP=0.9160 (WOW!) that is entirely explained
# by rs4969178 that has an HDL P=1.53e-12 and a TriG P=5.70e-6. Cluster 3 includes cis-e_AFMID, cis-e_SYNGR2, cis-e_LGALS3BP and
# cis-a_LGALS3BP with a PP=0.4650 that is almost entirely explained by rs12602116. This locus was of interest due to an FDR sig
# HDL MR trans-eQTL that overlapped the HDL GWAS signal that turned out to pairwise colocalize with HDL.

# 18q21.3, which is the 38th locus in the list
# Start with all traits included
locus38all=hyprcoloc(as.matrix(filt_betas[[38]]),as.matrix(filt_ses[[38]]),
                     trait.names=colnames(filt_betas[[38]]),snp.id=rownames(filt_betas[[38]]),ld.matrix = LD[[38]],
                     snpscores = T)
# There was only one cluster that included only BMI and HDL with PP=0.9784 that is entirely explained by rs6567160 that has a BMI
# P=1.8e-178 (that is in perfect LD with the top BMI SNP) and an HDL P=2.92e-9 (which is the top HDL SNP). This locus was of
# interest due to an FDR sig BMI MR trans-aQTL that overlapped the BMI GWAS signal. However, the BMI signal did not pairwise with 
# any QTLs at this locus.

# 19p13.1, which is the 39th locus in the list
# Start with all traits included
locus39all=hyprcoloc(as.matrix(filt_betas[[39]]),as.matrix(filt_ses[[39]]),
                     trait.names=colnames(filt_betas[[39]]),snp.id=rownames(filt_betas[[39]]),ld.matrix = LD[[39]],
                     snpscores = T)
# There was only 1 cluster that included only T2D and TriG with a PP=1 that is entirely explained by rs10401969 that is the top
# GWAS SNP for both T2D and TriG. However, this locus was of interest due to an overlap of these GWAS signals with ZNF101 cis-aQTL
# signal that was stronger than its eQTL, but it does not pairwise colocalize with either GWAS signal.

# 19q13.3_loc1, which is the 40th locus in the list
# Start with all traits included
locus40all=hyprcoloc(as.matrix(filt_betas[[40]]),as.matrix(filt_ses[[40]]),
                     trait.names=colnames(filt_betas[[40]]),snp.id=rownames(filt_betas[[40]]),ld.matrix = LD[[40]],
                     snpscores = T)
# Three clusters were found (1,2,6). The first cluster only included TriG and cis-e_APOE with a PP=0.9969 that is entirely explained
# by rs439401 that has a TriG P=1.42e-66. This is a striking colocalization for a very strong TriG signal, but it was not on my
# radar because the locus doesn't have any FDR sig cis-aQTL or TriG MR trans-QTL. The second cluster includes BMI, WHR and HDL with
# a PP=0.5668 that is entirely explained by rs2075650 that has a BMI P=1.25e-25, a WHR P=1.3e-4 and an HDL P=9.72e-26. Note that the
# pairwise colocalization between BMI and HDL was much stronger (PP=0.9975). This locus was of interest because the HDL signal
# overlapped ZNF221 and ZNF404 cis-aQTLs that were stronger than their eQTLs, but they did not pairwise colocalize. Cluster 6 was
# just 3 cis-eQTLs.

# 19q13.3_loc2, which is the 41st locus in the list
# Start with all traits included
locus41all=hyprcoloc(as.matrix(filt_betas[[41]]),as.matrix(filt_ses[[41]]),
                     trait.names=colnames(filt_betas[[41]]),snp.id=rownames(filt_betas[[41]]),ld.matrix = LD[[41]],
                     snpscores = T)
# Two clusters were found (1,3). The first included T2D, HDL and TriG with a PP=0.7279 that is best explained by rs2303108 with a
# T2D P=9.5e-5, an HDL P=6.32e-7 and a TriG P=9.81e-7. Note that this locus was of interest due to the overlap between a BMI signal
# and a BMI MR trans-eQTL, which did not pairwise colocalize in the BMI analysis. Also, although the BMI signal pairwise colocalizes
# with T2D, HDL and TriG, it is apparently different enough to be excluded from the mutli-trait cluster. The other cluster is
# simply 2 cis-eQTLs.

# Write results to text files
write.table(posColoc,"Pairwise_HyPrColoc_between_all_GWAS_and_cis-QTLs_for_select_loci.txt",sep = "\t",row.names = F,quote = F)
write.table(locus4all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_1q24_loc2.txt",sep = "\t",row.names = F,quote = F)
write.table(locus6all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_1q41.txt",sep = "\t",row.names = F,quote = F)
write.table(locus7all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_2p24.txt",sep = "\t",row.names = F,quote = F)
write.table(locus8all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_2p23.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus10all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_2q24.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus11all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_3p25.txt",sep = "\t",row.names = F,quote = F)
write.table(locus12all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_3q27.txt",sep = "\t",row.names = F,quote = F)
write.table(locus15all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_5q11.2.txt",sep = "\t",row.names = F,quote = F)
write.table(locus17all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_6p21_loc1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus18all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_6p21_loc2.txt",sep = "\t",row.names = F,quote = F)
write.table(locus19all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_6q22.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus20all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_7q32.txt",sep = "\t",row.names = F,quote = F)
write.table(locus28all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_11p11.2.txt",sep = "\t",row.names = F,quote = F)
write.table(locus29all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_11q13.txt",sep = "\t",row.names = F,quote = F)
write.table(locus34all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_15q24.1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus36all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_17q25.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus38all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_18q21.3.txt",sep = "\t",row.names = F,quote = F)
write.table(locus39all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_19p13.1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus40all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_19q13.3_loc1.txt",sep = "\t",row.names = F,quote = F)
write.table(locus41all[[1]],"HyPrColoc_of_all_GWAS_and_cis-QTLs_for_19q13.3_loc2.txt",sep = "\t",row.names = F,quote = F)

### Let's make some targetted LocusCompare plots

# 1q24
# First, grab the necessary P-values for the SNPs used in the HyPrColoc analyses for the traits of interest
bmi2=filt_bmi[na.omit(match(rownames(filt_betas[[2]]),filt_bmi$SNP)),c(3,9)]
whr2=filt_whr[na.omit(match(rownames(filt_betas[[2]]),filt_whr$SNP)),c(1,7)]
t2d2=filt_t2d[na.omit(match(rownames(filt_betas[[2]]),filt_t2d$rsID)),c(1,9)]
hdl2=filt_hdl[na.omit(match(rownames(filt_betas[[2]]),filt_hdl$SNP)),c(3,9)]
triG2=filt_triG[na.omit(match(rownames(filt_betas[[2]]),filt_triG$SNP)),c(3,9)]
ePRRC2C=filt_cis_eqtl[filt_cis_eqtl$gene=="PRRC2C",c(1,4)]
ePRRC2C=ePRRC2C[na.omit(match(rownames(filt_betas[[2]]),ePRRC2C$snps)),]
colnames(bmi2)=c("rsid","pval")
colnames(whr2)=c("rsid","pval")
colnames(t2d2)=c("rsid","pval")
colnames(hdl2)=c("rsid","pval")
colnames(triG2)=c("rsid","pval")
colnames(ePRRC2C)=c("rsid","pval")
rownames(bmi2)=bmi2$rsid
rownames(whr2)=whr2$rsid
rownames(t2d2)=t2d2$rsid
rownames(hdl2)=hdl2$rsid
rownames(triG2)=triG2$rsid
rownames(ePRRC2C)=ePRRC2C$rsid

# Make the LocusCompare plots
locuscompare(in_fn1=bmi2,in_fn2=whr2,title1 = "BMI GWAS", title2 = "WHR (BMI adj) GWAS",snp = "rs714515")
locuscompare(in_fn1=whr2,in_fn2=t2d2,title1 = "WHR (BMI adj) GWAS", title2 = "T2D (BMI adj) GWAS",snp = "rs714515")
locuscompare(in_fn1=whr2,in_fn2=hdl2,title1 = "WHR (BMI adj) GWAS", title2 = "HDL GWAS",snp = "rs714515")
locuscompare(in_fn1=whr2,in_fn2=triG2,title1 = "WHR (BMI adj) GWAS", title2 = "Triglycerides GWAS",snp = "rs714515")
locuscompare(in_fn1=bmi2,in_fn2=ePRRC2C,title1 = "BMI GWAS", title2 = "PRRC2C eQTL",snp = "rs16864515")

# 3q27
# First, grab the necessary P-values for the SNPs used in the HyPrColoc analyses for the traits of interest
bmi3=filt_bmi[na.omit(match(rownames(filt_betas[[3]]),filt_bmi$SNP)),c(3,9)]
whr3=filt_whr[na.omit(match(rownames(filt_betas[[3]]),filt_whr$SNP)),c(1,7)]
t2d3=filt_t2d[na.omit(match(rownames(filt_betas[[3]]),filt_t2d$rsID)),c(1,9)]
hdl3=filt_hdl[na.omit(match(rownames(filt_betas[[3]]),filt_hdl$SNP)),c(3,9)]
triG3=filt_triG[na.omit(match(rownames(filt_betas[[3]]),filt_triG$SNP)),c(3,9)]
eIGF2BP2=filt_cis_eqtl[filt_cis_eqtl$gene=="IGF2BP2",c(1,4)]
eIGF2BP2=eIGF2BP2[na.omit(match(rownames(filt_betas[[3]]),eIGF2BP2$snps)),]
colnames(bmi3)=c("rsid","pval")
colnames(whr3)=c("rsid","pval")
colnames(t2d3)=c("rsid","pval")
colnames(hdl3)=c("rsid","pval")
colnames(triG3)=c("rsid","pval")
colnames(eIGF2BP2)=c("rsid","pval")
rownames(bmi3)=bmi3$rsid
rownames(whr3)=whr3$rsid
rownames(t2d3)=t2d3$rsid
rownames(hdl3)=hdl3$rsid
rownames(triG3)=triG3$rsid
rownames(eIGF2BP2)=eIGF2BP2$rsid

# Make the LocusCompare plots
locuscompare(in_fn1=bmi3,in_fn2=whr3,title1 = "BMI GWAS", title2 = "WHR (BMI adj) GWAS",snp = "rs7646518")
locuscompare(in_fn1=whr3,in_fn2=t2d3,title1 = "WHR GWAS", title2 = "T2D (BMI adj) GWAS",snp = "rs7646518")
locuscompare(in_fn1=t2d3,in_fn2=hdl3,title1 = "T2D (BMI adj) GWAS", title2 = "HDL GWAS",snp = "rs7646518")
locuscompare(in_fn1=t2d3,in_fn2=triG3,title1 = "T2D (BMI adj) GWAS", title2 = "Triglycerides GWAS",snp = "rs7646518")
locuscompare(in_fn1=t2d3,in_fn2=eIGF2BP2,title1 = "T2D (BMI adj) GWAS", title2 = "IGF2BP2 eQTL",snp = "rs7646518")

# 5q11
# First, grab the necessary P-values for the SNPs used in the HyPrColoc analyses for the traits of interest
bmi4=filt_bmi[na.omit(match(rownames(filt_betas[[4]]),filt_bmi$SNP)),c(3,9)]
whr4=filt_whr[na.omit(match(rownames(filt_betas[[4]]),filt_whr$SNP)),c(1,7)]
t2d4=filt_t2d[na.omit(match(rownames(filt_betas[[4]]),filt_t2d$rsID)),c(1,9)]
hdl4=filt_hdl[na.omit(match(rownames(filt_betas[[4]]),filt_hdl$SNP)),c(3,9)]
triG4=filt_triG[na.omit(match(rownames(filt_betas[[4]]),filt_triG$SNP)),c(3,9)]
aFST=filt_cis_aqtl[filt_cis_aqtl$gene=="FST",c(1,4)]
aFST=aFST[na.omit(match(rownames(filt_betas[[4]]),aFST$snps)),]
colnames(bmi4)=c("rsid","pval")
colnames(whr4)=c("rsid","pval")
colnames(t2d4)=c("rsid","pval")
colnames(hdl4)=c("rsid","pval")
colnames(triG4)=c("rsid","pval")
colnames(aFST)=c("rsid","pval")
rownames(bmi4)=bmi4$rsid
rownames(whr4)=whr4$rsid
rownames(t2d4)=t2d4$rsid
rownames(hdl4)=hdl4$rsid
rownames(triG4)=triG4$rsid
rownames(aFST)=aFST$rsid

# Make the LocusCompare plots
locuscompare(in_fn1=bmi4,in_fn2=whr4,title1 = "BMI GWAS", title2 = "WHR (BMI adj) GWAS",snp = "rs6450176")
locuscompare(in_fn1=bmi4,in_fn2=t2d4,title1 = "BMI GWAS", title2 = "T2D (BMI adj) GWAS",snp = "rs6450176")
locuscompare(in_fn1=bmi4,in_fn2=hdl4,title1 = "BMI GWAS", title2 = "HDL GWAS",snp = "rs6450176")
locuscompare(in_fn1=bmi4,in_fn2=triG4,title1 = "BMI GWAS", title2 = "Triglycerides GWAS",snp = "rs6450176")
locuscompare(in_fn1=bmi4,in_fn2=aFST,title1 = "BMI GWAS", title2 = "FST cis-aQTL",snp = "rs6450176")

# 6p21
# First, grab the necessary P-values for the SNPs used in the HyPrColoc analyses for the traits of interest
bmi5=filt_bmi[na.omit(match(rownames(filt_betas[[5]]),filt_bmi$SNP)),c(3,9)]
whr5=filt_whr[na.omit(match(rownames(filt_betas[[5]]),filt_whr$SNP)),c(1,7)]
t2d5=filt_t2d[na.omit(match(rownames(filt_betas[[5]]),filt_t2d$rsID)),c(1,9)]
hdl5=filt_hdl[na.omit(match(rownames(filt_betas[[5]]),filt_hdl$SNP)),c(3,9)]
triG5=filt_triG[na.omit(match(rownames(filt_betas[[5]]),filt_triG$SNP)),c(3,9)]
eSNRPC=filt_cis_eqtl[filt_cis_eqtl$gene=="SNRPC",c(1,4)]
eSNRPC=eSNRPC[na.omit(match(rownames(filt_betas[[5]]),eSNRPC$snps)),]
eUHRF1BP1=filt_cis_eqtl[filt_cis_eqtl$gene=="UHRF1BP1",c(1,4)]
eUHRF1BP1=eUHRF1BP1[na.omit(match(rownames(filt_betas[[5]]),eUHRF1BP1$snps)),]
colnames(bmi5)=c("rsid","pval")
colnames(whr5)=c("rsid","pval")
colnames(t2d5)=c("rsid","pval")
colnames(hdl5)=c("rsid","pval")
colnames(triG5)=c("rsid","pval")
colnames(eSNRPC)=c("rsid","pval")
colnames(eUHRF1BP1)=c("rsid","pval")
rownames(bmi5)=bmi5$rsid
rownames(whr5)=whr5$rsid
rownames(t2d5)=t2d5$rsid
rownames(hdl5)=hdl5$rsid
rownames(triG5)=triG5$rsid
rownames(eSNRPC)=eSNRPC$rsid
rownames(eUHRF1BP1)=eUHRF1BP1$rsid

# Make the LocusCompare plots
locuscompare(in_fn1=bmi5,in_fn2=whr5,title1 = "BMI GWAS", title2 = "WHR (BMI adj) GWAS",snp = "rs2744957")
locuscompare(in_fn1=bmi5,in_fn2=t2d5,title1 = "BMI GWAS", title2 = "T2D (BMI adj) GWAS",snp = "rs2744957")
locuscompare(in_fn1=bmi5,in_fn2=hdl5,title1 = "BMI GWAS", title2 = "HDL GWAS",snp = "rs2744957")
locuscompare(in_fn1=bmi5,in_fn2=triG5,title1 = "BMI GWAS", title2 = "Triglycerides GWAS",snp = "rs2744957")
locuscompare(in_fn1=bmi5,in_fn2=eSNRPC,title1 = "BMI GWAS", title2 = "SNRPC cis-eQTL",snp = "rs2744957")
locuscompare(in_fn1=bmi5,in_fn2=eUHRF1BP1,title1 = "BMI GWAS", title2 = "UHRF1BP1",snp = "rs2744957")

# 7q32
# First, grab the necessary P-values for the SNPs used in the HyPrColoc analyses for the traits of interest
bmi6=filt_bmi[na.omit(match(rownames(filt_betas[[6]]),filt_bmi$SNP)),c(3,9)]
whr6=filt_whr[na.omit(match(rownames(filt_betas[[6]]),filt_whr$SNP)),c(1,7)]
t2d6=filt_t2d[na.omit(match(rownames(filt_betas[[6]]),filt_t2d$rsID)),c(1,9)]
hdl6=filt_hdl[na.omit(match(rownames(filt_betas[[6]]),filt_hdl$SNP)),c(3,9)]
triG6=filt_triG[na.omit(match(rownames(filt_betas[[6]]),filt_triG$SNP)),c(3,9)]
eLINC=filt_cis_eqtl[filt_cis_eqtl$gene=="LINC-PINT",c(1,4)]
eLINC=eLINC[na.omit(match(rownames(filt_betas[[6]]),eLINC$snps)),]
colnames(bmi6)=c("rsid","pval")
colnames(whr6)=c("rsid","pval")
colnames(t2d6)=c("rsid","pval")
colnames(hdl6)=c("rsid","pval")
colnames(triG6)=c("rsid","pval")
colnames(eLINC)=c("rsid","pval")
rownames(bmi6)=bmi6$rsid
rownames(whr6)=whr6$rsid
rownames(t2d6)=t2d6$rsid
rownames(hdl6)=hdl6$rsid
rownames(triG6)=triG6$rsid
rownames(eLINC)=eLINC$rsid

# Make the LocusCompare plots
locuscompare(in_fn1=bmi6,in_fn2=whr6,title1 = "BMI GWAS", title2 = "WHR (BMI adj) GWAS",snp = "rs4731702")
locuscompare(in_fn1=bmi6,in_fn2=t2d6,title1 = "BMI GWAS", title2 = "T2D (BMI adj) GWAS",snp = "rs4731702")
locuscompare(in_fn1=bmi6,in_fn2=hdl6,title1 = "BMI GWAS", title2 = "HDL GWAS",snp = "rs4731702")
locuscompare(in_fn1=bmi6,in_fn2=triG6,title1 = "BMI GWAS", title2 = "Triglycerides GWAS",snp = "rs4731702")
locuscompare(in_fn1=bmi6,in_fn2=eLINC,title1 = "BMI GWAS", title2 = "LINC-PINT cis-eQTL",snp = "rs4731702")

# 11q13
# First, grab the necessary P-values for the SNPs used in the HyPrColoc analyses for the traits of interest
bmi7=filt_bmi[na.omit(match(rownames(filt_betas[[7]]),filt_bmi$SNP)),c(3,9)]
whr7=filt_whr[na.omit(match(rownames(filt_betas[[7]]),filt_whr$SNP)),c(1,7)]
t2d7=filt_t2d[na.omit(match(rownames(filt_betas[[7]]),filt_t2d$rsID)),c(1,9)]
hdl7=filt_hdl[na.omit(match(rownames(filt_betas[[7]]),filt_hdl$SNP)),c(3,9)]
triG7=filt_triG[na.omit(match(rownames(filt_betas[[7]]),filt_triG$SNP)),c(3,9)]
eVEGFB=filt_cis_eqtl[filt_cis_eqtl$gene=="VEGFB",c(1,4)]
eVEGFB=eVEGFB[na.omit(match(rownames(filt_betas[[7]]),eVEGFB$snps)),]
eTEX40=filt_cis_eqtl[filt_cis_eqtl$gene=="TEX40",c(1,4)]
eTEX40=eTEX40[na.omit(match(rownames(filt_betas[[7]]),eTEX40$snps)),]
colnames(bmi7)=c("rsid","pval")
colnames(whr7)=c("rsid","pval")
colnames(t2d7)=c("rsid","pval")
colnames(hdl7)=c("rsid","pval")
colnames(triG7)=c("rsid","pval")
colnames(eVEGFB)=c("rsid","pval")
colnames(eTEX40)=c("rsid","pval")
rownames(bmi7)=bmi7$rsid
rownames(whr7)=whr7$rsid
rownames(t2d7)=t2d7$rsid
rownames(hdl7)=hdl7$rsid
rownames(triG7)=triG7$rsid
rownames(eVEGFB)=eVEGFB$rsid
rownames(eTEX40)=eTEX40$rsid

# Make the LocusCompare plots
locuscompare(in_fn1=bmi7,in_fn2=whr7,title1 = "BMI GWAS", title2 = "WHR (BMI adj) GWAS",snp = "rs3751120")
locuscompare(in_fn1=whr7,in_fn2=t2d7,title1 = "WHR (BMI adj) GWAS", title2 = "T2D (BMI adj) GWAS",snp = "rs3751120")
locuscompare(in_fn1=whr7,in_fn2=hdl7,title1 = "WHR (BMI adj) GWAS", title2 = "HDL GWAS",snp = "rs3751120")
locuscompare(in_fn1=whr7,in_fn2=triG7,title1 = "WHR (BMI adj) GWAS", title2 = "Triglycerides GWAS",snp = "rs3751120")
locuscompare(in_fn1=whr7,in_fn2=eVEGFB,title1 = "WHR (BMI adj) GWAS", title2 = "VEGFB cis-eQTL",snp = "rs3751120")
locuscompare(in_fn1=bmi7,in_fn2=eTEX40,title1 = "BMI GWAS", title2 = "TEX40 cis-eQTL",snp = "rs477895")

# 7q32
# First, grab the necessary P-values for the SNPs used in the HyPrColoc analyses for the traits of interest
bmi6=filt_bmi[na.omit(match(rownames(filt_betas[[6]]),filt_bmi$SNP)),c(3,9)]
whr6=filt_whr[na.omit(match(rownames(filt_betas[[6]]),filt_whr$SNP)),c(1,7)]
t2d6=filt_t2d[na.omit(match(rownames(filt_betas[[6]]),filt_t2d$rsID)),c(1,9)]
hdl6=filt_hdl[na.omit(match(rownames(filt_betas[[6]]),filt_hdl$SNP)),c(3,9)]
triG6=filt_triG[na.omit(match(rownames(filt_betas[[6]]),filt_triG$SNP)),c(3,9)]
eLINC=filt_cis_eqtl[filt_cis_eqtl$gene=="LINC-PINT",c(1,4)]
eLINC=eLINC[na.omit(match(rownames(filt_betas[[6]]),eLINC$snps)),]
colnames(bmi6)=c("rsid","pval")
colnames(whr6)=c("rsid","pval")
colnames(t2d6)=c("rsid","pval")
colnames(hdl6)=c("rsid","pval")
colnames(triG6)=c("rsid","pval")
colnames(eLINC)=c("rsid","pval")
rownames(bmi6)=bmi6$rsid
rownames(whr6)=whr6$rsid
rownames(t2d6)=t2d6$rsid
rownames(hdl6)=hdl6$rsid
rownames(triG6)=triG6$rsid
rownames(eLINC)=eLINC$rsid

# Make the LocusCompare plots
locuscompare(in_fn1=bmi6,in_fn2=whr6,title1 = "BMI GWAS", title2 = "WHR (BMI adj) GWAS",snp = "rs4731702")
locuscompare(in_fn1=bmi6,in_fn2=t2d6,title1 = "BMI GWAS", title2 = "T2D (BMI adj) GWAS",snp = "rs4731702")
locuscompare(in_fn1=bmi6,in_fn2=hdl6,title1 = "BMI GWAS", title2 = "HDL GWAS",snp = "rs4731702")
locuscompare(in_fn1=bmi6,in_fn2=triG6,title1 = "BMI GWAS", title2 = "Triglycerides GWAS",snp = "rs4731702")

# 7q32
# First, grab the necessary P-values for the SNPs used in the HyPrColoc analyses for the traits of interest
bmi6=filt_bmi[na.omit(match(rownames(filt_betas[[6]]),filt_bmi$SNP)),c(3,9)]
whr6=filt_whr[na.omit(match(rownames(filt_betas[[6]]),filt_whr$SNP)),c(1,7)]
t2d6=filt_t2d[na.omit(match(rownames(filt_betas[[6]]),filt_t2d$rsID)),c(1,9)]
hdl6=filt_hdl[na.omit(match(rownames(filt_betas[[6]]),filt_hdl$SNP)),c(3,9)]
triG6=filt_triG[na.omit(match(rownames(filt_betas[[6]]),filt_triG$SNP)),c(3,9)]
eLINC=filt_cis_eqtl[filt_cis_eqtl$gene=="LINC-PINT",c(1,4)]
eLINC=eLINC[na.omit(match(rownames(filt_betas[[6]]),eLINC$snps)),]
colnames(bmi6)=c("rsid","pval")
colnames(whr6)=c("rsid","pval")
colnames(t2d6)=c("rsid","pval")
colnames(hdl6)=c("rsid","pval")
colnames(triG6)=c("rsid","pval")
colnames(eLINC)=c("rsid","pval")
rownames(bmi6)=bmi6$rsid
rownames(whr6)=whr6$rsid
rownames(t2d6)=t2d6$rsid
rownames(hdl6)=hdl6$rsid
rownames(triG6)=triG6$rsid
rownames(eLINC)=eLINC$rsid

# Make the LocusCompare plots
locuscompare(in_fn1=bmi6,in_fn2=whr6,title1 = "BMI GWAS", title2 = "WHR (BMI adj) GWAS",snp = "rs4731702")
locuscompare(in_fn1=bmi6,in_fn2=t2d6,title1 = "BMI GWAS", title2 = "T2D (BMI adj) GWAS",snp = "rs4731702")
locuscompare(in_fn1=bmi6,in_fn2=hdl6,title1 = "BMI GWAS", title2 = "HDL GWAS",snp = "rs4731702")
locuscompare(in_fn1=bmi6,in_fn2=triG6,title1 = "BMI GWAS", title2 = "Triglycerides GWAS",snp = "rs4731702")

# 7q32
# First, grab the necessary P-values for the SNPs used in the HyPrColoc analyses for the traits of interest
bmi6=filt_bmi[na.omit(match(rownames(filt_betas[[6]]),filt_bmi$SNP)),c(3,9)]
whr6=filt_whr[na.omit(match(rownames(filt_betas[[6]]),filt_whr$SNP)),c(1,7)]
t2d6=filt_t2d[na.omit(match(rownames(filt_betas[[6]]),filt_t2d$rsID)),c(1,9)]
hdl6=filt_hdl[na.omit(match(rownames(filt_betas[[6]]),filt_hdl$SNP)),c(3,9)]
triG6=filt_triG[na.omit(match(rownames(filt_betas[[6]]),filt_triG$SNP)),c(3,9)]
eLINC=filt_cis_eqtl[filt_cis_eqtl$gene=="LINC-PINT",c(1,4)]
eLINC=eLINC[na.omit(match(rownames(filt_betas[[6]]),eLINC$snps)),]
colnames(bmi6)=c("rsid","pval")
colnames(whr6)=c("rsid","pval")
colnames(t2d6)=c("rsid","pval")
colnames(hdl6)=c("rsid","pval")
colnames(triG6)=c("rsid","pval")
colnames(eLINC)=c("rsid","pval")
rownames(bmi6)=bmi6$rsid
rownames(whr6)=whr6$rsid
rownames(t2d6)=t2d6$rsid
rownames(hdl6)=hdl6$rsid
rownames(triG6)=triG6$rsid
rownames(eLINC)=eLINC$rsid

# Make the LocusCompare plots
locuscompare(in_fn1=bmi6,in_fn2=whr6,title1 = "BMI GWAS", title2 = "WHR (BMI adj) GWAS",snp = "rs4731702")
locuscompare(in_fn1=bmi6,in_fn2=t2d6,title1 = "BMI GWAS", title2 = "T2D (BMI adj) GWAS",snp = "rs4731702")
locuscompare(in_fn1=bmi6,in_fn2=hdl6,title1 = "BMI GWAS", title2 = "HDL GWAS",snp = "rs4731702")
locuscompare(in_fn1=bmi6,in_fn2=triG6,title1 = "BMI GWAS", title2 = "Triglycerides GWAS",snp = "rs4731702")

















library(ggplot2)

# Set working dir
setwd("YOUR WORKING DIRECTORY")

# Read in data. The file locations are relative to your working directory, so adjust accordingly.
eds=read.table("./Wang_2020_EDS_scores.txt",sep = "\t",header = T) # There are duplicate ENSEMBL GeneIDs, so can't use row.names yet
id_map=read.table("./Eurobats_adipose_expressed_genes_with_names_TPM_with_ENSEMBL_IDs.txt",sep = "\t",header = T)
reg_list=read.table("./Full_Califano_regulator_list_plus_LINC-PINT_by_symbol.txt",sep = "\t",header = T)
cis_eBMI=read.table("./Eurobats_adipose_BMI_sig_cis-eQTLs_from_INT_logTPM_4213_regulators.txt",
                    sep = "\t",header = T)
cis_eWHR=read.table("./Eurobats_adipose_WHRadjBMI_sig_cis-eQTLs_from_INT_logTPM_of_4213_regulators.txt",
                    sep = "\t",header = T)
cis_eT2D=read.table("./Eurobats_adipose_BMIadjT2D_sig_cis-eQTLs_from_INT_logTPM_of_4213_regulators.txt",
                    sep = "\t",header = T)
cis_eHDL=read.table("./Eurobats_adipose_HDL_sig_cis-eQTLs_from_INT_logTPM_of_4213_regulators.txt",
                    sep = "\t",header = T)
cis_eTriG=read.table("./Eurobats_adipose_TriG_sig_cis-eQTLs_from_INT_logTPM_of_4213_regulators.txt",
                     sep = "\t",header = T)
cis_aBMI=read.table("./Eurobats_adipose_BMI_sig_cis-aQTLs_from_unnormalized_activities_4213_regulators.txt",
                    sep = "\t",header = T)
cis_aWHR=read.table("./Eurobats_adipose_WHRadjBMI_sig_cis-aQTLs_from_unnormalized_activities_of_4213_regulators.txt",
                    sep = "\t",header = T)
cis_aT2D=read.table("./Eurobats_adipose_BMIadjT2D_sig_cis-aQTLs_from_unnormalized_activities_of_4213_regulators.txt",
                    sep = "\t",header = T)
cis_aHDL=read.table("./Eurobats_adipose_HDL_sig_cis-aQTLs_from_unnormalized_activities_of_4213_regulators.txt",
                    sep = "\t",header = T)
cis_aTriG=read.table("./Eurobats_adipose_TriG_sig_cis-aQTLs_from_unnormalized_activities_of_4213_regulators.txt",
                     sep = "\t",header = T)
trans_eBMI=read.table("./Eurobats_adipose_BMI_sig_trans-eQTLs_from_INT_logTPM_of_BMI_RF_MRs.txt",
                      sep = "\t",header = T)
trans_eWHR=read.table("./Eurobats_adipose_WHRadjBMI_sig_trans-eQTLs_from_INT_logTPM_of_WHR_RF_MRs.txt",
                      sep = "\t",header = T)
trans_eT2D=read.table("./Eurobats_adipose_BMIadjT2D_sig_trans-eQTLs_from_INT_logTPM_of_HOMA-IR_RF_MRs.txt",
                      sep = "\t",header = T)
trans_eHDL=read.table("./Eurobats_adipose_HDL_sig_trans-eQTLs_from_INT_logTPM_of_HDL_RF_MRs.txt",
                      sep = "\t",header = T)
trans_eTriG=read.table("./Eurobats_adipose_TriG_sig_trans-eQTLs_from_INT_logTPM_of_TriG_RF_MRs.txt",
                       sep = "\t",header = T)
trans_aBMI=read.table("./Eurobats_adipose_BMI_sig_trans-aQTLs_from_unnormalized_activities_of_BMI_RF_MRs.txt",
                      sep = "\t",header = T)
trans_aWHR=read.table("./Eurobats_adipose_WHRadjBMI_sig_trans-aQTLs_from_unnormalized_activities_of_WHR_RF_MRs.txt",
                      sep = "\t",header = T)
trans_aT2D=read.table("./Eurobats_adipose_BMIadjT2D_sig_trans-aQTLs_from_unnormalized_activities_of_HOMA-IR_RF_MRs.txt",
                      sep = "\t",header = T)
trans_aHDL=read.table("./Eurobats_adipose_HDL_sig_trans-aQTLs_from_unnormalized_activities_of_HDL_RF_MRs.txt",
                      sep = "\t",header = T)
trans_aTriG=read.table("./Eurobats_adipose_TriG_sig_trans-aQTLs_from_unnormalized_activities_of_TriG_RF_MRs.txt",
                       sep = "\t",header = T)
all_snps_cis_eQTLs=read.table("./Eurobats_adipose_all_SNPs_cis-eQTLs_from_INT_logTPM_4213_regulators.txt",sep="\t",header=T)
all_snps_cis_aQTLs=read.table("./Eurobats_adipose_all_SNPs_cis-aQTLs_from_unnormalized_activities_4213_regulators.txt",sep="\t",header=T)
all_snps_trans_eQTLs=read.table("./Eurobats_adipose_all_SNPs_trans-eQTLs_from_INT_logTPM_of_4213_regulators.txt",sep="\t",header=T)
all_snps_trans_aQTLs=read.table("./Eurobats_adipose_all_SNPs_trans-aQTLs_from_unnormalized_activities_of_4213_regulators.txt",sep="\t",header=T)
full_cis_e=read.table("./Eurobats_adipose_all_SNPs_cis-eQTLs_from_INT_logTPM_of_all_expressed_genes.txt",
                      sep = "\t",header = T)
full_trans_e=read.table("./Eurobats_adipose_all_SNPs_trans-eQTLs_from_INT_logTPM_of_all_expressed_genes.txt",
                        sep = "\t",header = T)
# To be clear, the cis_eBMI, cis_aBMI, trans_eBMI, trans_aBMI, etc., are the original main QTL results for the regulators
# and MRs versus the indicated GWAS variants. The all_snps_cis-eQTLs, etc. are the QTLs for all regulators versus all SNPs.
# Finally, full_cis_e and full_trans_e are all variants versus all expressed genes, so there is no aQTL analysis for this set.

# First I need to map the ENSEMBL IDs to gene names, and filter and format the data
eds=eds[!(duplicated(eds$GeneSymbol)),]
eds=eds[na.omit(match(rownames(id_map),eds$GeneSymbol)),]
rownames(eds)=id_map[eds$GeneSymbol,1]

# In a preliminary analysis I found that the distributions for eGenes/aGenes across EDS bins depends on the P threshold
# applied. Consequently, I will try testing 3 different thresholds. Then, for each threshold, I will take the union of 
# the sets such that there is a single set each for cis-eGenes, cis-aGenes, trans-eGenes and trans-aGenes.
# P<0.05
cis_eBMI05=cis_eBMI[cis_eBMI$pvalue<=0.05,]
cis_eWHR05=cis_eWHR[cis_eWHR$pvalue<=0.05,]
cis_eT2D05=cis_eT2D[cis_eT2D$pvalue<=0.05,]
cis_eHDL05=cis_eHDL[cis_eHDL$pvalue<=0.05,]
cis_eTriG05=cis_eTriG[cis_eTriG$pvalue<=0.05,]
cis_aBMI05=cis_aBMI[cis_aBMI$pvalue<=0.05,]
cis_aWHR05=cis_aWHR[cis_aWHR$pvalue<=0.05,]
cis_aT2D05=cis_aT2D[cis_aT2D$pvalue<=0.05,]
cis_aHDL05=cis_aHDL[cis_aHDL$pvalue<=0.05,]
cis_aTriG05=cis_aTriG[cis_aTriG$pvalue<=0.05,]
trans_eBMI05=trans_eBMI[trans_eBMI$pvalue<=0.05,]
trans_eWHR05=trans_eWHR[trans_eWHR$pvalue<=0.05,]
trans_eT2D05=trans_eT2D[trans_eT2D$pvalue<=0.05,]
trans_eHDL05=trans_eHDL[trans_eHDL$pvalue<=0.05,]
trans_eTriG05=trans_eTriG[trans_eTriG$pvalue<=0.05,]
trans_aBMI05=trans_aBMI[trans_aBMI$pvalue<=0.05,]
trans_aWHR05=trans_aWHR[trans_aWHR$pvalue<=0.05,]
trans_aT2D05=trans_aT2D[trans_aT2D$pvalue<=0.05,]
trans_aHDL05=trans_aHDL[trans_aHDL$pvalue<=0.05,]
trans_aTriG05=trans_aTriG[trans_aTriG$pvalue<=0.05,]

cis_eQTLs05=rbind(cis_eBMI05,cis_eWHR05,cis_eT2D05,cis_eHDL05,cis_eTriG05)
cis_aQTLs05=rbind(cis_aBMI05,cis_aWHR05,cis_aT2D05,cis_aHDL05,cis_aTriG05)
trans_eQTLs05=rbind(trans_eBMI05,trans_eWHR05,trans_eT2D05,trans_eHDL05,trans_eTriG05)
trans_aQTLs05=rbind(trans_aBMI05,trans_aWHR05,trans_aT2D05,trans_aHDL05,trans_aTriG05)

cis_eGenes05=cis_eQTLs05[!duplicated(cis_eQTLs05$gene),]
cis_aGenes05=cis_aQTLs05[!duplicated(cis_aQTLs05$gene),]
trans_eGenes05=trans_eQTLs05[!duplicated(trans_eQTLs05$gene),]
trans_aGenes05=trans_aQTLs05[!duplicated(trans_aQTLs05$gene),]

all_snps_cis_e05=all_snps_cis_eQTLs # These cis-eQTLs were already filtered for P<0.05 by Matrix eQTL
all_snps_cis_a05=all_snps_cis_aQTLs # These cis-aQTLs were already filtered for P<0.05 by Matrix eQTL

all_snps_cis_eGenes05=all_snps_cis_e05[!duplicated(all_snps_cis_e05$gene),]
all_snps_cis_aGenes05=all_snps_cis_a05[!duplicated(all_snps_cis_a05$gene),]

full_cis_e05=full_cis_e # These cis-eQTLs were already filtered for P<0.05 by Matrix eQTL

full_cis_eGenes05=full_cis_e05[!duplicated(full_cis_e05$gene),]

# P<0.0005
cis_eBMI0005=cis_eBMI[cis_eBMI$pvalue<=0.0005,]
cis_eWHR0005=cis_eWHR[cis_eWHR$pvalue<=0.0005,]
cis_eT2D0005=cis_eT2D[cis_eT2D$pvalue<=0.0005,]
cis_eHDL0005=cis_eHDL[cis_eHDL$pvalue<=0.0005,]
cis_eTriG0005=cis_eTriG[cis_eTriG$pvalue<=0.0005,]
cis_aBMI0005=cis_aBMI[cis_aBMI$pvalue<=0.0005,]
cis_aWHR0005=cis_aWHR[cis_aWHR$pvalue<=0.0005,]
cis_aT2D0005=cis_aT2D[cis_aT2D$pvalue<=0.0005,]
cis_aHDL0005=cis_aHDL[cis_aHDL$pvalue<=0.0005,]
cis_aTriG0005=cis_aTriG[cis_aTriG$pvalue<=0.0005,]
trans_eBMI0005=trans_eBMI[trans_eBMI$pvalue<=0.0005,]
trans_eWHR0005=trans_eWHR[trans_eWHR$pvalue<=0.0005,]
trans_eT2D0005=trans_eT2D[trans_eT2D$pvalue<=0.0005,]
trans_eHDL0005=trans_eHDL[trans_eHDL$pvalue<=0.0005,]
trans_eTriG0005=trans_eTriG[trans_eTriG$pvalue<=0.0005,]
trans_aBMI0005=trans_aBMI[trans_aBMI$pvalue<=0.0005,]
trans_aWHR0005=trans_aWHR[trans_aWHR$pvalue<=0.0005,]
trans_aT2D0005=trans_aT2D[trans_aT2D$pvalue<=0.0005,]
trans_aHDL0005=trans_aHDL[trans_aHDL$pvalue<=0.0005,]
trans_aTriG0005=trans_aTriG[trans_aTriG$pvalue<=0.0005,]

cis_eQTLs0005=rbind(cis_eBMI0005,cis_eWHR0005,cis_eT2D0005,cis_eHDL0005,cis_eTriG0005)
cis_aQTLs0005=rbind(cis_aBMI0005,cis_aWHR0005,cis_aT2D0005,cis_aHDL0005,cis_aTriG0005)
trans_eQTLs0005=rbind(trans_eBMI0005,trans_eWHR0005,trans_eT2D0005,trans_eHDL0005,trans_eTriG0005)
trans_aQTLs0005=rbind(trans_aBMI0005,trans_aWHR0005,trans_aT2D0005,trans_aHDL0005,trans_aTriG0005)

cis_eGenes0005=cis_eQTLs0005[!duplicated(cis_eQTLs0005$gene),]
cis_aGenes0005=cis_aQTLs0005[!duplicated(cis_aQTLs0005$gene),]
trans_eGenes0005=trans_eQTLs0005[!duplicated(trans_eQTLs0005$gene),]
trans_aGenes0005=trans_aQTLs0005[!duplicated(trans_aQTLs0005$gene),]

all_snps_cis_e0005=all_snps_cis_eQTLs[all_snps_cis_eQTLs$pvalue<=0.0005,]
all_snps_cis_a0005=all_snps_cis_aQTLs[all_snps_cis_aQTLs$pvalue<=0.0005,]

all_snps_cis_eGenes0005=all_snps_cis_e0005[!duplicated(all_snps_cis_e0005$gene),]
all_snps_cis_aGenes0005=all_snps_cis_a0005[!duplicated(all_snps_cis_a0005$gene),]

full_cis_e0005=full_cis_e[full_cis_e$pvalue<=0.0005,]

full_cis_eGenes0005=full_cis_e0005[!duplicated(full_cis_e0005$gene),]

# P<0.000005
cis_eBMI000005=cis_eBMI[cis_eBMI$pvalue<=0.000005,]
cis_eWHR000005=cis_eWHR[cis_eWHR$pvalue<=0.000005,]
cis_eT2D000005=cis_eT2D[cis_eT2D$pvalue<=0.000005,]
cis_eHDL000005=cis_eHDL[cis_eHDL$pvalue<=0.000005,]
cis_eTriG000005=cis_eTriG[cis_eTriG$pvalue<=0.000005,]
cis_aBMI000005=cis_aBMI[cis_aBMI$pvalue<=0.000005,]
cis_aWHR000005=cis_aWHR[cis_aWHR$pvalue<=0.000005,]
cis_aT2D000005=cis_aT2D[cis_aT2D$pvalue<=0.000005,]
cis_aHDL000005=cis_aHDL[cis_aHDL$pvalue<=0.000005,]
cis_aTriG000005=cis_aTriG[cis_aTriG$pvalue<=0.000005,]
trans_eBMI000005=trans_eBMI[trans_eBMI$pvalue<=0.000005,]
trans_eWHR000005=trans_eWHR[trans_eWHR$pvalue<=0.000005,]
trans_eT2D000005=trans_eT2D[trans_eT2D$pvalue<=0.000005,]
trans_eHDL000005=trans_eHDL[trans_eHDL$pvalue<=0.000005,]
trans_eTriG000005=trans_eTriG[trans_eTriG$pvalue<=0.000005,]
trans_aBMI000005=trans_aBMI[trans_aBMI$pvalue<=0.000005,]
trans_aWHR000005=trans_aWHR[trans_aWHR$pvalue<=0.000005,]
trans_aT2D000005=trans_aT2D[trans_aT2D$pvalue<=0.000005,]
trans_aHDL000005=trans_aHDL[trans_aHDL$pvalue<=0.000005,]
trans_aTriG000005=trans_aTriG[trans_aTriG$pvalue<=0.000005,]

cis_eQTLs000005=rbind(cis_eBMI000005,cis_eWHR000005,cis_eT2D000005,cis_eHDL000005,cis_eTriG000005)
cis_aQTLs000005=rbind(cis_aBMI000005,cis_aWHR000005,cis_aT2D000005,cis_aHDL000005,cis_aTriG000005)
trans_eQTLs000005=rbind(trans_eBMI000005,trans_eWHR000005,trans_eT2D000005,trans_eHDL000005,trans_eTriG000005)
trans_aQTLs000005=rbind(trans_aBMI000005,trans_aWHR000005,trans_aT2D000005,trans_aHDL000005,trans_aTriG000005)

cis_eGenes000005=cis_eQTLs000005[!duplicated(cis_eQTLs000005$gene),]
cis_aGenes000005=cis_aQTLs000005[!duplicated(cis_aQTLs000005$gene),]
trans_eGenes000005=trans_eQTLs000005[!duplicated(trans_eQTLs000005$gene),]
trans_aGenes000005=trans_aQTLs000005[!duplicated(trans_aQTLs000005$gene),]

all_snps_cis_e000005=all_snps_cis_eQTLs[all_snps_cis_eQTLs$pvalue<=0.000005,]
all_snps_cis_a000005=all_snps_cis_aQTLs[all_snps_cis_aQTLs$pvalue<=0.000005,]

all_snps_cis_eGenes000005=all_snps_cis_e000005[!duplicated(all_snps_cis_e000005$gene),]
all_snps_cis_aGenes000005=all_snps_cis_a000005[!duplicated(all_snps_cis_a000005$gene),]

full_cis_e000005=full_cis_e[full_cis_e$pvalue<=0.000005,]

full_cis_eGenes000005=full_cis_e000005[!duplicated(full_cis_e000005$gene),]

# The trans analyses with all regulators and all SNPs, as well as with all genes and all SNPs were already filtered for P<1e-8
all_snps_trans_eGenes=all_snps_trans_eQTLs[!duplicated(all_snps_trans_eQTLs$gene),]
all_snps_trans_aGenes=all_snps_trans_aQTLs[!duplicated(all_snps_trans_aQTLs$gene),]
full_trans_eGenes=full_trans_e[!duplicated(full_trans_e$gene),]

# Now I need to split up genes into EDS bins and filter the regulator list
eds_bin1=eds[eds$EDS<=quantile(eds$EDS,0.2),]
eds_bin2=eds[(eds$EDS>quantile(eds$EDS,0.2) & eds$EDS<=quantile(eds$EDS,0.4)),]
eds_bin3=eds[(eds$EDS>quantile(eds$EDS,0.4) & eds$EDS<=quantile(eds$EDS,0.6)),]
eds_bin4=eds[(eds$EDS>quantile(eds$EDS,0.6) & eds$EDS<=quantile(eds$EDS,0.8)),]
eds_bin5=eds[eds$EDS>quantile(eds$EDS,0.8),]
filt_regs=reg_list[na.omit(match(rownames(eds),reg_list$Symbol)),] # 4218 regulators that are expressed and match the EDS list

# Now let's make the count table for the top and bottom EDS quintiles
count_genes=data.frame("EDS_quintile1"=2649,"EDS_quintile5"=2649,row.names = c("All_expressed_genes"))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% filt_regs),sum(rownames(eds_bin5) %in% filt_regs)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% cis_eGenes05$gene),sum(rownames(eds_bin5) %in% cis_eGenes05$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% cis_eGenes0005$gene),sum(rownames(eds_bin5) %in% cis_eGenes0005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% cis_eGenes000005$gene),sum(rownames(eds_bin5) %in% cis_eGenes000005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% cis_aGenes05$gene),sum(rownames(eds_bin5) %in% cis_aGenes05$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% cis_aGenes0005$gene),sum(rownames(eds_bin5) %in% cis_aGenes0005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% cis_aGenes000005$gene),sum(rownames(eds_bin5) %in% cis_aGenes000005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% all_snps_cis_eGenes05$gene),sum(rownames(eds_bin5) %in% all_snps_cis_eGenes05$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% all_snps_cis_eGenes0005$gene),sum(rownames(eds_bin5) %in% all_snps_cis_eGenes0005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% all_snps_cis_eGenes000005$gene),sum(rownames(eds_bin5) %in% all_snps_cis_eGenes000005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% all_snps_cis_aGenes05$gene),sum(rownames(eds_bin5) %in% all_snps_cis_aGenes05$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% all_snps_cis_aGenes0005$gene),sum(rownames(eds_bin5) %in% all_snps_cis_aGenes0005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% all_snps_cis_aGenes000005$gene),sum(rownames(eds_bin5) %in% all_snps_cis_aGenes000005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% full_cis_eGenes05$gene),sum(rownames(eds_bin5) %in% full_cis_eGenes05$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% full_cis_eGenes0005$gene),sum(rownames(eds_bin5) %in% full_cis_eGenes0005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% full_cis_eGenes000005$gene),sum(rownames(eds_bin5) %in% full_cis_eGenes000005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% trans_eGenes05$gene),sum(rownames(eds_bin5) %in% trans_eGenes05$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% trans_eGenes0005$gene),sum(rownames(eds_bin5) %in% trans_eGenes0005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% trans_eGenes000005$gene),sum(rownames(eds_bin5) %in% trans_eGenes000005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% trans_aGenes05$gene),sum(rownames(eds_bin5) %in% trans_aGenes05$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% trans_aGenes0005$gene),sum(rownames(eds_bin5) %in% trans_aGenes0005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% trans_aGenes000005$gene),sum(rownames(eds_bin5) %in% trans_aGenes000005$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% all_snps_trans_eGenes$gene),sum(rownames(eds_bin5) %in% all_snps_trans_eGenes$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% all_snps_trans_aGenes$gene),sum(rownames(eds_bin5) %in% all_snps_trans_aGenes$gene)))
count_genes=rbind(count_genes,c(sum(rownames(eds_bin1) %in% full_trans_eGenes$gene),sum(rownames(eds_bin5) %in% full_trans_eGenes$gene)))

rownames(count_genes)=c("All_expressed_genes","All_regulators","cis_eGenes05","cis_eGenes0005","cis_eGenes000005",
                        "cis_aGenes05","cis_aGenes0005","cis_aGenes000005","All_SNPs_cis_eGenes05","All_SNPs_cis_eGenes0005",
                        "All_SNPs_cis_eGenes000005","All_SNPs_cis_aGenes05","All_SNPs_cis_aGenes0005","All_SNPs_cis_aGenes000005",
                        "Full_cis_eGenes05","Full_cis_eGenes0005","Full_cis_eGenes000005","trans_eMRs05","trans_eMRs0005","trans_eMRs000005",
                        "trans_aMRs05","trans_aMRs0005","trans_aMRs000005","All_SNPs_trans_eGenes","All_SNPs_trans_aGenes",
                        "Full_trans_eGenes")

# For graphing the changes in proportions in this table, let's calculate the fold difference of EDS_quintile5 percentage relative to All_expressed_genes.
gene_ratio=count_genes$EDS_quintile5/count_genes$EDS_quintile1

par(mar=c(10,4,4,2))
barplot(gene_ratio,names.arg = rownames(count_genes),las=2,cex.names = 0.75)

### Let's compare distributions of EDS scores among QTLs with various P thresholds using a t-test.

# First grab the EDS scores for the genes in the various gene sets
eds_dists=list(NULL)
eds_dists[[1]]=eds$EDS
eds_dists[[2]]=eds[na.omit(match(reg_list$Symbol,rownames(eds))),"EDS"]
eds_dists[[3]]=eds[na.omit(match(cis_eGenes05$gene,rownames(eds))),"EDS"]
eds_dists[[4]]=eds[na.omit(match(cis_eGenes0005$gene,rownames(eds))),"EDS"]
eds_dists[[5]]=eds[na.omit(match(cis_eGenes000005$gene,rownames(eds))),"EDS"]
eds_dists[[6]]=eds[na.omit(match(cis_aGenes05$gene,rownames(eds))),"EDS"]
eds_dists[[7]]=eds[na.omit(match(cis_aGenes0005$gene,rownames(eds))),"EDS"]
eds_dists[[8]]=eds[na.omit(match(cis_aGenes000005$gene,rownames(eds))),"EDS"]
eds_dists[[9]]=eds[na.omit(match(all_snps_cis_eGenes05$gene,rownames(eds))),"EDS"]
eds_dists[[10]]=eds[na.omit(match(all_snps_cis_eGenes0005$gene,rownames(eds))),"EDS"]
eds_dists[[11]]=eds[na.omit(match(all_snps_cis_eGenes000005$gene,rownames(eds))),"EDS"]
eds_dists[[12]]=eds[na.omit(match(all_snps_cis_aGenes05$gene,rownames(eds))),"EDS"]
eds_dists[[13]]=eds[na.omit(match(all_snps_cis_aGenes0005$gene,rownames(eds))),"EDS"]
eds_dists[[14]]=eds[na.omit(match(all_snps_cis_aGenes000005$gene,rownames(eds))),"EDS"]
eds_dists[[15]]=eds[na.omit(match(full_cis_eGenes05$gene,rownames(eds))),"EDS"]
eds_dists[[16]]=eds[na.omit(match(full_cis_eGenes0005$gene,rownames(eds))),"EDS"]
eds_dists[[17]]=eds[na.omit(match(full_cis_eGenes000005$gene,rownames(eds))),"EDS"]
eds_dists[[18]]=eds[na.omit(match(trans_eGenes05$gene,rownames(eds))),"EDS"]
eds_dists[[19]]=eds[na.omit(match(trans_eGenes0005$gene,rownames(eds))),"EDS"]
eds_dists[[20]]=eds[na.omit(match(trans_eGenes000005$gene,rownames(eds))),"EDS"]
eds_dists[[21]]=eds[na.omit(match(trans_aGenes05$gene,rownames(eds))),"EDS"]
eds_dists[[22]]=eds[na.omit(match(trans_aGenes0005$gene,rownames(eds))),"EDS"]
eds_dists[[23]]=eds[na.omit(match(trans_aGenes000005$gene,rownames(eds))),"EDS"]
eds_dists[[24]]=eds[na.omit(match(all_snps_trans_eGenes$gene,rownames(eds))),"EDS"]
eds_dists[[25]]=eds[na.omit(match(all_snps_trans_aGenes$gene,rownames(eds))),"EDS"]
eds_dists[[26]]=eds[na.omit(match(full_trans_eGenes$gene,rownames(eds))),"EDS"]
names(eds_dists)=c("All_expressed_genes","All_regulators","cis_eGenes05","cis_eGenes0005","cis_eGenes000005",
                   "cis_aGenes05","cis_aGenes0005","cis_aGenes000005","All_SNPs_cis_eGenes05","All_SNPs_cis_eGenes0005",
                   "All_SNPs_cis_eGenes000005","All_SNPs_cis_aGenes05","All_SNPs_cis_aGenes0005","All_SNPs_cis_aGenes000005",
                   "Full_cis_eGenes05","Full_cis_eGenes0005","Full_cis_eGenes000005","trans_eMRs05","trans_eMRs0005","trans_eMRs000005",
                   "trans_aMRs05","trans_aMRs0005","trans_aMRs000005","All_SNPs_trans_eGenes","All_SNPs_trans_aGenes",
                   "Full_trans_eGenes")

# Let's plot these with a boxplots
par(mar=c(13,4,2,2))
boxplot(eds_dists,las=2,cex.names = 0.75)

# t-tests
t.test(eds_dists$All_expressed_genes,eds_dists$All_regulators) # P=1.20e-32

# cis-QTLs vs all expressed genes
t.test(eds_dists$All_expressed_genes,eds_dists$cis_eGenes05) # P=8.23e-19
t.test(eds_dists$All_expressed_genes,eds_dists$cis_eGenes0005) # P=0.0001661
t.test(eds_dists$All_expressed_genes,eds_dists$cis_eGenes000005) # P=0.09787
t.test(eds_dists$All_expressed_genes,eds_dists$cis_aGenes05) # P=1.08e-14
t.test(eds_dists$All_expressed_genes,eds_dists$cis_aGenes0005) # P=0.02258
t.test(eds_dists$All_expressed_genes,eds_dists$cis_aGenes000005) # P=0.4001 very poorly powered given only 6 aQTLs at this threshold
t.test(eds_dists$All_expressed_genes,eds_dists$All_SNPs_cis_eGenes05) # P=4.81e-31
t.test(eds_dists$All_expressed_genes,eds_dists$All_SNPs_cis_eGenes0005) # P=1.29e-25
t.test(eds_dists$All_expressed_genes,eds_dists$All_SNPs_cis_eGenes000005) # P=0.0002032
t.test(eds_dists$All_expressed_genes,eds_dists$All_SNPs_cis_aGenes05) # P=4.12e-31
t.test(eds_dists$All_expressed_genes,eds_dists$All_SNPs_cis_aGenes0005) # P=2.03e-18
t.test(eds_dists$All_expressed_genes,eds_dists$All_SNPs_cis_aGenes000005) # P=0.05396
t.test(eds_dists$All_expressed_genes,eds_dists$Full_cis_eGenes05) # P=0.7931
t.test(eds_dists$All_expressed_genes,eds_dists$Full_cis_eGenes0005) # P=0.9869
t.test(eds_dists$All_expressed_genes,eds_dists$Full_cis_eGenes000005) # P=9.38e-5

# cis-QTLs vs all regulators
t.test(eds_dists$All_regulators,eds_dists$cis_eGenes05) # P=0.08145
t.test(eds_dists$All_regulators,eds_dists$cis_eGenes0005) # P=0.7253
t.test(eds_dists$All_regulators,eds_dists$cis_eGenes000005) # P=0.2969
t.test(eds_dists$All_regulators,eds_dists$cis_aGenes05) # P=0.1597
t.test(eds_dists$All_regulators,eds_dists$cis_aGenes0005) # P=0.2858
t.test(eds_dists$All_regulators,eds_dists$cis_aGenes000005) # P=0.7648 very poorly powered given only 6 aQTLs at this threshold
t.test(eds_dists$All_regulators,eds_dists$All_SNPs_cis_eGenes05) # P=0.9432
t.test(eds_dists$All_regulators,eds_dists$All_SNPs_cis_eGenes0005) # P=0.4648
t.test(eds_dists$All_regulators,eds_dists$All_SNPs_cis_eGenes000005) # P=0.002172
t.test(eds_dists$All_regulators,eds_dists$All_SNPs_cis_aGenes05) # P=0.9583
t.test(eds_dists$All_regulators,eds_dists$All_SNPs_cis_aGenes0005) # P=0.1066
t.test(eds_dists$All_regulators,eds_dists$All_SNPs_cis_aGenes000005) # P=0.6494
t.test(eds_dists$All_regulators,eds_dists$Full_cis_eGenes05) # P=2.66e-33
t.test(eds_dists$All_regulators,eds_dists$Full_cis_eGenes0005) # P=6.83e-29
t.test(eds_dists$All_regulators,eds_dists$Full_cis_eGenes000005) # P=2.16e-38

# trans-QTLs vs all expressed genes
t.test(eds_dists$All_expressed_genes,eds_dists$trans_eMRs05) # P=5.13e-22 This is essentially all MRs vs all expressed genes
t.test(eds_dists$All_expressed_genes,eds_dists$trans_eMRs0005) # P=8.47e-19
t.test(eds_dists$All_expressed_genes,eds_dists$trans_eMRs000005) # P=0.0003069
t.test(eds_dists$All_expressed_genes,eds_dists$trans_aMRs05) # P=5.13e-22 This is essentially all MRs vs all expressed genes
t.test(eds_dists$All_expressed_genes,eds_dists$trans_aMRs0005) # P=2.09e-18
t.test(eds_dists$All_expressed_genes,eds_dists$trans_aMRs000005) # P=1.57e-8
t.test(eds_dists$All_expressed_genes,eds_dists$All_SNPs_trans_eGenes) # P=3.14e-9
t.test(eds_dists$All_expressed_genes,eds_dists$All_SNPs_trans_aGenes) # P=5.38e-12
t.test(eds_dists$All_expressed_genes,eds_dists$Full_trans_eGenes) # P=0.0007574

# trans-QTLs vs all regulators
t.test(eds_dists$All_regulators,eds_dists$trans_eMRs05) # P=4.25e-11 This is essentially all MRs vs all regulators
t.test(eds_dists$All_regulators,eds_dists$trans_eMRs0005) # P=2.64e-10
t.test(eds_dists$All_regulators,eds_dists$trans_eMRs000005) # P=0.01171
t.test(eds_dists$All_regulators,eds_dists$trans_aMRs05) # P=4.25e-11 This is essentially all MRs vs all regulators
t.test(eds_dists$All_regulators,eds_dists$trans_aMRs0005) # P=1.93e-9
t.test(eds_dists$All_regulators,eds_dists$trans_aMRs000005) # P=1.71e-5
t.test(eds_dists$All_regulators,eds_dists$All_SNPs_trans_eGenes) # P=0.009747
t.test(eds_dists$All_regulators,eds_dists$All_SNPs_trans_aGenes) # P=0.0001215
t.test(eds_dists$All_regulators,eds_dists$Full_trans_eGenes) # P=0.00598

# Targeted comparisons
t.test(eds_dists$cis_eGenes05,eds_dists$cis_eGenes000005) # P=0.09222
t.test(eds_dists$cis_aGenes05,eds_dists$cis_aGenes000005) # P=0.8706 very poorly powered given only 6 aQTLs at 0.000005 P threshold
t.test(eds_dists$cis_eGenes05,eds_dists$cis_aGenes05) # P=0.8825
t.test(eds_dists$cis_eGenes0005,eds_dists$cis_aGenes0005) # P=0.3756
t.test(eds_dists$cis_eGenes000005,eds_dists$cis_aGenes000005) # P=0.6116 very poorly powered given only 6 aQTLs at 0.000005 P threshold
t.test(eds_dists$All_SNPs_cis_eGenes05,eds_dists$All_SNPs_cis_eGenes000005) # P=0.002698
t.test(eds_dists$All_SNPs_cis_aGenes05,eds_dists$All_SNPs_cis_aGenes000005) # P=0.6436
t.test(eds_dists$All_SNPs_cis_eGenes05,eds_dists$All_SNPs_cis_aGenes05) # P=985
t.test(eds_dists$All_SNPs_cis_eGenes0005,eds_dists$All_SNPs_cis_aGenes0005) # P=0.332
t.test(eds_dists$All_SNPs_cis_eGenes000005,eds_dists$All_SNPs_cis_aGenes000005) # P=0.2537
t.test(eds_dists$Full_cis_eGenes05,eds_dists$Full_cis_eGenes000005) # P=0.0002145
t.test(eds_dists$trans_eMRs05,eds_dists$trans_eMRs000005) # P=0.9858
t.test(eds_dists$trans_aMRs05,eds_dists$trans_aMRs000005) # P=0.3333
t.test(eds_dists$trans_eMRs05,eds_dists$trans_aMRs05) # P=1 They are the same full set of MRs
t.test(eds_dists$trans_eMRs0005,eds_dists$trans_aMRs0005) # P=0.6166
t.test(eds_dists$trans_eMRs000005,eds_dists$trans_aMRs000005) # P=0.5294
t.test(eds_dists$trans_eMRs000005,eds_dists$All_SNPs_trans_eGenes) # P=0.1315
t.test(eds_dists$trans_aMRs000005,eds_dists$All_SNPs_trans_aGenes) # P=0.0313
t.test(eds_dists$Full_trans_eGenes,eds_dists$trans_eMRs000005) # P=0.002587
t.test(eds_dists$Full_trans_eGenes,eds_dists$trans_aMRs000005) # P=1.07e-6
t.test(eds_dists$Full_trans_eGenes,eds_dists$All_SNPs_trans_eGenes) # P=0.0001629
t.test(eds_dists$Full_trans_eGenes,eds_dists$All_SNPs_trans_aGenes) # P=1.10e-6
t.test(eds_dists$All_SNPs_trans_eGenes,eds_dists$All_SNPs_trans_aGenes) # P=0.2839

# Here I try splitting the cis-e/aGenes into P<5e-6, 5e-6<P<5e-4 and 5e-4<P<5e-2 bins. I will keep the other trans analyses as well even though 
# they are at a different P threshold.

cis_eGenes_bin1=cis_eGenes05[cis_eGenes05$pvalue>0.0005,]
cis_eGenes_bin2=cis_eGenes05[(cis_eGenes05$pvalue<0.0005) & (cis_eGenes05$pvalue>0.000005),]
cis_eGenes_bin3=cis_eGenes05[cis_eGenes05$pvalue<0.000005,]
cis_aGenes_bin1=cis_aGenes05[cis_aGenes05$pvalue>0.0005,]
cis_aGenes_bin2=cis_aGenes05[(cis_aGenes05$pvalue<0.0005) & (cis_aGenes05$pvalue>0.000005),]
cis_aGenes_bin3=cis_aGenes05[cis_aGenes05$pvalue<0.000005,]
all_snps_cis_eGenes_bin1=all_snps_cis_eGenes05[all_snps_cis_eGenes05$pvalue>0.0005,]
all_snps_cis_eGenes_bin2=all_snps_cis_eGenes05[(all_snps_cis_eGenes05$pvalue<0.0005) & (all_snps_cis_eGenes05$pvalue>0.000005),]
all_snps_cis_eGenes_bin3=all_snps_cis_eGenes05[all_snps_cis_eGenes05$pvalue<0.000005,]
all_snps_cis_aGenes_bin1=all_snps_cis_aGenes05[all_snps_cis_aGenes05$pvalue>0.0005,]
all_snps_cis_aGenes_bin2=all_snps_cis_aGenes05[(all_snps_cis_aGenes05$pvalue<0.0005) & (all_snps_cis_aGenes05$pvalue>0.000005),]
all_snps_cis_aGenes_bin3=all_snps_cis_aGenes05[all_snps_cis_aGenes05$pvalue<0.000005,]
full_cis_eGenes_bin1=full_cis_eGenes05[full_cis_eGenes05$pvalue>0.0005,]
full_cis_eGenes_bin2=full_cis_eGenes05[(full_cis_eGenes05$pvalue<0.0005) & (full_cis_eGenes05$pvalue>0.000005),]
full_cis_eGenes_bin3=full_cis_eGenes05[full_cis_eGenes05$pvalue<0.000005,]
trans_eGenes_bin1=trans_eGenes05[trans_eGenes05$pvalue>0.0005,]
trans_eGenes_bin2=trans_eGenes05[(trans_eGenes05$pvalue<0.0005) & (trans_eGenes05$pvalue>0.000005),]
trans_eGenes_bin3=trans_eGenes05[trans_eGenes05$pvalue<0.000005,]
trans_aGenes_bin1=trans_aGenes05[trans_aGenes05$pvalue>0.0005,]
trans_aGenes_bin2=trans_aGenes05[(trans_aGenes05$pvalue<0.0005) & (trans_aGenes05$pvalue>0.000005),]
trans_aGenes_bin3=trans_aGenes05[trans_aGenes05$pvalue<0.000005,]

strat_eds_dists=list("All expressed genes"=eds$EDS,"All regulators"=eds[na.omit(match(reg_list$Symbol,rownames(eds))),"EDS"],
                    "cis-eGenes bin1"=cis_eGenes_bin1$EDS,"cis-eGenes bin2"=cis_eGenes_bin2$EDS,"cis-eGenes bin3"=cis_eGenes_bin3$EDS,
                    "cis-aGenes bin1"=cis_aGenes_bin1$EDS,"cis-aGenes bin2"=cis_aGenes_bin2$EDS,"cis-aGenes bin3"=cis_aGenes_bin3$EDS,
                    "All SNPs cis-eGenes bin1"=all_snps_cis_eGenes_bin1$EDS,"All SNPs cis-eGenes bin2"=all_snps_cis_eGenes_bin2$EDS,"All SNPs cis-eGenes bin3"=all_snps_cis_eGenes_bin3$EDS,
                    "All SNPs cis-aGenes bin1"=all_snps_cis_aGenes_bin1$EDS,"All SNPs cis-aGenes bin2"=all_snps_cis_aGenes_bin2$EDS,"All SNPs cis-aGenes bin3"=all_snps_cis_aGenes_bin3$EDS,
                    "Full cis-eGenes bin1"=full_cis_eGenes_bin1$EDS,"Full cis-eGenes bin2"=full_cis_eGenes_bin2$EDS,"Full cis-eGenes bin3"=full_cis_eGenes_bin3$EDS,
                    "trans-eMRs bin1"=trans_eGenes_bin1$EDS,"trans-eMRs bin2"=trans_eGenes_bin2$EDS,"trans-eMRs bin3"=trans_eGenes_bin3$EDS,
                    "trans-aMRs bin1"=trans_aGenes_bin1$EDS,"trans-aMRs bin2"=trans_aGenes_bin2$EDS,"trans-aMRs bin3"=trans_aGenes_bin3$EDS,
                    "All SNPs trans-eGenes"=all_snps_trans_eGenes$EDS,"All SNPs trans-aGenes"=all_snps_trans_aGenes$EDS,"Full trans-eGenes"=full_trans_eGenes$EDS)

par(mar=c(13,4,2,2))
boxplot(strat_eds_dists,las=2,cex.names = 0.75)

# I need the data in strat_eds_dists in a data.frame format for ggplot
strat_eds_df=data.frame("Gene_set"=NULL, "EDS_scores"=NULL, "Analysis_set"=NULL,"Set_type"=NULL)
for(i in names(strat_eds_dists)){
  if(grepl('-e',i)){
    set_type="eGenes"
  } else{
      if(grepl('-a',i)){
        set_type="aGenes"
      } else{
          set_type="Reference genes"
      }
    }
  analysis_set=gsub(' bin\\d','',i)
  temp=data.frame("Gene_set"=rep(i,length(strat_eds_dists[[i]])),
                  "EDS_scores"=strat_eds_dists[[i]],
                  "Analysis_set"=rep(analysis_set,length(strat_eds_dists[[i]])),
                  "Set_type"=rep(set_type,length(strat_eds_dists[[i]])))
  strat_eds_df=rbind(strat_eds_df,temp)
}
strat_eds_df$Gene_set=factor(strat_eds_df$Gene_set,levels = unique(strat_eds_df$Gene_set)) # This prevents ggplot from reordering my categorical variables alphabetically
strat_eds_df$Analysis_set=factor(strat_eds_df$Analysis_set,levels = unique(strat_eds_df$Analysis_set)) # This prevents ggplot from reordering my categorical variables alphabetically
strat_eds_df$Set_type=factor(strat_eds_df$Set_type,levels = unique(strat_eds_df$Set_type)) # This prevents ggplot from reordering my categorical variables alphabetically

test=ggplot(strat_eds_df,aes(x=Gene_set,y=EDS_scores,fill=Set_type)) +
  geom_violin(trim=F,position = position_dodge(0.8),width = 0.7) +
  scale_fill_manual(breaks=unique(strat_eds_df$Set_type),values = c("#B2BEB5","#0099FF","#FF6666")) +
  labs(x="Gene set",y="EDS score") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90,vjust = 0.5,hjust = 1))
test

test2=ggplot(strat_eds_df[(strat_eds_df$Analysis_set!="cis-eGenes" & strat_eds_df$Analysis_set!="cis-aGenes"),],aes(x=Gene_set,y=EDS_scores,fill=Set_type)) +
  geom_boxplot(notch = T,outlier.size = 1) +
  scale_fill_manual(values = c("#B2BEB5","#0099FF","#FF6666")) +
  labs(x="Gene set",y="EDS score") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90,vjust = 0.5,hjust = 1))
test2

cis_strat_eds_df=strat_eds_df[strat_eds_df$Gene_set %in% 
                                 c("All expressed genes","All regulators","cis-eGenes bin3","cis-aGenes bin3",
                                   "All SNPs cis-eGenes bin3","All SNPs cis-aGenes bin3","Full cis-eGenes bin3"),]
cis_strat_eds_df$Gene_set=factor(cis_strat_eds_df$Gene_set,levels = 
                                    c("All expressed genes","All regulators","Full cis-eGenes bin3","All SNPs cis-eGenes bin3",
                                      "All SNPs cis-aGenes bin3","cis-eGenes bin3","cis-aGenes bin3")) # This will set the ordering of my categorical variables as I want them in the ggplot
trans_strat_eds_df=strat_eds_df[strat_eds_df$Gene_set %in% 
                                c("All expressed genes","All regulators","trans-eMRs bin3","trans-aMRs bin3",
                                  "All SNPs trans-eGenes","All SNPs trans-aGenes","Full trans-eGenes"),]
trans_strat_eds_df$Gene_set=factor(trans_strat_eds_df$Gene_set,levels = 
                                   c("All expressed genes","All regulators","Full trans-eGenes","All SNPs trans-eGenes",
                                     "All SNPs trans-aGenes","trans-eMRs bin3","trans-aMRs bin3")) # This will set the ordering of my categorical variables as I want them in the ggplot

# Build some final boxplots with ggplot2
cis_box=ggplot(cis_strat_eds_df,aes(x=Gene_set,y=EDS_scores,fill=Set_type)) +
  geom_boxplot(outlier.size = 1) +
  scale_fill_manual(values = c("#B2BEB5","#0099FF","#FF6666")) +
  labs(x="Gene set",y="EDS score") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90,vjust = 0.5,hjust = 1))


trans_box=ggplot(trans_strat_eds_df,aes(x=Gene_set,y=EDS_scores,fill=Set_type)) +
  geom_boxplot(outlier.size = 1) +
  scale_fill_manual(values = c("#B2BEB5","#0099FF","#FF6666")) +
  labs(x="Gene set",y="EDS score") + 
  stat_summary(fun.data = "median_hilow",fun.args = list(conf.int=0.5),position = position_dodge(0.8)) +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 90,vjust = 0.5,hjust = 1))

# Write the plots to files
pdf("EDS_scores_for_cis-genes_boxplot.pdf",width = 8)
cis_box
dev.off()

pdf("EDS_scores_for_trans-genes_boxplot.pdf",width = 8)
trans_box
dev.off()

# t-tests
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`All regulators`) # P=1.20e-32

# cis-QTLs vs all expressed genes
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`cis-eGenes bin1`) # P=2.88e-16, mean of y=0.5356
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`cis-eGenes bin2`) # P=0.000671, mean of y=0.5464
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`cis-eGenes bin3`) # P=0.09718, mean of y=0.5160
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`cis-aGenes bin1`) # P=2.25e-13, mean of y=0.5324
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`cis-aGenes bin2`) # P=0.01179, mean of y=0.5706
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`cis-aGenes bin3`) # P=0.2741, mean of y=0.5657, very poorly powered given only 5 aGenes in this bin
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`All SNPs cis-eGenes bin1`) # P=2.24e-10, mean of y=0.5211
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`All SNPs cis-eGenes bin2`) # P=5.83e-29, mean of y=0.5439
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`All SNPs cis-eGenes bin3`) # P=0.0002032, mean of y=0.5119
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`All SNPs cis-aGenes bin1`) # P=3.63e-18, mean of y=0.5223
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`All SNPs cis-aGenes bin2`) # P=2.03e-18, mean of y=0.5400
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`All SNPs cis-aGenes bin3`) # P=0.05396, mean of y=0.5357
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`Full cis-eGenes bin1`) # P=0.5682, mean of y=0.4945
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`Full cis-eGenes bin2`) # P=0.8.81e-5, mean of y=0.5059
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`Full cis-eGenes bin3`) # P=9.38e-5, mean of y=0.4864

# cis-QTLs vs all regulators
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`cis-eGenes bin1`) # P=0.07609, mean of y=0.5355
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`cis-eGenes bin2`) # P=0.1738, mean of y=0.5464
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`cis-eGenes bin3`) # P=0.395, mean of y=0.5160
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`cis-aGenes bin1`) # P=0.2613, mean of y=0.5324
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`cis-aGenes bin2`) # P=0.1244, mean of y=0.5706
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`cis-aGenes bin3`) # P=0.5159, mean of y=0.5657, very poorly powered given only 5 aGenes in this bin
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`All SNPs cis-eGenes bin1`) # P=0.2291, mean of y=0.5211
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`All SNPs cis-eGenes bin2`) # P=0.0001678, mean of y=0.5439
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`All SNPs cis-eGenes bin3`) # P=0.002172, mean of y=0.5119
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`All SNPs cis-aGenes bin1`) # P=0.2475, mean of y=0.5223
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`All SNPs cis-aGenes bin2`) # P=0.1161, mean of y=0.5400
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`All SNPs cis-aGenes bin3`) # P=0.6494, mean of y=0.5357
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`Full cis-eGenes bin1`) # P=1.69e-25, mean of y=0.4945
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`Full cis-eGenes bin2`) # P=1.30e-10, mean of y=0.5059
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`Full cis-eGenes bin3`) # P=2.16e-38, mean of y=0.4864

# trans-QTLs vs all expressed genes
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`trans-eMRs bin1`) # P=2.12e-5, mean of y=0.5648, This is no longer all the MRs due to the stratification
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`trans-eMRs bin2`) # P=2.44e-15, mean of y=0.5985
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`trans-eMRs bin3`) # P=0.0003081, mean of y=0.5950
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`trans-aMRs bin1`) # P=0.0001418, mean of y=0.5758, This is no longer all the MRs due to the stratification
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`trans-aMRs bin2`) # P=2.29e-12, mean of y=0.5861
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`trans-aMRs bin3`) # P=3.21e-8, mean of y=0.6032
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`All SNPs trans-eGenes`) # P=3.14e-9, mean of y=0.5502
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`All SNPs trans-aGenes`) # P=5.38e-12, mean of y=0.5640
t.test(strat_eds_dists$`All expressed genes`,strat_eds_dists$`Full trans-eGenes`) # P=0.0007574, mean of y=0.5122

# trans-QTLs vs all regulators
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`trans-eMRs bin1`) # P=0.01488, mean of y=0.5648, This is no longer all the MRs due to the stratification
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`trans-eMRs bin2`) # P=8.51e-9, mean of y=0.5985
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`trans-eMRs bin3`) # P=0.009181, mean of y=0.5950
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`trans-aMRs bin1`) # P=0.01501, mean of y=0.5758, This is no longer all the MRs due to the stratification
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`trans-aMRs bin2`) # P=1.77e-6, mean of y=0.5861
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`trans-aMRs bin3`) # P=3.29e-5, mean of y=0.6032
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`All SNPs trans-eGenes`) # P=0.009747, mean of y=0.5502
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`All SNPs trans-aGenes`) # P=0.0001215, mean of y=0.5640
t.test(strat_eds_dists$`All regulators`,strat_eds_dists$`Full trans-eGenes`) # P=0.00598, mean of y=0.5122

# Targeted comparisons
t.test(strat_eds_dists$`cis-eGenes bin1`,strat_eds_dists$`cis-eGenes bin3`) # P=0.1303, mean of x=0.5356 and y=0.5160
t.test(strat_eds_dists$`cis-aGenes bin1`,strat_eds_dists$`cis-aGenes bin3`) # P=0.5788, mean of x=0.5324 and y=0.5657, very poorly powered given only 5 aGenes in this bin
t.test(strat_eds_dists$`cis-eGenes bin1`,strat_eds_dists$`cis-aGenes bin1`) # P=0.6278, mean of x=0.5356 and y=0.5324
t.test(strat_eds_dists$`cis-eGenes bin2`,strat_eds_dists$`cis-aGenes bin2`) # P=0.4445, mean of x=0.5464 and y=0.5706
t.test(strat_eds_dists$`cis-eGenes bin3`,strat_eds_dists$`cis-aGenes bin3`) # P=0.4242, mean of x=0.5160 and y=0.5657 very poorly powered given only 5 aGenes in this bin
t.test(strat_eds_dists$`All SNPs cis-eGenes bin1`,strat_eds_dists$`All SNPs cis-eGenes bin3`) # P=0.101, mean of x=0.5648 and y=0.5119
t.test(strat_eds_dists$`All SNPs cis-aGenes bin1`,strat_eds_dists$`All SNPs cis-aGenes bin3`) # P=0.5131, mean of x=0.5223 and y=0.5357
t.test(strat_eds_dists$`All SNPs cis-eGenes bin1`,strat_eds_dists$`All SNPs cis-aGenes bin1`) # P=0.8056, mean of x=0.5211 and y=0.5223
t.test(strat_eds_dists$`All SNPs cis-eGenes bin2`,strat_eds_dists$`All SNPs cis-aGenes bin2`) # P=0.09032, mean of x=0.5439 and y=0.5400
t.test(strat_eds_dists$`All SNPs cis-eGenes bin3`,strat_eds_dists$`All SNPs cis-aGenes bin3`) # P=0.2537, mean of x=0.5119 and y=0.5357
t.test(strat_eds_dists$`Full cis-eGenes bin1`,strat_eds_dists$`Full cis-eGenes bin3`) # P=0.005884, mean of x=0.4945 and y=0.4864
t.test(strat_eds_dists$`trans-eMRs bin1`,strat_eds_dists$`trans-eMRs bin3`) # P=0.3049, mean of x=0.5648 and y=0.5950
t.test(strat_eds_dists$`trans-aMRs bin1`,strat_eds_dists$`trans-aMRs bin3`) # P=0.2943, mean of x=0.5758 and y=0.6032
t.test(strat_eds_dists$`trans-eMRs bin1`,strat_eds_dists$`trans-aMRs bin1`) # P=0.6589, mean of x=0.5648 and y=0.5758
t.test(strat_eds_dists$`trans-eMRs bin2`,strat_eds_dists$`trans-aMRs bin2`) # P=0.4567, mean of x=0.5985 and y=0.5861
t.test(strat_eds_dists$`trans-eMRs bin3`,strat_eds_dists$`trans-aMRs bin3`) # P=0.7852, mean of x=0.5950 and y=0.6032
t.test(strat_eds_dists$`All SNPs trans-eGenes`,strat_eds_dists$`Full trans-eGenes`) # P=0.0001629, mean of x=0.5502 and y=0.5122
t.test(strat_eds_dists$`All SNPs trans-aGenes`,strat_eds_dists$`Full trans-eGenes`) # P=1.10e-6, mean of x=0.5640 and y=0.5122
t.test(strat_eds_dists$`All SNPs trans-eGenes`,strat_eds_dists$`All SNPs trans-aGenes`) # P=0.2839, mean of x=0.5502 and y=0.5640
t.test(strat_eds_dists$`All SNPs trans-eGenes`,strat_eds_dists$`trans-eMRs bin3`) # P=0.09564, mean of x=0.5502 and y=0.5950
t.test(strat_eds_dists$`All SNPs trans-aGenes`,strat_eds_dists$`trans-aMRs bin3`) # P=0.04682, mean of x=0.5640 and y=0.6032
t.test(strat_eds_dists$`Full trans-eGenes`,strat_eds_dists$`trans-eMRs bin3`) # P=0.002203, mean of x=0.5122 and y=0.5950
t.test(strat_eds_dists$`Full trans-eGenes`,strat_eds_dists$`trans-aMRs bin3`) # P=2.13e-6, mean of x=0.5122 and y=0.6032
t.test(strat_eds_dists$`trans-eMRs bin3`,strat_eds_dists$`trans-aMRs bin3`) # P=0.7852, mean of x=0.5950 and y=0.6032

# Ultimately, for stringency and simplicity of reporting, I have decided to stick to the bin3 analyses (i.e., p<0.000005).

# As a supplemental table I want the summaries on some of these distributions
strat_eds_sum=rbind("All expressed genes"=summary(eds$EDS),"All regulators"=summary(eds[na.omit(match(reg_list$Symbol,rownames(eds))),"EDS"]),
                    "cis-eGenes"=summary(full_cis_eGenes_bin3$EDS),
                    "cis-eRegulators"=summary(all_snps_cis_eGenes_bin3$EDS),"cis-aRegulators"=summary(all_snps_cis_aGenes_bin3$EDS),
                    "cis-eRegulators (GWAS SNPs)"=summary(cis_eGenes_bin3$EDS),"cis-aRegulators (GWAS SNPs)"=summary(cis_aGenes_bin3$EDS),
                    "trans-eGenes"=summary(full_trans_eGenes$EDS),
                    "trans-eRegulators"=summary(all_snps_trans_eGenes$EDS),"trans-aRegulators"=summary(all_snps_trans_aGenes$EDS),
                    "trans-eMRs (GWAS SNPs)"=summary(trans_eGenes_bin3$EDS),"trans-aMRs (GWAS SNPs)"=summary(trans_aGenes_bin3$EDS)
)
strat_eds_sum=as.data.frame(strat_eds_sum[,-7]) # There were apparently NAs in 2 of these sets, which made for an extra column.

write.table(cbind(rownames(strat_eds_sum),strat_eds_sum),"Summary_stats_on_gene_set_EDS_distributions.txt",
            sep = "\t",row.names = F,col.names = T,quote = F)



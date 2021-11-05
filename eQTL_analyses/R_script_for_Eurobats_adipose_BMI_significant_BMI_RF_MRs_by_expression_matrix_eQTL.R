#install.packages("MatrixEQTL")

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = "../BMI_significant_SNPs_QCd_filtered.dosage";
snps_location_file_name = "../BMI_significant_SNPs_locations_QCd.txt";

# Gene expression file name
expression_file_name = "../../Adipose\ expression\ data/FINAL_logTPMs_and_activities/Filtered_Eurobats_adipose_qnorm_INT_logTPMs_for_BMI_RF_MRs_by_expression.txt";
gene_location_file_name = "../../Adipose\ expression\ data/FINAL_logTPMs_and_activities/Hg19_gene_map_for_13776_expressed_genes_in_Eurobats_adipose.map";

# Covariates file name
# Set to character() for no covariates
covariates_file_name = "../Filtered_Eurobats_adipose_covars_no_PEER.txt";

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 1;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);


## Normal quantile transformation of gene expression data

for( sl in 1:length(gene) ) {
  mat = gene[[sl]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(gene)+1));
  gene[[sl]] = mat;
}
rm(sl, mat);


## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
cis_eqtls<-me$cis$eqtls
cat('Detected distant eQTLs:', '\n');
trans_eqtls<-me$trans$eqtls

## Plot the Q-Q plot of local and distant p-values
jpeg("Eurobats_adipose_BMI_sig_eQTLs_from_INT_logTPM_of_BMI_RF_MRs_by_expression_plot.jpg")
plot(me)
dev.off()

write.table(cis_eqtls,"Eurobats_adipose_BMI_sig_cis-eQTLs_from_INT_logTPM_of_BMI_RF_MRs_by_expression.txt",sep="\t",quote = FALSE,row.names=FALSE)
write.table(trans_eqtls,"Eurobats_adipose_BMI_sig_trans-eQTLs_from_INT_logTPM_of_BMI_RF_MRs_by_expression.txt",sep="\t",quote = FALSE,row.names=FALSE)

q(save="no")
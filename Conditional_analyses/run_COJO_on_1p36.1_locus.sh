#!/bin/bash

# Analysis requires gcta to be installed and available.
module load gcta

# A few quick notes on running COJO. First, make sure to include the full summary statistics as the COJO file rather than a subset for a locus or loci of interest. 
# The analysis can then be focused to specific SNPs with --extract. The reference set data can be restricted to just the SNPs of interest. 
# Also, the alleles need to match between the reference data and the GWAS summary statistics according to Jian Yang, one of COJO's creators. Consequently, I had to edit the summary statistics for my SNPs of interest.

awk '{print $3}' Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED_locus_1p36.1_with_locations.txt > chr1p36.1_snplist.txt
gcta64 --bfile /DCEG/Branches/LTG/Amundadottir/PanScan2019/delivery.10.24.2019/imputation/HRC/P12/chr1p36.1 --cojo-file Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED_flipped_1p36.1_alleles_no_positions.txt --extract chr1p36.1_snplist.txt --cojo-cond rs6692586.txt --cojo-p 5e-2 --out Meta-analysis_Locke_et_al_1p36.1_COJO_conditioned_on_rs6692586
gcta64 --bfile /DCEG/Branches/LTG/Amundadottir/PanScan2019/delivery.10.24.2019/imputation/HRC/P12/chr1p36.1 --cojo-file Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED_flipped_1p36.1_alleles_no_positions.txt --extract chr1p36.1_snplist.txt --cojo-cond rs4654828.txt --cojo-p 5e-2 --out Meta-analysis_Locke_et_al_1p36.1_COJO_conditioned_on_rs4654828
gcta64 --bfile /DCEG/Branches/LTG/Amundadottir/PanScan2019/delivery.10.24.2019/imputation/HRC/P12/chr1p36.1 --cojo-file Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED_flipped_1p36.1_alleles_no_positions.txt --extract chr1p36.1_snplist.txt --cojo-cond rs12408468.txt --cojo-p 5e-2 --out Meta-analysis_Locke_et_al_1p36.1_COJO_conditioned_on_rs12408468
gcta64 --bfile /DCEG/Branches/LTG/Amundadottir/PanScan2019/delivery.10.24.2019/imputation/HRC/P12/chr1p36.1 --cojo-file Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED_flipped_1p36.1_alleles_no_positions.txt --extract chr1p36.1_snplist.txt --cojo-cond rs6692586+rs4654828.txt --cojo-p 5e-2 --out Meta-analysis_Locke_et_al_1p36.1_COJO_conditioned_on_rs6692586+rs4654828
gcta64 --bfile /DCEG/Branches/LTG/Amundadottir/PanScan2019/delivery.10.24.2019/imputation/HRC/P12/chr1p36.1 --cojo-file Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED_flipped_1p36.1_alleles_no_positions.txt --extract chr1p36.1_snplist.txt --cojo-cond rs6692586+rs12408468.txt --cojo-p 5e-2 --out Meta-analysis_Locke_et_al_1p36.1_COJO_conditioned_on_rs6692586+rs12408468
gcta64 --bfile /DCEG/Branches/LTG/Amundadottir/PanScan2019/delivery.10.24.2019/imputation/HRC/P12/chr1p36.1 --cojo-file Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED_flipped_1p36.1_alleles_no_positions.txt --extract chr1p36.1_snplist.txt --cojo-cond rs4654828+rs12408468.txt --cojo-p 5e-2 --out Meta-analysis_Locke_et_al_1p36.1_COJO_conditioned_on_rs4654828+rs12408468
gcta64 --bfile /DCEG/Branches/LTG/Amundadottir/PanScan2019/delivery.10.24.2019/imputation/HRC/P12/chr1p36.1 --cojo-file Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED_flipped_1p36.1_alleles_no_positions.txt --extract chr1p36.1_snplist.txt --cojo-cond rs6692586+rs4654828+rs12408468.txt --cojo-p 5e-2 --out Meta-analysis_Locke_et_al_1p36.1_COJO_conditioned_on_rs6692586+rs4654828+rs12408468

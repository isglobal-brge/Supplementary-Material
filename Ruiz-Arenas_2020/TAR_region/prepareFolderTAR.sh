#'#################################################################################
#'#################################################################################
#' Prepare folder
#'#################################################################################
#'#################################################################################

mkdir data
mkdir results
mkdir results/preproc
mkdir results/models
mkdir results/figures

preproc=results/preproc

## Copy to data the imputed and phased files
cp ~/PublicData/STUDY/1000GENOME/Samples_Pop1GK.Rdata data
ln -s ~/PublicData/STUDY/1000GENOME/VCF/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz data/1000G_chr1.vcf.gz
ln -s ~/PublicData/STUDY/1000GENOME/VCF/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi data/1000G_chr1.vcf.gz.tbi

## GTEx data
cp  ~/data/CarlosRuiz/Inversions/LDclassifier/CancerRegion/GTEX/GTEX.het_check.csv data/GTEX.het_check.csv
ln -s ~/data/CarlosRuiz/Inversions/GTEX/GTEX_annotated.vcf.gz data/GTEX.vcf.gz
ln -s ~/data/CarlosRuiz/Inversions/GTEX/GTEX_RSE_covariates_inversions_goodTissues.RData data/GTEX_SE.Rdata
cp  ~/data/CarlosRuiz/Inversions/GTEX/GTEx_Analysis_v7_eQTL_covariates/Whole_Blood.v7.covariates.txt data/whole_blood_covars.txt
cp  ~/data/CarlosRuiz/Inversions/GTEX/GTEx_Analysis_v7_eQTL_covariates/Cells_EBV-transformed_lymphocytes.v7.covariates.txt data/lymphos_covars.txt


ln -s ~/data/CarlosRuiz/Inversions/LDclassifier/CancerRegion/ImputationGWAS/chr1.dose.idsrs.vcf.gz data/CGEMS.vcf.gz
ln -s ~/data/CarlosRuiz/Inversions/LDclassifier/CancerRegion/ImputationGWAS/chr1.dose.idsrs.vcf.gz.tbi data/CGEMS.vcf.gz.tbi

cp ~/data/CarlosRuiz/Inversions/Cancer/GWAS_Breast/final.fam data/cgems.fam
cp ~/data/CarlosRuiz/Inversions/Cancer/GWAS_Breast/phs000147.v3.pht000663.v3.p1.c1.CGEMS_Breast_Cancer_Subject_Phenotypes.GRU.txt data/cgems.phenos.txt
cp ~/data/CarlosRuiz/Inversions/LDclassifier/CancerRegion/imputed.adjust.assoc.logistic data/cgems.gwas_res.txt



ln -s ~/data/CarlosRuiz/Inversions/LDclassifier/CancerRegion/ImputationTCGA/chr1.dose.idsrs.vcf.gz data/TCGA.vcf.gz
ln -s ~/data/CarlosRuiz/Inversions/LDclassifier/CancerRegion/ImputationTCGA/chr1.dose.idsrs.vcf.gz.tbi data/TCGA.vcf.gz.tbi
cp ~/data/CarlosRuiz/Inversions/LDclassifier/CancerRegion/ImputationTCGA/TCGA_ancestry.Rdata data/


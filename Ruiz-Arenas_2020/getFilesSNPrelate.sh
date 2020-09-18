#'#################################################################################
#'#################################################################################
#' Prepare folder for runSNPrelate.Rmd
#'#################################################################################
#'#################################################################################

mkdir data
mkdir results
mkdir results/preproc
mkdir results/plotFiles

preproc=results/preproc

## Add Drosophila raw data (downloaded from http://dgrp2.gnets.ncsu.edu/data.html)
ln -s ~/data/CarlosRuiz/RecombClust/Drosophila/data/dgrp2.vcf.gz data/dgrp2.vcf.gz
ln -s ~/data/CarlosRuiz/RecombClust/Drosophila/data/dgrp2.vcf.gz.tbi data/dgrp2.vcf.gz.tbi
ln -s ~/data/CarlosRuiz/RecombClust/Drosophila/data/ranges.csv data/ranges.csv
ln -s ~/data/CarlosRuiz/RecombClust/Drosophila/data/inv_genos.xlsx data/inv_genos.xlsx

## Add Human raw data (Data from 1000 Genomes)
ln -s ~/PublicData/STUDY/1000GENOME/VCF/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz data/chr1.vcf.gz
ln -s ~/PublicData/STUDY/1000GENOME/VCF/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi data/chr1.vcf.gz.tbi
ln -s ~/PublicData/STUDY/1000GENOME/VCF/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz data/chr2.vcf.gz
ln -s ~/PublicData/STUDY/1000GENOME/VCF/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi data/chr2.vcf.gz.tbi
ln -s ~/PublicData/STUDY/1000GENOME/VCF/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz data/chr8.vcf.gz
ln -s ~/PublicData/STUDY/1000GENOME/VCF/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi data/chr8.vcf.gz.tbi
ln -s ~/PublicData/STUDY/1000GENOME/VCF/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz data/chr17.vcf.gz
ln -s ~/PublicData/STUDY/1000GENOME/VCF/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi data/chr17.vcf.gz.tbi
ln -s ~/data/CarlosRuiz/Inversions/Samples_Pop1GK.Rdata data/Samples_Pop1GK.Rdata

## Add recombClust classifications
ln -s ~/data/CarlosRuiz/RecombClust/HumanPopScan/results/models/reg24recombRes.Rdata data/LCTmodels.Rdata
ln -s ~/data/CarlosRuiz/RecombClust/Chr1_region/results/models/recombClust_1000G_short.Rdata data/TAR_reg_models.Rdata
ln -s ~/data/CarlosRuiz/Inversions/LDclassifier/EuropeanSamples/PCsinvs8_17_EUR_pruned0.4.Rdata data/HumanInvsPCs.Rdata

## recombClust results obtained with runRecombClustDrosophila.Rmd
ln -s ~/data/CarlosRuiz/RecombClust/Drosophila/results/models/In\(2L\)trecombRes.Rdata data/Dros2Lmodels.Rdata
ln -s ~/data/CarlosRuiz/RecombClust/Drosophila/results/models/In\(2R\)NSrecombRes.Rdata data/Dros2Rmodels.Rdata
ln -s ~/data/CarlosRuiz/RecombClust/Drosophila/results/models/In\(3R\)MorecombRes.Rdata data/Dros3Rmodels.Rdata


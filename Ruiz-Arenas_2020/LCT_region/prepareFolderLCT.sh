#'#################################################################################
#'#################################################################################
#' Prepare folder for LCT analysis
#'#################################################################################
#'#################################################################################

mkdir data
mkdir data/recomb
mkdir results
mkdir results/models

## Add individuals ancestry
cp ~/PublicData/STUDY/1000GENOME/Samples_Pop1GK.Rdata data/

## Download recombination maps (in data)
wget https://pophuman.uab.cat/files/wig/recomb_Bherer2017_sexavg_10kb.bw

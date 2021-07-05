---
title: "a Bioconductor package to identify outliers in rare diseases DNA methylation data"
author: "Leire Abarrategui, Carlos Ruiz-Arenas, Carles Hernandez-Ferrer, Juan R. Gonzalez, Patricia Ryser-Welch"
date: "2021-07-05"
header-includes:
   - \usepackage{float}
output:
  pdf_document: 
    latex_engine: xelatex
    
subtitle: "Supplementary material" 
vignette: >
  %\VignetteIndexEntry{Detection of epimutations with state of the art methods in methylation data} 
  %\VignetteEngine{knitr::rmarkdown} 
  %\VignetteEncoding{UTF-8}
bibliography: refs.bib
---



# Introduction

## Background

Rare diseases are pathologies with a low prevalence (< 1 per 2,000 people) [@EU_RD]. Most of these pathologies have an onset during childhood and a strong genetic etiology [@lopez2016social]. Consequently, rare disease diagnosis has relied on identifying genetic and genomic mutations that can cause the disease [@aref2019]. Although these variants have provided a diagnosis for many patients and families, around 60% of the cases remained undiagnosed [@lionel2018improved]. Aberrant methylation can be an underlying cause of undiagnosed patients, either as a primary event (a.k.a. epimutation) or as a functional consequence of chromatin dysregulation by genetic or environmental agents (a.k.a. episignature)  [@aref2019]. Epimutations are the cause of some rare diseases, such as Prader-Willi, Angelman or Beckwith-Wiedemann syndromes [@aref2019] and some human malformations [@serra2015dna]. Syndrome-specific episignatures are increasingly defined as biomarkers for a growing number of disorders [@aref2019; @garg2020survey]. Therefore, tools to detect epimutations and episignatures should be made available to the rare disease community and included in standardized analysis workflows.


This manual describes the `epimutacions` package tools to identify epivariants using multiple outlier detection  approaches. Also,  includes functions to plot and annotate the epimutations. The full `epimutacions ` user´s guide is available in this vignette. 

The name of the package is `epimutacions` (pronounced `ɛ pi mu ta 'sj ons`) which means epimutations in Catalan, a language from the northeast of Spain.

## Methodology

The `epimutacions` package computes a genome-wide DNA methylation analysis to detect the epigenetic variants to be considered as biomarkers for samples with rare diseases (epimutations). The method compares a case sample with suspected rare disease against a reference panel.  The package focused on the detection of outliers in DNA methylation patterns associated with the diseases as proposed by [@aref2019].  

The identification of relevant genomic methylation regions for a given sample having a rare disease will be driven by detecting differentially methylated CpG sites when comparing beta values of all control samples with the given proband. Firstly, bump-hunter [@jaffe2012bump] approach is used to identify the Differentially Methylated Regions (DMRs). After that, CpGs in the proband sample are tested in those DMRs in order to identify regions with CpGs being outliers when comparing with the reference panel. To this end, different anomaly detection statistical approaches are used. These include Multivariate Analysis of Variance (MANOVA) [@friedrich2017package], Multivariate Linear Model [@martin2020multivariate], isolation forest [@cortes2021package] and robust mahalanobis distance [@maechler2021package]. However, Barbosa [@barbosa2018] and Beta methods do not use bump-hunter output. Barbosa [@barbosa2018] checks for each CpG, if the proband’s measurement is an outlier. Then, it calls an epimutation to those regions where 3 contiguous CpGs are outliers, and they are separated by less than 500 base pairs. Beta approach models the DNA methylation data using a beta distribution.  

## Input data

The package allows two different types of inputs: 

 * (1) Case samples `IDAT` files (raw microarray intensities) together with `RGChannelSet` class object as reference panel. The reference panel can be supplied by the user or can be selected through the example datasets that the package provides. 
 
* (2) `GenomicRatioSet` class object containing case and control samples. 

The input data should contain information about \beta values of CpG sites, phenotype and feature data. 

Normalization through `epi_preprocess()` function is highly recommended when  combining data from different sources.  In order to remove the unwanted variation caused by the batch effect when combining data from different sources.   

Finally, a `GenomicRatioSet` class object  the input of the main function,   `epimutations()` function. It should be mentioned that  case samples and reference panel are introduced separately. 

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021/vignette/fig/workflow} 

}

\caption{Allowed data formats, normalization and input types}\label{fig:workflow}
\end{figure}



# Getting started

The `epimutacions` package is installed by executing: 

```r
install_github("isglobal-brge/epimutacions")
```

The package is it loaded in R as usual:


```r
library(epimutacions)
```

The document has the following dependencies: 


```r
library(Knitr)
library(kableExtra)
```

# Datasets

## Candidate regions

Epimutations detection has two main steps: (1) definition of candidate regions and (2) evaluation of outlier significance. Although there are different algorithms to define epimutations regions, they share common features. In general, we define an epimutation as at least 3 contiguous  CpGs with a maximum distance of 1kb between them.  

In Illumina 450K array, probes are unequally distributed along the genome, limiting the number of regions that can fulfil the requirements to be considered an epimutation. So, we have computed a dataset containing the regions that are candidates to become an epimutation.  

To define the candidate epimutations, we relied on the clustering from bumphunter. We defined a primary dataset with all the CpGs from the Illumina 450K array.  Then, we run bumphunter and selected those regions with at least 3 CpGs. As a result, we found 40408 candidate epimutations which are available in `candRegsGR` dataset. The code for generating these regions can be found in epimutacion package. 

In addition, we converted the candidate region from hg19 to hg38 coordinates, using NCBI remap. We selected regions that mapped to one region in hg38 with the same length. This yielded a total of 39944, the 98.85% of total hg19 regions. After converting to hg38, we can use these ranges to be annotated to ENCODE cREs. Overall, we mapped 30163 candidate regions to cREs, representing 74.65% of total candidate regions. 



```r
data("candRegsGR")
candRegsGR
```

```
GRanges object with 40408 ranges and 9 metadata columns:
                 seqnames              ranges strand |     value      area
                    <Rle>           <IRanges>  <Rle> | <numeric> <numeric>
   chr6_32128101     chr6   32128101-32173532      * |         1       381
   chr6_33156164     chr6   33156164-33181870      * |         1       291
   chr6_32034322     chr6   32034322-32059605      * |         1       239
   chr6_31618987     chr6   31618987-31639143      * |         1       234
   chr6_33279563     chr6   33279563-33292029      * |         1       233
             ...      ...                 ...    ... .       ...       ...
  chr9_140652685     chr9 140652685-140652743      * |         1         3
  chr9_140656200     chr9 140656200-140657381      * |         1         3
  chr9_140680393     chr9 140680393-140681206      * |         1         3
  chr9_140732731     chr9 140732731-140733980      * |         1         3
  chr9_141012312     chr9 141012312-141013537      * |         1         3
                   cluster indexStart  indexEnd         L  clusterL
                 <numeric>  <integer> <integer> <numeric> <integer>
   chr6_32128101    133070     165174    165554       381       381
   chr6_33156164    133204     167451    167741       291       291
   chr6_32034322    133058     164512    164750       239       239
   chr6_31618987    132987     162583    162816       234       234
   chr6_33279563    133221     168282    168514       233       233
             ...       ...        ...       ...       ...       ...
  chr9_140652685    162642     247198    247200         3         3
  chr9_140656200    162643     247201    247203         3         3
  chr9_140680393    162649     247214    247216         3         3
  chr9_140732731    162660     247252    247254         3         3
  chr9_141012312    162683     247294    247296         3         3
                                    CRE               CRE_type
                            <character>            <character>
   chr6_32128101 EH38E2459822,EH38E24.. pELS,CTCF-bound;PLS;..
   chr6_33156164 EH38E2460436,EH38E24.. PLS;pELS,CTCF-bound;..
   chr6_32034322 EH38E2459711,EH38E24.. dELS;dELS;dELS;pELS,..
   chr6_31618987 EH38E2459340,EH38E24.. pELS,CTCF-bound;PLS,..
   chr6_33279563 EH38E2460551,EH38E24.. pELS,CTCF-bound;pELS..
             ...                    ...                    ...
  chr9_140652685                                              
  chr9_140656200 EH38E2738315,EH38E27.. pELS,CTCF-bound;pELS..
  chr9_140680393 EH38E2738332,EH38E27..   dELS,CTCF-bound;dELS
  chr9_140732731                                              
  chr9_141012312                                              
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
```

## GenomicRatioSet

The package includes a small `GenomicRatioSet` class dataset (`methy` ) containing the DNA methylation profiles from a total of   individuals, 3 cases and 48 controls. The DNA methylation profiles were generated using the Illumina 450k Human Methylation BeadChip. The data were obtained from [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) and adapted for the package usage. 



```r
data("methy")
methy
```

```
class: GenomicRatioSet 
dim: 80731 51 
metadata(0):
assays(3): Beta M CN
rownames(80731): cg00725145 cg16080333 ... cg07468397 cg08821909
rowData names(0):
colnames(51): GSM2808239 GSM2808240 ... GSM2562700 GSM2562701
colData names(4): sampleID age sex status
Annotation
  array: IlluminaHumanMethylation450k
  annotation: ilmn12.hg19
Preprocessing
  Method: NA
  minfi version: NA
  Manifest version: NA
```

```r
table(methy$status)
```

```

   case control 
      3      48 
```

We are going to create two different datasets for further analysis, `case_samples` and `control_panel`: 

```r
case_samples <- methy[,methy$status == "case"]
control_samples <- methy[,methy$status == "control"]
```

# Preprocessing
The preprocessing in `epimutacions` package is done by `epi_preprocess()` function. It contains 6 preprocessing methods corresponding to minfi package that can be selected by the user:

+---------------+-----------------------------+----------------------------------------------------------------------+
| Method        |  Function                   | Description                                                          |
+===============+=============================+======================================================================+
| `raw`         | `preprocessRaw`             | Converts the Red/Green channel for an Illumina methylation           |  
|               |                             |  array into methylation signal, without using any normalization      |
+---------------+-----------------------------+----------------------------------------------------------------------+
| `illumina  `  | `preprocessIllumina`        |  Implements preprocessing for Illumina methylation                   |
|               |                             |      microarrays as used in Genome Studio                            |
+---------------+-----------------------------+----------------------------------------------------------------------+
|`swan`         | `preprocessSWAN `           |Subset-quantile Within Array Normalisation (SWAN). It allows Infinium |
|               |                             | I and II type probes on a single array to be normalized together     |
+---------------+-----------------------------+----------------------------------------------------------------------+
|  `quantile`   |`preprocessQuantile`         |Implements stratified quantile normalization preprocessing for        |                
|               |                             |Illumina methylation microarrays                                      |
+---------------+-----------------------------+----------------------------------------------------------------------+
| `noob`        | `preprocessNoob`            | Noob (normal-exponential out-of-band) is a background correction     | 
|               |                             | method with dye-bias normalization for                               |
|               |                             | Illumina Infinium methylation arrays                                 |
+---------------+-----------------------------+----------------------------------------------------------------------+
| `funnorm`     | `preprocessFunnorm `        | Functional normalization (FunNorm) is a between-array                |
|               |                             | normalization method for the Illumina Infinium                       |
|               |                             | HumanMethylation450 platform                                         | 
+---------------+-----------------------------+----------------------------------------------------------------------+

In addition, the unique parameters for each normalization approach are defined through `norm_parameters()`: 

+----------+-----------------------+-----------------------------------------------------------------------------+
| Method   |Parameters             |Description                                                                  |
+==========+=======================+=============================================================================+
|`illumina`|`bg.correct`\          |Performs background correction\                                              |
|          |`normalize`\           |Performs controls normalization\                                             |        
|          |`reference`            |The reference array for control normalization                                |
+----------+-----------------------+-----------------------------------------------------------------------------+
|`quantile`|`fixOutliers`\         | Low outlier Meth and Unmeth signals will be fixed\                          |                       
|          |`removeBadSamples`\    | Remove bad samples\                                                         |
|          |`badSampleCutoff`\     | The cutoff to label samples as 'bad'\                                       |
|          |`quantileNormalize`\   | Performs quantile normalization\                                            |
|          | `stratified`          | Performs quantile normalization within region strata\                       |
|          |`mergeManifest`\       | Merged to the output the information in the associated manifest package\    |
|          | \                     |                                                                             |
|          |`sex`                  | Sex of the samples                                                          |
+----------+-----------------------+-----------------------------------------------------------------------------+
|`noob`    | `offset`\             | Offset for the normexp background correct\                                  |
|          | `dyeCorr`\            | Performs dye normalization\                                                 |  
|          | `dyeMethod`           | Dye bias correction to be done                                              |
+---------------+------------------+-----------------------------------------------------------------------------+
|`funnorm` | `nPCs `\              | The number of principal components from the control probes\                 |       
|          | `sex`\                | Sex of the samples\                                                         |      
|          | `bgCorr`\             | Performs NOOB background correction prior to functional normalization\      |
|          | \                     |                                                                             |
|          | `dyeCorr`\            | Performs dye normalization\                                                 |
|          | `keepCN`              | Keeps copy number estimates                                                 |
+---------------+-----------------------------+------------------------------------------------------------------+

The default settings for each method can be obtained by invoking  the function `norm_parameters()` with no arguments:


```r
norm_parameters()
```

```
$illumina
$illumina$bg.correct
[1] TRUE

$illumina$normalize
[1] "controls" "no"      

$illumina$reference
[1] 1


$quantile
$quantile$fixOutliers
[1] TRUE

$quantile$removeBadSamples
[1] FALSE

$quantile$badSampleCutoff
[1] 10.5

$quantile$quantileNormalize
[1] TRUE

$quantile$stratified
[1] TRUE

$quantile$mergeManifest
[1] FALSE

$quantile$sex
NULL


$noob
$noob$offset
[1] 15

$noob$dyeCorr
[1] TRUE

$noob$dyeMethod
[1] "single"    "reference"


$funnorm
$funnorm$nPCs
[1] 2

$funnorm$sex
NULL

$funnorm$bgCorr
[1] TRUE

$funnorm$dyeCorr
[1] TRUE

$funnorm$keepCN
[1] FALSE
```

However, to modify the parameters related to a  method you can do as the following example for `illumina` approach:


```r
parameters <- norm_parameters(illumina = list("bg.correct" = FALSE))
parameters$illumina$bg.correct
```

```
[1] FALSE
```

# Epimutations

## Epimutations detection

The `epimutacions` package includes 6 methods for epivariants identification: (1) Multivariate Analysis of variance (`manova`), (2) Multivariate Linear Model (`mlm`), (3) isolation forest (`isoforest`), (4) robust mahalanobis distance  (`mahdistmcd`) (5) `barbosa` and (6) `beta`. 


In the mentioned first 4 methods, firstly, Differentially Methylated Regions (DMRs) are identified using bump-hunter method [@jaffe2012bump]. Then, those DMRs are tested to identify regions with CpGs being outliers when comparing with the reference panel. However, `barbosa` and `beta` do not identify outliers by filtering the DMRs. `barbosa`  utilized a sliding window approach to individually compare the methylation value  in each proband against the reference panel. `Beta` used beta distribution to identify epivariants in the case sample. 



```r
epi_mvo <- epimutations(case_samples, control_samples, method = "manova")
epi_ml <- epimutations(case_samples, control_samples, method = "mlm")
epi_iso <- epimutations(case_samples, control_samples, method = "isoforest")
epi_mcd <- epimutations(case_samples, control_samples, method = "mahdistmcd")
```



```r
epi_brb <- epimutations(case_samples, control_samples, method = "barbosa")
epi_beta <- epimutations(case_samples, control_samples, method = "beta")
```


## Unique parameters

The `epi_parameters()` function  is  useful to set the individual parameters for each   approach. The arguments are described in the following table: 

+---------------+-----------------------------+----------------------------------------------------------------------+
| Method        | Parameter                   | Description                                                          |
+===============+=============================+======================================================================+
| `manova`\     | `pvalue_cutoff`             | The threshold p-value to select which CpG regions are outliers       |
|  `mlm`\       |                             |                                                                      |
|  `beta`       |                             |                                                                      |
+---------------+-----------------------------+----------------------------------------------------------------------+
| `iso.forest`  | `outlier_score_cutoff`\     | The threshold to select which CpG regions are outliers\              |
|               | `ntrees`                    | The number of binary trees to build for the model                    |
+---------------+-----------------------------+----------------------------------------------------------------------+
|`mahdist.mcd`  | `nsamp`                     | The number of subsets used for initial estimates in the MCD          |
+---------------+-----------------------------+----------------------------------------------------------------------+
|  `barbosa`    |`window_sz`\                 |The maximum distance between CpGs to be considered in the same DMR\   |                
|               | \                           |                                                                      |
|               |`offset_mean`/`offset_abs`   | The upper and lower threshold to consider a CpG an outlier           |
+---------------+-----------------------------+----------------------------------------------------------------------+
| `beta`        | `pvalue_cutoff`\            | The minimum p-value to consider a CpG an outlier\                    |
|               | `diff_threshold`            | The minimum methylation difference between the CpG and the           |
|               |                             |  mean methylation to consider a position an outlier                  |
+---------------+-----------------------------+----------------------------------------------------------------------+

Invoking `epi_parameters()` with no arguments returns a list of the default settings for each method: 


```r
epi_parameters()
```

```
$manova
$manova$pvalue_cutoff
[1] 0.05


$mlm
$mlm$pvalue_cutoff
[1] 0.05


$isoforest
$isoforest$outlier_score_cutoff
[1] 0.5

$isoforest$ntrees
[1] 100


$mahdistmcd
$mahdistmcd$nsamp
[1] "deterministic"


$barbosa
$barbosa$window_sz
[1] 10

$barbosa$offset_mean
[1] 0.15

$barbosa$offset_abs
[1] 0.1


$beta
$beta$pvalue_cutoff
[1] 1e-06

$beta$diff_threshold
[1] 0.1
```

The set up of any parameter can be done as the following example of p-value cut-off for  `manova`: 


```r
parameters <- epi_parameters(manova = list("pvalue_cutoff" = 0.01))
parameters$manova$pvalue_cutoff
```

```
[1] 0.01
```


## Results description

The `epimutations` function returns a tibble containing all the epivariants identified in the given case sample. In case no epimutation is found, a row containing the case sample information and missing values for each argument is returned. The following table describes each argument in the result data frame: 


+-----------------------+-----------------------------------------------------------------------------------------------------------------------+
| Column name           | Description                                                                                                           |
+=======================+=======================================================================================================================+
| `epi_id`              | Systematic name for each epimutation identified                                                                       |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------+
| `sample`              | The name of the sample containing that epimutation                                                                    |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------+
| `chromosome`          | The location of the epimutation                                                                                       |
|`start`                |                                                                                                                       |
|`end`                  |                                                                                                                       |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------+
| `sz`                  | The window's size of the event                                                                                        |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------+
| `cpg_n`               | The number of CpGs in the epimutation                                                                                 |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------+
| `cpg_n`               | The names of CpGs in the epimutation                                                                                  |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------+
| `outlier_score`       | For method `manova` it provides the approximation to F-test and the Pillai score, separated by `/`\                   |
|                       | For method `mlm` it provides the approximation to F-test and the R2 of the model, separated by `/`\                   |
|                       | For method `isoforest` it provides the magnitude of the outlier score.\                                               |
|                       | For method `beta` it provides the mean p-value of all GpGs in that DMR\                                               |
|                       | For methods `barbosa` and `mahdistmcd` it is filled with `NA`.                                                        |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------+
| `pvalue`              | For methods `manova` and  `mlm` it provides the p-value obtained from the model.\                                     |
|                       | For method `barbosa`, `isoforest`, `beta` and `mahdistmcd` it is filled with `NA`.                                    |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------+ 
|`outlier_direction`    |    Indicates the direction of the outlier with "hypomethylation" and "hypermethylation".\                             | 
|                       | For `manova`, `mlm`, `isoforest`, and `mahdistmcd` it is computed from the values obtained from `bumphunter`.\        |
|                       | For `beta` is computed from the p value for each CpG using `diff_threshold` and `pvalue_threshold` arguments.\        |
|                       | For `barbosa` it is computed from the location of the sample in the reference distribution (left vs. right outlier).  |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------+
| `adj_pvalue`          | For methods `manova` and  `mlm` it provides the adjusted p-value with                                                 |
|                       | Benjamini-Hochberg based on the total number of regions detected by Bumphunter.     \                                 |
|                       | For method `barbosa`, `isoforest`, `mahdistmcd` and `beta` it is filled with `NA`.                                    |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------+ 
| `epi_region_id`       | Name of the epimutation region as defined in `candRegsGR`.                                                            |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------+
| `CRE`                 | cREs (cis-Regulatory Elements) as defined by ENCODE overlapping the epimutation region.                               |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------+
| `CRE_type`            | Type of cREs (cis-Regulatory Elements) as defined by ENCODE.                                                          |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------+

## Epimutations annotations

The `epimutacions` package also includes the `annotate_epimutations` function dedicated to enriching the epimutations identified by the previously described methods:


```r
rst_mvo <- annotate_epimutations(epi_mvo)
```


```r
rst_mvo[1:2, c(1, 12:14)]
```


Table: epimutations annotation

|epi_id       | adj_pvalue|epi_region_id  |CRE                                                                                                                                            |
|:------------|----------:|:--------------|:----------------------------------------------------------------------------------------------------------------------------------------------|
|epi_manova_1 |          0|chr19_12776725 |EH38E1939817,EH38E1939818,EH38E1939819                                                                                                         |
|epi_manova_2 |          0|chr7_90892836  |EH38E2570884,EH38E2570885,EH38E2570886,EH38E2570887,EH38E2570888,EH38E2570889,EH38E2570890,EH38E2570891,EH38E2570892,EH38E2570893,EH38E2570894 |


## Epimutation visualization

The  visualization approach locates the epimutations along the genome. The function `plot_epimutations` plots the methylation values of the  individual with the epimutation in red, the control samples in dashed black lines and population mean in blue:  

```r
plot_epimutations(as.data.frame(epi_mvo[1,]), methy)
```

![](sup_mat_files/figure-latex/unnamed-chunk-12-1.pdf)<!-- --> 

Furthermore, it includes the gene annotations in the regions in which the epivariation is located.  This can be achieved by using the argument  `gene_annot == TRUE`: 


```r
plot_epimutations(as.data.frame(epi_mvo[1,]), methy, genes_annot = TRUE)
```

![](sup_mat_files/figure-latex/plot_mvo_genes_annot-1.pdf)<!-- --> 


Also, it is possible to plot the chromatin marks H3K4me3, H3K27me3  and H3K27ac by setting the argument `regulation = TRUE`:

* **H3K4me3**: commonly associated with the activation of transcription of nearby genes.
* **H3K27me3**: is used in epigenetics to look for inactive genes.
* **H3K27ac**: is associated with the higher activation of transcription and therefore defined as an active enhancer mark



```r
plot_epimutations(as.data.frame(epi_mvo[1,]), methy, regulation = TRUE)
```

![](sup_mat_files/figure-latex/plot_mvo_regulation-1.pdf)<!-- --> 

# Method validation

## Data collection 

The data were obtained for the studies previously described [@garg2020survey].  The datasets were downloaded from [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/).  We accessed DNA methylation data from a total 1, 417 individuals from [GSE51032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51032) and 
[GSE111629](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) cohorts. The DNA methylation profiles were generated using the Illumina 450k Human Methylation BeadChip.  

The GSE51032 study analysed primary cancers samples:  424 cancer free, 235 primary breast cancer, 166 
primary colorectal cancer and 20 other primary cancers. The GSE111629 cohort 335 Parkinson’s disease and 237 control samples.  

## Validation

We evaluated the performance of the method using TPR (True Positive Rate), False Positive Rate (FPR) and accuracy. We use the TPR to measure the proportion of detected epivariations by the `epimutations` approach  present in the validated [@garg2020survey]. FPR to calculate the identified epimutations outside the once found in [@garg2020survey], whether validated or not. The accuracy measures the closeness of the detected epimutation to the validated regions. 

We select samples differently depending on the study group and measure to compute. Control samples were selected randomly using different sample size: 20, 30, 40, 50, 60, 70, 80, 90 and 100. However, case samples were selected considering validated epimutations (for TPR and accuracy)  or excluding epivariations found (for FPR) [@garg2020survey].  

The validated epimutations in table 1 were only present on 5 individuals: GSM1235784 from GSE51032 cohort and GSM3035933, GSM3035791, GSM3035807 and GSM3035685 from GSE111629. Therefore, they were established as case samples when computing TPR and accuracy. Nevertheless, we compute FPR excluding the samples containing at least one epimutation found by [@garg2020survey].  For the remaining case samples, 4 were selected randomly in each execution.

We execute 100 times the same process for each control sample size. We define for the analysis regions of $\approx$ 20 kb containing $\ge$ 3 GpGs. 


Table: validated epimutations (Garg et al. 2020).

|Chromosome |    Start|      End| Width|Strand |Samples               |
|:----------|--------:|--------:|-----:|:------|:---------------------|
|chr17      | 46018653| 46019185|   533|*      |GSM1235784/GSM3035791 |
|chr19      | 11199850| 11200147|   298|*      |GSM3035685            |
|chr5       | 10249760| 10251253|  1494|*      |GSM3035933            |
|chr5       | 67583971| 67584381|   411|*      |GSM3035791/GSM3035807 |

Additionally, we have plotted the methylation values of the samples in the regions where the validated epimutations were found.
 

```
[1] "C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021"
```

\begin{figure}[H]

{\centering \includegraphics[width=28.82in]{C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021/vignette/fig/GSM1235784_chr17-46018654-46019184} 

}

\caption{GSE51032 samples in the region chr17:46018654-46019184}\label{fig:graph_GSM1235784}
\end{figure}

\begin{figure}[H]

{\centering \includegraphics[width=28.82in]{C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021/vignette/fig/GSM3035685_chr19-11199851-11200146} 

}

\caption{GSE111629 samples in the region chr19:11199851-11200146}\label{fig:graph_GSM3035685}
\end{figure}

\begin{figure}[H]

{\centering \includegraphics[width=28.82in]{C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021/vignette/fig/GSM3035791_chr5-67584194-67584380} 

}

\caption{GSE111629 samples in the region chr5:67584194-67584380}\label{fig:graph_GSM3035791}
\end{figure}

\begin{figure}[H]

{\centering \includegraphics[width=28.82in]{C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021/vignette/fig/GSM3035791_chr17-46018654-46019184} 

}

\caption{GSE111629 samples in the region chr17:46018654-46019184}\label{fig:graph_GSM3035791.2}
\end{figure}

\begin{figure}[H]

{\centering \includegraphics[width=28.82in]{C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021/vignette/fig/GSM3035807_chr5-67583972-67584380} 

}

\caption{GSE111629 samples in the region chr5:67583972-67584380}\label{fig:graph_GSM3035807}
\end{figure}

\begin{figure}[H]

{\centering \includegraphics[width=28.82in]{C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021/vignette/fig/GSM3035933_chr5-10249761-10251252} 

}

\caption{GSE111629 samples in the region chr5:10249761-10251252}\label{fig:graph_GSM3035933}
\end{figure}
## Results

We compare GSM1235784 case sample against randomly selected control samples from GSE51032 and GSM3035933, GSM3035791, GSM3035807 and GSM3035685 case samples against controls from GSE111629 specifying a region of 20 kb and $\ge$ 3 GpGs. 

We obtained similar results in both cohorts. We observed that the methods manova, mahalanobis distance and multivariate linear models identified the validated epimutations with a TPR of $>99\%$ even if the control sample is small. However,  the TPR in isolation forest increases together with the number of control samples obtaining a TPR $\ge$ 75 with 50 control samples or more. The TPR in barbosa and beta approaches for GSE51032 dataset is small ($<50\%$).  Nonetheless, for GSE111629 the TPR value increases considerably $>99\%$. Regarding the accuracy, all the statistical  approaches detect the epivariants with $>80\%$ of closeness to the validated epimutations. 

We detected possible epivariations outside the epimutations found by [@garg2020survey] selecting control and case samples randomly. For the analysis, we selected  regions of 20 kb and $\ge$ 3 GpGs. We compared each case sample  individually against control samples.  We observed that in both cohorts and for every approach the FPR value is very small $<0.01\%$. 

\begin{table}
\centering
\begin{tabular}[t]{c|c|c|c}
\hline
method & TPR & accuracy & FPR\\
\hline
\multicolumn{4}{l}{\textbf{n20}}\\
\hline
\hspace{1em}manova & 98 & 99.6 & 0.00\\
\hline
\hspace{1em}mlm & 100 & 99.6 & \vphantom{8} 0.00\\
\hline
\hspace{1em}mahdistmcd & 100 & 99.6 & \vphantom{8} 0.00\\
\hline
\hspace{1em}isoforest & 28 & 99.6 & 0.00\\
\hline
\hspace{1em}barbosa & 100 & 94.4 & 0.25\\
\hline
\hspace{1em}beta & 99 & 89.0 & 0.25\\
\hline
\multicolumn{4}{l}{\textbf{n30}}\\
\hline
\hspace{1em}manova & 100 & 99.6 & \vphantom{7} 0.00\\
\hline
\hspace{1em}mlm & 100 & 99.6 & \vphantom{7} 0.00\\
\hline
\hspace{1em}mahdistmcd & 100 & 99.6 & \vphantom{7} 0.00\\
\hline
\hspace{1em}isoforest & 63 & 99.6 & 0.00\\
\hline
\hspace{1em}barbosa & 100 & 94.4 & 0.00\\
\hline
\hspace{1em}beta & 100 & 95.4 & 0.00\\
\hline
\multicolumn{4}{l}{\textbf{n40}}\\
\hline
\hspace{1em}manova & 100 & 99.6 & \vphantom{6} 0.00\\
\hline
\hspace{1em}mlm & 100 & 99.6 & \vphantom{6} 0.00\\
\hline
\hspace{1em}mahdistmcd & 100 & 99.6 & \vphantom{6} 0.00\\
\hline
\hspace{1em}isoforest & 65 & 99.6 & 0.00\\
\hline
\hspace{1em}barbosa & 100 & 92.4 & 0.00\\
\hline
\hspace{1em}beta & 100 & 92.6 & \vphantom{1} 0.00\\
\hline
\multicolumn{4}{l}{\textbf{n50}}\\
\hline
\hspace{1em}manova & 100 & 99.6 & \vphantom{5} 0.00\\
\hline
\hspace{1em}mlm & 100 & 99.6 & \vphantom{5} 0.00\\
\hline
\hspace{1em}mahdistmcd & 100 & 99.6 & \vphantom{5} 0.00\\
\hline
\hspace{1em}isoforest & 83 & 99.6 & 0.00\\
\hline
\hspace{1em}barbosa & 100 & 88.0 & 0.00\\
\hline
\hspace{1em}beta & 100 & 93.8 & 0.00\\
\hline
\multicolumn{4}{l}{\textbf{n60}}\\
\hline
\hspace{1em}manova & 100 & 99.6 & \vphantom{4} 0.00\\
\hline
\hspace{1em}mlm & 100 & 99.6 & \vphantom{4} 0.00\\
\hline
\hspace{1em}mahdistmcd & 100 & 99.6 & \vphantom{4} 0.00\\
\hline
\hspace{1em}isoforest & 88 & 99.6 & 0.00\\
\hline
\hspace{1em}barbosa & 100 & 89.2 & 0.00\\
\hline
\hspace{1em}beta & 100 & 92.6 & 0.00\\
\hline
\multicolumn{4}{l}{\textbf{n70}}\\
\hline
\hspace{1em}manova & 100 & 99.6 & \vphantom{3} 0.00\\
\hline
\hspace{1em}mlm & 100 & 99.6 & \vphantom{3} 0.00\\
\hline
\hspace{1em}mahdistmcd & 100 & 99.6 & \vphantom{3} 0.00\\
\hline
\hspace{1em}isoforest & 94 & 99.6 & 0.00\\
\hline
\hspace{1em}barbosa & 100 & 90.0 & 0.00\\
\hline
\hspace{1em}beta & 100 & 93.2 & 0.00\\
\hline
\multicolumn{4}{l}{\textbf{n80}}\\
\hline
\hspace{1em}manova & 100 & 99.6 & \vphantom{2} 0.00\\
\hline
\hspace{1em}mlm & 100 & 99.6 & \vphantom{2} 0.00\\
\hline
\hspace{1em}mahdistmcd & 100 & 99.6 & \vphantom{2} 0.00\\
\hline
\hspace{1em}isoforest & 97 & 99.6 & 0.00\\
\hline
\hspace{1em}barbosa & 100 & 86.0 & 0.00\\
\hline
\hspace{1em}beta & 100 & 92.9 & \vphantom{1} 0.00\\
\hline
\multicolumn{4}{l}{\textbf{n90}}\\
\hline
\hspace{1em}manova & 100 & 99.6 & \vphantom{1} 0.00\\
\hline
\hspace{1em}mlm & 100 & 99.6 & \vphantom{1} 0.00\\
\hline
\hspace{1em}mahdistmcd & 100 & 99.6 & \vphantom{1} 0.00\\
\hline
\hspace{1em}isoforest & 98 & 99.6 & 0.00\\
\hline
\hspace{1em}barbosa & 100 & 81.3 & \vphantom{1} 0.00\\
\hline
\hspace{1em}beta & 100 & 92.9 & 0.00\\
\hline
\multicolumn{4}{l}{\textbf{n100}}\\
\hline
\hspace{1em}manova & 100 & 99.6 & 0.00\\
\hline
\hspace{1em}mlm & 100 & 99.6 & 0.00\\
\hline
\hspace{1em}mahdistmcd & 100 & 99.6 & 0.00\\
\hline
\hspace{1em}isoforest & 99 & 99.6 & 0.00\\
\hline
\hspace{1em}barbosa & 100 & 81.3 & 0.00\\
\hline
\hspace{1em}beta & 100 & 90.3 & 0.00\\
\hline
\end{tabular}
\end{table}



\begin{figure}[H]

{\centering \includegraphics[width=27.99in]{C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021/result_files/4-Graph/GSE51032_chr17_46018653_46019185} 

}

\caption{epimutations performance for GSE51032 cohort detecting the epivariation located in chr5:10249760-10251253}\label{fig:graph_GSE51032}
\end{figure}

\begin{table}
\centering
\begin{tabular}[t]{c|c|c|c|c}
\hline
method & n & TPR & accuracy & FPR\\
\hline
\multicolumn{5}{l}{\textbf{n20}}\\
\hline
\hspace{1em}barbosa & 100 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}beta & 100 & 100.000 & 92.825 & 0.25\\
\hline
\hspace{1em}isoforest & 100 & 98.625 & 93.000 & 0.00\\
\hline
\hspace{1em}mahdistmcd & 100 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}manova & 100 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}mlm & 100 & 100.000 & 92.825 & 0.00\\
\hline
\multicolumn{5}{l}{\textbf{n30}}\\
\hline
\hspace{1em}barbosa & 20 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}beta & 20 & 87.500 & 92.800 & 0.00\\
\hline
\hspace{1em}isoforest & 20 & 11.500 & 86.800 & 0.00\\
\hline
\hspace{1em}mahdistmcd & 20 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}manova & 20 & 98.000 & 92.825 & 0.00\\
\hline
\hspace{1em}mlm & 20 & 100.000 & 92.850 & 0.00\\
\hline
\multicolumn{5}{l}{\textbf{n40}}\\
\hline
\hspace{1em}barbosa & 30 & 100.000 & 92.825 & 0.25\\
\hline
\hspace{1em}beta & 30 & 87.500 & 92.775 & 0.25\\
\hline
\hspace{1em}isoforest & 30 & 28.000 & 93.150 & 0.00\\
\hline
\hspace{1em}mahdistmcd & 30 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}manova & 30 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}mlm & 30 & 100.000 & 92.825 & 0.00\\
\hline
\multicolumn{5}{l}{\textbf{n50}}\\
\hline
\hspace{1em}barbosa & 40 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}beta & 40 & 87.500 & 92.800 & 0.00\\
\hline
\hspace{1em}isoforest & 40 & 46.375 & 93.950 & 0.00\\
\hline
\hspace{1em}mahdistmcd & 40 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}manova & 40 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}mlm & 40 & 100.000 & 92.825 & 0.00\\
\hline
\multicolumn{5}{l}{\textbf{n60}}\\
\hline
\hspace{1em}barbosa & 50 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}beta & 50 & 87.500 & 92.825 & 0.00\\
\hline
\hspace{1em}isoforest & 50 & 70.125 & 93.500 & 0.00\\
\hline
\hspace{1em}mahdistmcd & 50 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}manova & 50 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}mlm & 50 & 100.000 & 92.825 & 0.00\\
\hline
\multicolumn{5}{l}{\textbf{n70}}\\
\hline
\hspace{1em}barbosa & 60 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}beta & 60 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}isoforest & 60 & 78.750 & 93.500 & 0.00\\
\hline
\hspace{1em}mahdistmcd & 60 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}manova & 60 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}mlm & 60 & 100.000 & 92.825 & 0.00\\
\hline
\multicolumn{5}{l}{\textbf{n80}}\\
\hline
\hspace{1em}barbosa & 70 & 100.000 & 92.825 & 0.25\\
\hline
\hspace{1em}beta & 70 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}isoforest & 70 & 90.125 & 93.575 & 0.00\\
\hline
\hspace{1em}mahdistmcd & 70 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}manova & 70 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}mlm & 70 & 100.000 & 92.825 & 0.00\\
\hline
\multicolumn{5}{l}{\textbf{n90}}\\
\hline
\hspace{1em}barbosa & 80 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}beta & 80 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}isoforest & 80 & 96.375 & 93.075 & 0.00\\
\hline
\hspace{1em}mahdistmcd & 80 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}manova & 80 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}mlm & 80 & 100.000 & 92.825 & 0.00\\
\hline
\multicolumn{5}{l}{\textbf{n100}}\\
\hline
\hspace{1em}barbosa & 90 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}beta & 90 & 100.000 & 92.825 & 0.25\\
\hline
\hspace{1em}isoforest & 90 & 97.125 & 93.000 & 0.00\\
\hline
\hspace{1em}mahdistmcd & 90 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}manova & 90 & 100.000 & 92.825 & 0.00\\
\hline
\hspace{1em}mlm & 90 & 100.000 & 92.825 & 0.00\\
\hline
\end{tabular}
\end{table}

\begin{figure}[H]

{\centering \includegraphics[width=27.99in]{C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021/result_files/4-Graph/GSE111629_chr5_10249760_10251253} 

}

\caption{epimutations performance using GSE111629 cohort to detect the epivariation located in chr5:10249760-10251253}\label{fig:graph_GSE111629_1}
\end{figure}

\begin{figure}

{\centering \includegraphics[width=27.99in]{C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021/result_files/4-Graph/GSE111629_chr5_67583971_67584381} 

}

\caption{epimutations performance using GSE111629 cohort to detect the epivariation located in chr5:67583971-67584381}\label{fig:graph_GSE111629_2}
\end{figure}

\begin{figure}

{\centering \includegraphics[width=27.99in]{C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021/result_files/4-Graph/GSE111629_chr17_46018653_46019185} 

}

\caption{epimutations performance using GSE111629 cohort to detect the epivariation located in chr17:46018653-46019185}\label{fig:graph_GSE111629_3}
\end{figure}

\begin{figure}

{\centering \includegraphics[width=27.99in]{C:/Users/nla94/Documents/GitHub/Supplementary-Material/Abarrategui_2021/result_files/4-Graph/GSE111629_chr19_11199850_11200147} 

}

\caption{epimutations performance using GSE111629 cohort to detect the epivariation located in chr5:11199850-11200147}\label{fig:graph_GSE111629_4}
\end{figure}

\newpage

# Acknowledgements

We acknowledge the organizers of the [European BioHackathon 2020](https://www.biohackathon-europe.org/) for their support.

All the team members of *Project #5* for the contribution to this package: 

| Name | Surname | ORCID | Affiliation | Team |
|:-----|:--------|:-----:|:------------|:-----|
| Leire | Abarrategui | 0000-0002-1175-038X | Faculty of Medical Sciences, Newcastle University, Newcastle-Upon-Tyne, UK; Autonomous University of Barcelona (UAB), Barcelona, Spain | Development |
| Lordstrong | Akano | 0000-0002-1404-0295 | College of Medicine, University of Ibadan | Development |
| James | Baye | 0000-0002-0078-3688 | Wellcome/MRC Cambridge Stem Cell Institute, University of Cambridge, Cambridge CB2 0AW, UK; Department of Physics, University of Cambridge, Cambridge CB2 3DY, UK | Development |
| Alejandro | Caceres | - | ISGlobal, Barcelona Institute for Global Health, Dr Aiguader 88, 08003 Barcelona, Spain; Centro de Investigación Biomédica en Red en Epidemiología y Salud Pública (CIBERESP), Madrid, Spain | Development |
| Carles | Hernandez-Ferrer | 0000-0002-8029-7160 | Centro Nacional de Análisis Genómico (CNAG-CRG), Center for Genomic, Regulation; Barcelona Institute of Science and Technology (BIST), Barcelona, Catalonia, Spain | Development |		
| Pavlo | Hrab | 0000-0002-0742-8478 | Department of Genetics and Biotechnology, Biology faculty, Ivan Franko National University of Lviv | Validation |
| Raquel | Manzano | 0000-0002-5124-8992 | Cancer Research UK Cambridge Institute; University of Cambridge, Cambridge, United Kingdom | Reporting |
| Margherita | Mutarelli | 0000-0002-2168-5059 | Institute of Applied Sciences and Intelligent Systems (ISASI-CNR) | Validation |
| Carlos | Ruiz-Arenas | 0000-0002-6014-3498 | Centro de Investigación Biomédica en Red de Enfermedades Raras (CIBERER), Barcelona, Spain; Universitat Pompeu Fabra (UPF), Barcelona, Spain | Reporting |

# References

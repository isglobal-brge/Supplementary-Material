# Defective X chromosome inactivation and cancer risk in women 

Alejandro Cáceres, Luis A. Pérez-Jurado, Albert Alegret-García, Varun B. Dwaraka, Ryan Smith ans, Juan R. González

### Computer Code and Data

These resources comprise all the computer code and necessary data to replicate results in Caceres et al 2024. 

The computer code is written in R/Bioconductor and it is organized following the results section of the manuscript. 

The processed data folder is available at https://doi.org/10.6084/m9.figshare.25886674



**Section 1: X-Ra as a measure of quantitative XCI**  

- TCGA_XRA_CALLING
- TCGA_UNDISEASED_TISSUES

**Section 2: X-Ra as a marker of female cancer** 

- TCGA_XRA_CALLING
- TCGA_TUMORS

**Section 3: X-Ra as a marker of breast cancer stratification**

- TCGA_GEO_BRCA

**Section 4: X-Ra associations with somatic mutations in breast cancer**

- TCGA_BRCA_CNVs

**Section 5: X-Ra in blood is associated with aging and cancer** 

- TRUDIAGNOSTIC_BLOOD
- GENOA_BLOOD_AGE_META
- MESA_MONOCYTES
- GSE142536_BLOOD


Within the code files, figures and tables of the manuscript are indicated.

**Data content:**

Download data at https://figshare.com/account/articles/25886674

- Annotation data files
- Files with selection of X chromosome methylation CpGs for different studies
- Supplementary material used from other studies (cited within the code files)

chrXRa software can be obtained in: [https://xra.isglobal.org] with documentation and terms of licence.


Seesion info is as follows:

```
R version 4.3.2 (2023-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=Spanish_Spain.utf8  LC_CTYPE=Spanish_Spain.utf8    LC_MONETARY=Spanish_Spain.utf8
[4] LC_NUMERIC=C                   LC_TIME=Spanish_Spain.utf8    

time zone: Europe/Madrid
tzcode source: internal

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] methylclock_1.8.0           quadprog_1.5-8              devtools_2.4.5             
 [4] usethis_2.2.2               methylclockData_1.10.0      futile.logger_1.4.3        
 [7] irr_0.84.1                  lpSolve_5.6.19              visreg_2.7.0               
[10] arm_1.13-1                  lme4_1.1-35.1               Matrix_1.6-2               
[13] MASS_7.3-60                 minfi_1.48.0                bumphunter_1.44.0          
[16] locfit_1.5-9.8              iterators_1.0.14            foreach_1.5.2              
[19] Biostrings_2.70.1           XVector_0.42.0              meta_6.5-0                 
[22] GEOquery_2.70.0             XRa_1.0                     TCGAbiolinks_2.30.0        
[25] curatedTCGAData_1.24.0      MultiAssayExperiment_1.28.0 SummarizedExperiment_1.32.0
[28] MatrixGenerics_1.14.0       matrixStats_1.1.0           EnhancedVolcano_1.20.0     
[31] ggrepel_0.9.4               sva_3.50.0                  BiocParallel_1.36.0        
[34] genefilter_1.84.0           mgcv_1.9-0                  nlme_3.1-163               
[37] limma_3.58.1                AF_0.1.5                    ivtools_2.3.0              
[40] data.table_1.14.8           stdReg_3.4.1                drgee_1.1.10               
[43] pROC_1.18.5                 cvAUC_1.1.4                 survival_3.5-7             
[46] ggforestplot_0.1.0          GenomicRanges_1.54.1        GenomeInfoDb_1.38.1        
[49] ggplot2_3.4.4               org.Hs.eg.db_3.18.0         AnnotationDbi_1.64.1       
[52] IRanges_2.36.0              S4Vectors_0.40.1            biomaRt_2.58.0             
[55] clusterProfiler_4.10.0      chrXRa_1.0                  Biobase_2.62.0             
[58] BiocGenerics_0.48.1        

loaded via a namespace (and not attached):
  [1] R.methodsS3_1.8.2             progress_1.2.2                ggpp_0.5.5                   
  [4] urlchecker_1.0.1              RPMM_1.25                     HDF5Array_1.30.0             
  [7] vctrs_0.6.4                   digest_0.6.33                 png_0.1-8                    
 [10] BiocBaseUtils_1.4.0           PerformanceAnalytics_2.0.4    reshape_0.8.9                
 [13] reshape2_1.4.4                httpuv_1.6.12                 qvalue_2.34.0                
 [16] withr_2.5.2                   xfun_0.41                     ggfun_0.1.3                  
 [19] ggpubr_0.6.0                  ellipsis_0.3.2                doRNG_1.8.6                  
 [22] memoise_2.0.1                 MatrixModels_0.5-3            gson_0.1.0                   
 [25] profvis_0.3.8                 systemfonts_1.0.5             zoo_1.8-12                   
 [28] ragg_1.2.6                    tidytree_0.4.5                R.oo_1.25.0                  
 [31] prettyunits_1.2.0             KEGGREST_1.42.0               promises_1.2.1               
 [34] httr_1.4.7                    downloader_0.4                rstatix_0.7.2                
 [37] restfulr_0.0.15               rhdf5filters_1.14.1           ps_1.7.5                     
 [40] rhdf5_2.46.0                  rstudioapi_0.15.0             miniUI_0.1.1.1               
 [43] generics_0.1.3                DOSE_3.28.1                   processx_3.8.2               
 [46] curl_5.1.0                    zlibbioc_1.48.0               ggraph_2.1.0                 
 [49] polyclip_1.10-6               GenomeInfoDbData_1.2.11       ExperimentHub_2.10.0         
 [52] SparseArray_1.2.2             RBGL_1.78.0                   interactiveDisplayBase_1.40.0
 [55] xtable_1.8-4                  stringr_1.5.0                 evaluate_0.23                
 [58] S4Arrays_1.2.0                BiocFileCache_2.10.1          preprocessCore_1.64.0        
 [61] hms_1.1.3                     colorspace_2.1-0              filelock_1.0.2               
 [64] polynom_1.4-1                 ROCR_1.0-11                   magrittr_2.0.3               
 [67] readr_2.1.4                   later_1.3.1                   viridis_0.6.4                
 [70] ggtree_3.10.0                 lattice_0.22-5                SparseM_1.81                 
 [73] XML_3.99-0.15                 shadowtext_0.1.2              cowplot_1.1.1                
 [76] xts_0.13.1                    pillar_1.9.0                  compiler_4.3.2               
 [79] stringi_1.7.12                minqa_1.2.6                   GenomicAlignments_1.38.0     
 [82] plyr_1.8.9                    crayon_1.5.2                  abind_1.4-5                  
 [85] BiocIO_1.12.0                 metadat_1.2-0                 gridGraphics_0.5-1           
 [88] graphlayouts_1.0.2            bit_4.0.5                     mathjaxr_1.6-0               
 [91] dplyr_1.1.3                   fastmatch_1.1-4               codetools_0.2-19             
 [94] textshaping_0.3.7             openssl_2.1.1                 multtest_2.58.0              
 [97] mime_0.12                     splines_4.3.2                 Rcpp_1.0.11                  
[100] quantreg_5.97                 BiocCheck_1.38.0              dbplyr_2.4.0                 
[103] sparseMatrixStats_1.14.0      TCGAbiolinksGUI.data_1.22.0   HDO.db_0.99.1                
[106] knitr_1.45                    blob_1.2.4                    utf8_1.2.4                   
[109] BiocVersion_3.18.0            fs_1.6.3                      DelayedMatrixStats_1.24.0    
[112] pkgbuild_1.4.2                ggsignif_0.6.4                ggplotify_0.1.2              
[115] tibble_3.2.1                  callr_3.7.3                   statmod_1.5.0                
[118] tzdb_0.4.0                    tweenr_2.0.2                  pkgconfig_2.0.3              
[121] tools_4.3.2                   cachem_1.0.8                  RSQLite_2.3.3                
[124] viridisLite_0.4.2             rvest_1.0.3                   DBI_1.1.3                    
[127] numDeriv_2016.8-1.1           impute_1.76.0                 rmarkdown_2.25               
[130] fastmap_1.1.1                 scales_1.2.1                  grid_4.3.2                   
[133] sdamr_0.2.0                   Rsamtools_2.18.0              broom_1.0.5                  
[136] AnnotationHub_3.10.0          patchwork_1.1.3               coda_0.19-4                  
[139] BiocManager_1.30.22           carData_3.0-5                 graph_1.80.0                 
[142] scrime_1.3.5                  farver_2.1.1                  tidygraph_1.2.3              
[145] scatterpie_0.2.1              AnnotationForge_1.44.0        biocViews_1.70.0             
[148] yaml_2.3.7                    rtracklayer_1.62.0            illuminaio_0.44.0            
[151] cli_3.6.1                     purrr_1.0.2                   siggenes_1.76.0              
[154] ahaz_1.15                     lifecycle_1.0.4               askpass_1.2.0                
[157] ggpmisc_0.5.5                 sessioninfo_1.2.2             lambda.r_1.2.4               
[160] backports_1.4.1               annotate_1.80.0               gtable_0.3.4                 
[163] rjson_0.2.21                  metafor_4.4-0                 ape_5.7-1                    
[166] jsonlite_1.8.7                edgeR_4.0.1                   bitops_1.0-7                 
[169] bit64_4.0.5                   yulab.utils_0.1.0             base64_2.0.1                 
[172] futile.options_1.0.1          GOSemSim_2.28.0               R.utils_2.12.2               
[175] lazyeval_0.2.2                dynamicTreeCut_1.63-1         shiny_1.7.5.1                
[178] htmltools_0.5.7               enrichplot_1.22.0             GO.db_3.18.0                 
[181] formatR_1.14                  rappdirs_0.3.3                glue_1.6.2                   
[184] RCurl_1.98-1.13               treeio_1.26.0                 mclust_6.0.1                 
[187] gridExtra_2.3                 boot_1.3-28.1                 igraph_1.5.1                 
[190] R6_2.5.1                      tidyr_1.3.0                   CompQuadForm_1.4.3           
[193] labeling_0.4.3                GenomicFeatures_1.54.1        cluster_2.1.4                
[196] rngtools_1.5.2                pkgload_1.3.3                 Rhdf5lib_1.24.0              
[199] AnnotationHubData_1.32.0      stringdist_0.9.12             beanplot_1.3.1               
[202] aplot_0.2.2                   nloptr_2.0.3                  DelayedArray_0.28.0          
[205] tidyselect_1.2.0              nleqslv_3.3.5                 ggforce_0.4.1                
[208] xml2_1.3.5                    car_3.1-2                     munsell_0.5.0                
[211] nor1mix_1.3-2                 htmlwidgets_1.6.2             fgsea_1.28.0                 
[214] RColorBrewer_1.1-3            rlang_1.1.2                   tidyverse_2.0.0              
[217] remotes_2.4.2.1               lmerTest_3.1-3                RUnit_0.4.32                 
[220] fansi_1.0.5                   ExperimentHubData_1.28.0      OrganismDbi_1.44.0

```



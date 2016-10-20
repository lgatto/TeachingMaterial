# Wrap-up



## References and resources

* [Visualisation of proteomics data using R and Bioconductor](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4510819/)
* [Using R and Bioconductor for proteomics data analysis](http://arxiv.org/pdf/1305.6559v1.pdf)
* `RforProteomics`: http://bioconductor.org/packages/RforProteomics
* [R/Bioconductor work-flow](http://bioconductor.org/help/workflows/proteomics/)

## Other relevant packages/pipelines

- Analysis of post translational modification with *[isobar](http://bioconductor.org/packages/isobar)*.
- Processing and analysis or isobaric tagging mass spectrometry with
  *[isobar](http://bioconductor.org/packages/isobar)* and *[MSnbase](http://bioconductor.org/packages/MSnbase)*.
- Analysis of spatial proteomics data with *[pRoloc](http://bioconductor.org/packages/pRoloc)*.
- Analysis of MALDI data with the *[MALDIquant](http://bioconductor.org/packages/MALDIquant)* package.
- Access to the Proteomics Standard Initiative Common QUery InterfaCe
  with the *[PSICQUIC](http://bioconductor.org/packages/PSICQUIC)* package.
- *[Cardinal](http://bioconductor.org/packages/Cardinal)*: A mass spectrometry imaging toolbox for
  statistical analysis.
- *[protViz](http://cran.fhcrc.org/web/packages/protViz/index.html)*: Visualising and Analysing Mass Spectrometry
  Related Data in Proteomics
- *[aLFQ](http://cran.fhcrc.org/web/packages/aLFQ/index.html)*: Estimating Absolute Protein Quantities from
  Label-Free LC-MS/MS Proteomics Data.
- *[protiq](http://cran.fhcrc.org/web/packages/protiq/index.html)*: Protein (identification and) quantification
  based on peptide evidence.
- *[MSstats](http://bioconductor.org/packages/MSstats)*: Protein Significance Analysis in DDA, SRM
  and DIA for Label-free or Label-based Proteomics Experiments


### DIA

- Analysis of label-free data from a Synapt G2 (including ion
  mobility) with *[synapter](http://bioconductor.org/packages/synapter)*.
- *[SWATH2stats](http://bioconductor.org/packages/SWATH2stats)*: Transform and Filter SWATH Data for
  Statistical Packages and
- *[specL](http://bioconductor.org/packages/specL)*: Prepare Peptide Spectrum Matches for Use in
  Targeted Proteomics
- *[SwathXtend](http://bioconductor.org/packages/SwathXtend)*: SWATH extended library generation and
  statistical data analysis


### Session info


```r
sessionInfo()
```

```
## R version 3.3.1 Patched (2016-08-02 r71022)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.5 LTS
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    parallel  methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] BiocStyle_2.1.33      qvalue_2.5.2          MSnID_1.7.3          
##  [4] msmsTests_1.11.0      msmsEDA_1.11.0        limma_3.29.21        
##  [7] multtest_2.29.0       RColorBrewer_1.1-2    ggplot2_2.1.0        
## [10] magrittr_1.5          hexbin_1.27.1         dplyr_0.5.0          
## [13] readxl_0.1.1          gridExtra_2.2.1       RforProteomics_1.11.2
## [16] mzID_1.11.2           msdata_0.12.4         lattice_0.20-34      
## [19] pRolocdata_1.11.9     pRoloc_1.13.17        MLInterfaces_1.53.1  
## [22] cluster_2.0.5         annotate_1.51.1       XML_3.98-1.4         
## [25] AnnotationDbi_1.35.4  IRanges_2.7.17        S4Vectors_0.11.19    
## [28] MSnbase_1.99.7        ProtGenerics_1.5.1    BiocParallel_1.7.9   
## [31] mzR_2.7.12            Rcpp_0.12.7           Biobase_2.33.4       
## [34] BiocGenerics_0.19.2   gplots_3.0.1          knitr_1.14           
## 
## loaded via a namespace (and not attached):
##   [1] plyr_1.8.4                    GSEABase_1.35.5              
##   [3] splines_3.3.1                 ggvis_0.4.3                  
##   [5] digest_0.6.10                 foreach_1.4.3                
##   [7] BiocInstaller_1.23.9          htmltools_0.3.5              
##   [9] gdata_2.17.0                  doParallel_1.0.10            
##  [11] sfsmisc_1.1-0                 rda_1.0.2-2                  
##  [13] R.utils_2.4.0                 lpSolve_5.6.13               
##  [15] colorspace_1.2-7              RCurl_1.95-4.8               
##  [17] jsonlite_1.1                  graph_1.51.0                 
##  [19] genefilter_1.55.2             lme4_1.1-12                  
##  [21] impute_1.47.0                 survival_2.39-5              
##  [23] iterators_1.0.8               gtable_0.2.0                 
##  [25] zlibbioc_1.19.0               MatrixModels_0.4-1           
##  [27] R.cache_0.12.0                car_2.1-3                    
##  [29] kernlab_0.9-25                prabclus_2.2-6               
##  [31] DEoptimR_1.0-6                SparseM_1.72                 
##  [33] scales_0.4.0                  vsn_3.41.5                   
##  [35] mvtnorm_1.0-5                 edgeR_3.15.6                 
##  [37] DBI_0.5-1                     xtable_1.8-2                 
##  [39] proxy_0.4-16                  mclust_5.2                   
##  [41] preprocessCore_1.35.0         htmlwidgets_0.7              
##  [43] sampling_2.7                  threejs_0.2.2                
##  [45] FNN_1.1                       fpc_2.1-10                   
##  [47] modeltools_0.2-21             R.methodsS3_1.7.1            
##  [49] flexmix_2.3-13                nnet_7.3-12                  
##  [51] locfit_1.5-9.1                RJSONIO_1.3-0                
##  [53] caret_6.0-71                  reshape2_1.4.1               
##  [55] munsell_0.4.3                 mlbench_2.1-1                
##  [57] biocViews_1.41.9              tools_3.3.1                  
##  [59] RSQLite_1.0.0                 pls_2.5-0                    
##  [61] evaluate_0.10                 stringr_1.1.0                
##  [63] robustbase_0.92-6             caTools_1.17.1               
##  [65] randomForest_4.6-12           dendextend_1.3.0             
##  [67] RBGL_1.49.3                   nlme_3.1-128                 
##  [69] whisker_0.3-2                 mime_0.5                     
##  [71] quantreg_5.29                 formatR_1.4                  
##  [73] R.oo_1.20.0                   biomaRt_2.29.2               
##  [75] pbkrtest_0.4-6                interactiveDisplayBase_1.11.3
##  [77] e1071_1.6-7                   affyio_1.43.0                
##  [79] tibble_1.2                    stringi_1.1.2                
##  [81] rpx_1.9.4                     trimcluster_0.1-2            
##  [83] Matrix_1.2-7.1                nloptr_1.0.4                 
##  [85] gbm_2.1.1                     RUnit_0.4.31                 
##  [87] MALDIquant_1.15               data.table_1.9.6             
##  [89] bitops_1.0-6                  httpuv_1.3.3                 
##  [91] R6_2.2.0                      pcaMethods_1.65.0            
##  [93] affy_1.51.1                   hwriter_1.3.2                
##  [95] KernSmooth_2.23-15            gridSVG_1.5-0                
##  [97] codetools_0.2-15              MASS_7.3-45                  
##  [99] gtools_3.5.0                  assertthat_0.1               
##  [ reached getOption("max.print") -- omitted 11 entries ]
```

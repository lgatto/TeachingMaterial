# Wrap-up



## References and resources

* [Visualisation of proteomics data using R and Bioconductor](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4510819/)
* [Using R and Bioconductor for proteomics data analysis](http://arxiv.org/pdf/1305.6559v1.pdf)
* `RforProteomics`: http://bioconductor.org/packages/RforProteomics
* [R/Bioconductor workflow](http://bioconductor.org/help/workflows/proteomics/)

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
## [1] parallel  methods   stats     graphics  grDevices utils     datasets 
## [8] base     
## 
## other attached packages:
##  [1] BiocStyle_2.1.33      qvalue_2.5.2          MSnID_1.7.3          
##  [4] msmsTests_1.11.0      msmsEDA_1.11.0        limma_3.29.21        
##  [7] multtest_2.29.0       RColorBrewer_1.1-2    ggplot2_2.1.0        
## [10] magrittr_1.5          hexbin_1.27.1         readxl_0.1.1         
## [13] dplyr_0.5.0           gridExtra_2.2.1       RforProteomics_1.11.2
## [16] mzID_1.11.2           MSnbase_1.99.7        ProtGenerics_1.5.1   
## [19] BiocParallel_1.7.9    mzR_2.7.12            Rcpp_0.12.7          
## [22] Biobase_2.33.4        BiocGenerics_0.19.2   msdata_0.12.4        
## [25] lattice_0.20-34       knitr_1.14           
## 
## loaded via a namespace (and not attached):
##  [1] Category_2.39.0               bitops_1.0-6                 
##  [3] doParallel_1.0.10             R.cache_0.12.0               
##  [5] tools_3.3.1                   R6_2.2.0                     
##  [7] affyio_1.43.0                 KernSmooth_2.23-15           
##  [9] DBI_0.5-1                     colorspace_1.2-7             
## [11] preprocessCore_1.35.0         chron_2.3-47                 
## [13] graph_1.51.0                  formatR_1.4                  
## [15] caTools_1.17.1                scales_0.4.0                 
## [17] genefilter_1.55.2             affy_1.51.1                  
## [19] RBGL_1.49.3                   stringr_1.1.0                
## [21] digest_0.6.10                 R.utils_2.4.0                
## [23] htmltools_0.3.5               RSQLite_1.0.0                
## [25] impute_1.47.0                 BiocInstaller_1.23.9         
## [27] shiny_0.14.1                  gtools_3.5.0                 
## [29] R.oo_1.20.0                   RCurl_1.95-4.8               
## [31] MALDIquant_1.15               Matrix_1.2-7.1               
## [33] munsell_0.4.3                 S4Vectors_0.11.19            
## [35] R.methodsS3_1.7.1             vsn_3.41.5                   
## [37] stringi_1.1.2                 edgeR_3.15.6                 
## [39] MASS_7.3-45                   RJSONIO_1.3-0                
## [41] zlibbioc_1.19.0               gplots_3.0.1                 
## [43] biocViews_1.41.9              plyr_1.8.4                   
## [45] rpx_1.9.4                     grid_3.3.1                   
## [47] gdata_2.17.0                  splines_3.3.1                
## [49] annotate_1.51.1               locfit_1.5-9.1               
## [51] interactiveDisplay_1.11.2     RUnit_0.4.31                 
## [53] reshape2_1.4.1                codetools_0.2-15             
## [55] stats4_3.3.1                  XML_3.98-1.4                 
## [57] evaluate_0.10                 pcaMethods_1.65.0            
## [59] data.table_1.9.6              httpuv_1.3.3                 
## [61] foreach_1.4.3                 gtable_0.2.0                 
## [63] assertthat_0.1                mime_0.5                     
## [65] xtable_1.8-2                  survival_2.39-5              
## [67] tibble_1.2                    iterators_1.0.8              
## [69] AnnotationDbi_1.35.4          IRanges_2.7.17               
## [71] interactiveDisplayBase_1.11.3 gridSVG_1.5-0                
## [73] GSEABase_1.35.5
```

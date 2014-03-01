**This is your third and final homework due March 13.**
  

You will have to analyze the RNA-seq data presented in:
Henn, A. D. et al. High-resolution temporal response patterns to influenza vaccine reveal a distinct human plasma cell gene signature. Scientific Reports 3, 2327 (2013).
  
1. Get the data from GEO. Please look at the class lecture slides as we've done it already
2. Use voom and limma to find genes that are differentially expressed at each time point compared to baseline (day 0). Use an FDR cutoff of 0.01.
Display your results using pheatmap showing the log fold-change of the differentially expressed genes grouped by time point.  
3. Perform a GSEA analysis using camera and the MSigDB Reactome pathway gene signatures. Display your results using pheatmap, again group by timepoint. This is similar to what we've done in class.
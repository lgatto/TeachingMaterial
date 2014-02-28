# Biostat578 Midterm Exam

-------

## Name:

## Student Number:

---------
This is an **open** book/computer exam. Feel free to use whatever resources you need, but work on your own. There are 12 questions for a total of  100 points. You have until 11:30am to complete this exam. Good luck!
---------

- 1) [5 points] `X` and `Y` are two datatables in `R`, how would do an inner join of `X` and `Y`.

**Solution:** Your first need to use `setkey(X,"key")`, `setkey(Y,"key")`. Then you can do `X[Y,match=0]`. This is faster than `merge`.

- 2) [10 points] Let `Z` be a datatable in `R` with three columns: `student_name`, `grade`, and `assignment_type`. Grades are numerical values from 0 to 100. `assignment_type` can be "assignment" or "in-class exam". Write a one line command to compute the average grade for each student broken down by `assignment_type`, and sort by grade from lowest to highest.

**Solution:** Z[,list(mgrade=mean(grade)), by="assignment_type"][order(mgrade)]

- 3) [10 points] Let `E` be an `ExpressionSet` in `R`. How would you retrieve the expression matrix, the probe information and the sample information? [State the name of the methods (i.e. R function) you would use].

exprs(E), fData(E), pData(E)

- 4) [10 points] State three advantages of next generation sequencing over microarrays. Explain your answer. 

**Solution:** No cross-hybridization, more accurate detection of low expression transcripts, can be used to study many biological processes.

- 5) [5 points] What is the main idea behind quantile normalization?

**Solution:** To make the distribution of two (or more) microarrays the same, which is valid when we assume most gene are not changing. 

- 6) [5 points] What is the main idea behind lowess normalization?

**Solution:** To remove the non-linear dependence between the M and A values. 

- 7) [10 points] You've performed a gene expression experiment using microarrays. Unfortunately, some of your samples were processed in June, while others were processed in October. You suspect a batch effect and would like to correct for it in your limma analysis. How would you do that? **You may assume that your contrast of interest is not confounded with batch order**

**Solution:** Simply include the batch variable (time information) in your limma analysis.

- 8) [5 points] What is the false discovery rate? Why is it preferable to control this versus the family-wise error rate?

**Solution:** The false discovery rate (FDR) is the proportion of false positive among your discoveries. In large scale genomics, it is preferable to the FWER because it is less conservative. 

- 9) [15 points] You have a gene expression experiment looking at expression changes in a cohort of 40 subjects before and after drug treatment. You follow the subjects over time and have samples at day 0 (day of drug treatment), day 14 and day 28. How would do set up your design matrix in limma (write the actual R command), what contrasts would you test?

**Solution:** 
`mm <- model.matrix(~day+subject, eSet)`

Then test for the contrasts: "day28-day0", "day14-day0"

- 10) [10 points] What is the main idea behind `limma`? Why is it preferable over a traditional linear model?


**Solution:** Limma uses an empirical Bayes approach to regularize variances that are then used to form a modified t-statistics. This approach leads to better performance when the number of replicates is small, which is often the case in genomics.

- 11) [5 points] `limma` was derived for gene microarrays. It's been shown recently that it could also be applied to next generation sequencing data after proper data transformation using `voom`. What is `voom` actually doing?

**Solution:** `voom` is used to prepare RNA-seq count data for `limma` analysis. The `voom` method estimates the mean-variance relationship of the log-counts, generates a precision weight for each observation, and then enters these into a limma empirical Bayes analysis pipeline.

- 12) [10 points] Do you need to normalize next generation sequencing data? If yes, what procedure would you use?

**Solution:** Yes, you need to normalize RNA-seq data due to variation in sequencing depth. Without normalization the results would not be comparable across libraries.

- 13) [Bonus, 10 points]. Prove that the Bonferroni mutiple testing procedure controls the familywise error rate at the desired level. 

**Solution:** This comes directly from Boole's inequality. 
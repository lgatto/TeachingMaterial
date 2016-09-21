# R OO Programming and Package Development

This is the material for the 1-day course taught at the University of
Cambridge.

For the course IT preparation, we need to be able to build and check
the [`sequences`](https://github.com/lgatto/sequences). package. The
best way to test if this is possible is to run the following code
chunk.

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("lgatto/sequences", quick = FALSE, dependencies = TRUE)
```

For more material, see the
[TeachingMaterial](http://lgatto.github.io/TeachingMaterial/)
site/repository.

All the material is available under a
[Creative Commons Attribution-ShareAlike 3.0 License](http://creativecommons.org/licenses/by-sa/3.0/).

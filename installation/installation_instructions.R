packages <- c("ggplot2", "reshape2", "psych", "pwr",
              "MSstats", "limma", "MSnbase", "isobar",
	      "aLFQ", "MSnID", "RforProteomics")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(packages)
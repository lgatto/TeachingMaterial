# ----------------------------------------------------------------------
# installallpkgs.R
# ----------------------------------------------------------------------
#
# This script attempts to find all of the R packages used in all
# of the R presentation files in the current working folder and then 
# attempts to install them, either with "install.packages" or "biocLite",
# then tries to load them into memory with "require". The intent is to
# ensure that any packages required for the scripts to run have been
# installed previously. 
#
# Only those packages which are not already installed will be installed. 
# The "require" function is used to check those packages after 
# they are installed to make sure they can load into memory. Therefore, 
# even if packages have already been loaded with "library", this script 
# should run very quickly once the packages have been initially installed.
#
# This script was initially developed for use with course materials of
# "BIOSTAT 578A: Bioinformatics for Big Omics Data", specifically the  
# files (*.Rpres) found at https://github.com/raphg/Biostat-578/
#
# The idea is to download the course materials with "git clone", then
# run this script to install all of the packages used in the Rpres files
# so that the slides will display properly when run using, e.g., RStudio.
# Since many of the presentaions do not include installation commands 
# for all of the required packages, the slides may not display properly
# without installing packages beforehand.
#
# Usage: 
#
#   1. Use setwd() to set the working directory to the folder containing
#      the R script and R presentation files.
#
#   2. Run this command: source("installallpkgs.R")
#      ... where you should use the actual file path to this script.
#
# Author: Brian High
# Date: 2014-11-29
# License: http://creativecommons.org/licenses/by-sa/3.0/deed.en_US

# ----------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------

# tryinstall() function
#     Conditionally install using install.packages
#            or biocLite()
#     Usage: tryinstall(c("package1", "package2", ...))
tryinstall <- function(p) {
    n <- p[!(p %in% installed.packages()[,"Package"])]
    if(length(n)) {
        install.packages(n, repos="http://cran.fhcrc.org") | {
            source("http://bioconductor.org/biocLite.R")
            biocLite(n, ask = FALSE)
        }
    }
}

# ----------------------------------------------------------------------
# Main Routine
# ----------------------------------------------------------------------

# Compile a list of R and Rpres filenames in the current directory
# Tip: Use pattern="*.Rpres|.R" to also check *.R scripts.
filenames <- list.files(".", pattern="*.Rpres", full.names=FALSE)

# Parse each file to find the packages used and compile into a list
allpkgs <- c()
for (filename in filenames) {
    pkgs <- readLines(filename, warn = FALSE)
    pkgs <- unlist(strsplit(x = pkgs, split = ";[ ]*"))
    pkgs <- pkgs[grepl("(library|require|install\\.packages|biocLite)\\(", pkgs)]
    pkgs <- gsub(".*\\((.*)\\).*", "\\1", pkgs)
    pkgs <- unlist(strsplit(x = pkgs, split = ",[ ]*"))
    pkgs <- gsub('["()]', "", pkgs)
    pkgs <- unique(pkgs[!grepl("=", pkgs)])
    allpkgs <- c(allpkgs,pkgs)
}

# Remove duplicates
allpkgs <- unique(allpkgs)

# Save a copy of the package list
write(allpkgs, "packages_list.txt")

# Attempt to install each package using the tryinstall function
for (pkg in allpkgs) tryinstall(pkg)

# Uncomment the next line to install the development version of data.table.
# install.packages("data.table", repos="http://R-forge.R-project.org")

# Check that all packages will load
for (pkg in allpkgs) require(pkg, character.only = TRUE)

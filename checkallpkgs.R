# ----------------------------------------------------------------------
# checkallpkgs.R
# ----------------------------------------------------------------------
#
# This script attempts to find all of the R packages used in all
# of the R presentation files in the current working folder and then 
# attempts to install them with "biocLite" then tries to load them 
# into memory with "require". The intent is to ensure that any packages 
# required for the scripts to run have been installed previously. 
#
# Only those packages which are not already installed will be installed. 
# The "require" function is used to check those packages after 
# they are installed to make sure they can load into memory. Therefore, 
# even if packages have already been loaded with "library", this script 
# should run very quickly once the packages have been initially installed.
#
# This script was initially developed for use with course materials of
# "BIOSTAT 578A: Bioinformatics for Big Omics Data", specifically the  
# files (*.Rpres, *.Rmd) found at https://github.com/raphg/Biostat-578/
#
# The idea is to download the course materials with "git clone", then
# run this script to install all of the packages used in presentation files
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
#   2. Run this command: source("checkallpkgs.R")
#      ... where you should use the actual file path to this script.
#
# Author: Brian High with contributions from Rafael Gottardo
# Date: 2015-10-06
# License: http://creativecommons.org/licenses/by-sa/3.0/deed.en_US

# ----------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------

# tryinstall() function
#     Conditionally install using biocLite()
#     Usage: tryinstall(c("package1", "package2", ...))
tryinstall <- function(p) {
    n <- p[!(p %in% installed.packages()[,"Package"])]
    if(length(n)) {
        biocLite(n, ask = FALSE)
    }
}

# ----------------------------------------------------------------------
# Main Routine
# ----------------------------------------------------------------------

# Source biocLite so that it we can use biocLite() to install packages.
source("http://bioconductor.org/biocLite.R")

# Compile a list of Rpres and Rmd filenames in the current directory
# Tip: Use pattern="*.(Rpres|R|Rmd)" to also check *.R scripts.
filenames <- list.files(".", pattern="*.(Rpres|Rmd)", full.names=FALSE)

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

# Remove non-packages
allpkgs <- allpkgs[allpkgs!="txdbpkg"]

# Save a copy of the package list
write(allpkgs, "packages_list.txt")

# Attempt to load or install the packages using the tryinstall function
tryinstall(allpkgs)

# Check that all packages will load
for (pkg in allpkgs) require(pkg, character.only = TRUE)

# Review any warnings
warnings()

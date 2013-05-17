#################################################################
## 03_toRscript.R                                              ##
## Extract R script from Sweave files by package name argument ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###

# Load functions from tools library.  Used for parsing R code from Sweave docs
library("tools")

# Define user function for R script extraction
# Function has one argument: pkgn = package name
toRscript<-function(pkgn){

# vigSrc is a vector of vignette sweave file paths, obtained using list files.  File filter is *.Rnw, directory to search
# is doc of given package [pkgn].  List files will return full paths.  System.file function ensures paths are created appropriate to OS  
	vigSrc <- list.files(pattern = "Rnw$",
		system.file("doc", package = pkgn), full.names = TRUE)

# Output to console file paths
	vigSrc 

# Loop through items in vigSrc.
		for (v in vigSrc) {

# Strangle is the tools function responsible for extracting R code
			Stangle(v)

		}
# Function doesn't return a value.
}

### END ###

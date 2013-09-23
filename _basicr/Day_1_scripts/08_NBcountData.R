#################################################################
## 08_NBcountData.R                                            ##
## Exercise in reading in data, returning basic information    ##
## and exporting objects to spreesheets as CSV files           ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###
## Read in tab delimitted text file of results
## Remove index column
rawData <- read.delim("08_NBcountData.txt")
rawData[1:10,]
rawData <- rawData[,-1] # Drop the patient index column

## Getting basics summary data
nrow(rawData) 			# number of rows
ncol(rawData)  			# number of columns
dim(rawData)   			# data frame dimensions (rows X columns)‏
colMeans(rawData) 		# average stats

## Reorder table by decreasing nuclei count
rawDataOrder <- rawData[order(rawData[,1],decreasing=T),]
rawDataOrder

#  Identify patients with >33% NB amplification
prop <- rawData$NB_Amp / rawData$Nuclei
amp <- which(prop > 0.33)

plot(prop, ylim=c(0,1.2)) # plot a simple chart of NB amplifications
abline(h=0.33, lwd=1.5, lty=2) # Add a dotted line at 33%

# Write out results table as comma separated values file
write.csv(rawData[amp,],file="selectedSamples.csv")

#  In case you've forgotten where they've been saved
getwd() # prints working directory path, so you know where to look for the results file

### Solution to problem
### Which samples are near normal i.e. NBamp% < 0.33 & NBdel==0)
norm <- which((rawData$NB_Amp / rawData$Nuclei) < 0.33 & rawData$NB_Del==0)
write.csv(rawData[norm,],file="My_NB_output.csv")
### END ###


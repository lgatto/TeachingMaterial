#################################################################
## 10_heightExercise.R                                         ##
## Exercise in perfomring classical tests                      ##
## Simple t test                                               ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ##

# Read in data Could also use read.delim
heightData<-read.csv("10_heightData.csv",header=T)


# Display variance of Male and Female columns on console
var(heightData$Male) ; var(heightData$Female)

# Assign Fisher's to var.height
var.height <- var(heightData$Female) / var(heightData$Male)

# Display it on console
var.height

# Compute quantiles of F distribution (critical F) for above
qf(0.975,99,99)

# Determine 2 sided P Value associated with above F stat
2*pf(var.height, 99, 99, lower.tail=FALSE)

# What is the critical T for this
qt(0.975,198)

# Do an inbuilt t test
t.test(heightData$Male, heightData$Female)

# Do an inbuilt t test
wilcox.test(heightData$Male, heightData$Female)

### END ###


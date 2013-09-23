#################################################################
## 11_CFAcountData.R                                           ##
## Exercise in writing scripts with multiple steps             ##
## Script assigns several factor objects.  Reads in 3 data     ##
## files.  Undertakes some analysis and generates results      ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###
###############################

###########################
## Part 1.  Read in data ##
###########################

# load in the data from the three into three separate data frames
t1 = read.csv("11_CFA_Run1Counts.csv")
t2 = read.csv("11_CFA_Run2Counts.csv")
t3 = read.csv("11_CFA_Run3Counts.csv")

# concatenate the three data frames
colony = rbind(t1, t2, t3)

# verify we loaded factors
colony$Plate

# add the missing Run column - factors are stored as numbers ! 
runNum <- factor(rep(1:3,each=36),labels=c("Run1","Run2","Run3")) 
colony <- cbind( "Run" = runNum, colony )

# reorder factor levels in their natural order (instead of alphabetical)
colony$Treatment <- factor(colony$Treatment, c("OZ", "LZ", "MZ", "HZ"))
colony$Plate <- factor(colony$Plate, c("KDX","SCX","CNX"))

# show the full table
colony

###################################
### Part 2.  Investigating data ###
###################################
tapply(colony$Count, list(colony$Run, colony$Plate, colony$Treatment), mean)

#jpeg(file="fig1.jpg",width=675,height=900,res=150)
#par(oma=c(4,2,2,2))
boxplot(Count~Run*Plate*Treatment, las=2, cex=0.2, data=colony)
#dev.off()

#jpeg(file="fig2.jpg",width=675,height=900,res=150)
barplot(tapply(colony$Count, list(colony$Plate, colony$Treatment), mean), beside=T)
#dev.off()

#################################
### Part 3.  Summarizing data ###
#################################
result <- tapply(colony$Count, list(colony$Treatment, colony$Plate), mean)
result <- data.frame(result) # result of tapply is a matrix, convert back to data frame
result

# calculate kdx and scx values after background correction
kdx <- result$KDX - result$CNX
scx <- result$SCX - result$CNX

result <- cbind(kdx, scx)
# remove the 0Z entry
result <- result[-1,]

result
#jpeg(file="fig3.jpg",width=675,height=900,res=150)
barplot(result,beside=T)
#dev.off()

wilcox.test(result[,1],result[,2],paired=T)
cor.test(result[,1],result[,2],paired=T)
write.csv(result,"CFAresults.csv")

### END ###

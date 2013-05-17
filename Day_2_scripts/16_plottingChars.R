#################################################################
## 15_plottingChars.R                                          ##
## Graphics example, generates a jpg image of the 25 std       ##
## plotting characters in R                                    ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###

# Initiate file to save the plot in
jpeg("plotChar.jpg",width=1600,height=1600,res=300)

# x-y coordinates, plotting character index counter
xCounter<-1
yCounter<-1

# Plotting character index counter
plotChar<-0

# Create empty plot with the correct dimensions etc.
plot(NULL, xlim=c(0,8), ylim=c(0,5),xaxt="n",yaxt="n",ylab="",xlab="",main="26 standard plotting characters")

# Fill out the plot
while (plotChar <26){

	if(xCounter<7){
		xCounter<-xCounter+1
	} else {
		xCounter<-1
		yCounter<-yCounter+1
	}

	points(xCounter,yCounter,pch=plotChar, cex=2)	# Makes the point
	text(xCounter,(yCounter-0.3),plotChar)		# Adds the number just below the point
	plotChar<-plotChar+1				# Prepares for next iteration
}

dev.off()


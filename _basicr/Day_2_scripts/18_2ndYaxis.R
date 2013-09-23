#################################################################
## 17_2ndYaxis.R                                               ##
## Graphics example, how to add a second Y axis using          ##
## par(new=T)                                                  ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###

# Generate some fake data
x1<-1:20
y1<-sample(1000,20)
y2<-runif(20)

# labels for y2 axis
y2axis<-seq(0,1,.2)

# Set up figure margins
par(mar=c(4,4,4,4))

# Draw the first plot x1,y1 - big red crossed open dots
plot(x1,y1,type="p",pch=10,cex=2,col="red",
	main="A second axis example",
	ylab="Big values",ylim=c(0,1100),
	xlab="Ordered units")

# Connect the dots with a nice dashed green line
points(x1,y1,type="l",lty=3,lwd=2,col="green")

# Overlay a new plot
par(new=T)

# Draw the second plot x1,y2 - bigger black dots
# Remember axes=FALSE, no box
plot(x1,y2,type="p",pch=20,cex=2,col="black",axes=FALSE,bty="n",
	xlab="",ylab="")

# Connect the dots with a grey dash
points(x1,y2,type="l",lty=2,lwd=2,col="grey")

# Create the second Y axis on the right (4) side
axis(side=4,at=pretty(y2axis))

# Generate an axis label
mtext("Little values",side=4,line=2.5)

# Add a legend - could use locator rather than coordinates
legend(15,0.2,c("Big Y","Little Y"),lty=1,lwd=2,col=c("green","grey"))


#################### END #############################################

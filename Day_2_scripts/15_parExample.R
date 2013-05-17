#################################################################
## 14_parExample.R                                             ##
## Graphics example, getting to grips with PAR                 ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###

# Produce a two by two plot
par(mfrow=c(2,2))

# Set outer margins
# Bottom, left, top, right
par(oma=c(1,0,1,0))

# Set figure margins
par(mar=c(4,2,4,2))

# Set the plot background
par(bg="lightblue",fg="darkgrey")

# Set default plot character and scale factor
par(pch=16,cex=1.4)

# Plot 4 figures, and enclose in a dashed box
plot(1:10)
box("figure",lty=3,lwd=1,col="blue")
plot(1:10)
box("figure",lty=3,lwd=2,col="yellow")
plot(1:10)
box("figure",lty=3,lwd=3,col="black")
plot(1:10)
box("figure",lty=3,lwd=4,col="magenta")

# Draw an outer green box
box("outer",lty=1,lwd=3, col="green")


### END ###

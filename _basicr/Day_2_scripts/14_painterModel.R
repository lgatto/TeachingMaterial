#################################################################
## 13_painterModel.R                                           ##
## Graphics example, demo notion of painter's model in R       ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###

#  Generate some fake plotting data
x <- seq(-2, 2, 0.1)
y <- sin(x)

#  Plot big red points with axes limits
plot(y~x, ylim=c(-1.5,1.5), xlim=c(-2.5,2.5), 
	col="red" ,pch=16, cex=1.4)

#  Overlay thin blue line
lines(y~x, ylim=c(-1.5,1.5), xlim=c(-2.5,2.5), type="l", 
	col="blue", lty=1,lwd=2)

# Draw a white box to obscure bottom half of graph
rect(-2.5,0,2.5,-1.5, col="white", border="white")


### END ###

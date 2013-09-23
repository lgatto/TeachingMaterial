#################################################################
## 16_buildingAplot.R                                          ##
## Graphics example, a fully annotated plot with margin text   ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###
# 1 Figure per page, white background, black objects, scale factor 1
par(mfrow=c(1,1))
par(bg="white",fg="black",cex=1)

# 1 line gap for outer margins
par(oma=c(1,1,1,1))

# Default settings for figure margins.
# Note use 0f +.1 to offset text from axes labels a little bit
par(mar=c(5,4,4,2)+0.1)

# Draw a plot
plot(1:10,main="The plot title", sub="A subtitle", 
	xlab="Numbers",
	ylab="More numbers")

# Add some margin text
mtext(c("Bottom", "Left", "Top", 	"Right"),
		c(1,2,3,4),  	line=.5)

# Add an in plot text
text(2,10,"Text at X=2,Y=10")

# Click to add a legend
legend(locator(1),"Some Legend", fill="red")


### END ###

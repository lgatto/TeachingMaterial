#################################################################
## 19_6PanelPlotscript.R                                       ##
## Graphics exercise.  Generate some high level plots. Display ##
## In various ways.  Get user to write a function              ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###

#Make a function to draw the plots
plotXYBarBox <- function(){
	#Some arbitrary data for plotting
	
	# for x/y plot
	x <- 1:20
	y <- sample(1000, 20) # randomly pick y values from 1:1000

	# for bar plots, a vector of values - 10 random numbers from 1:100
	bar.data <- sample(100, 10)
	
	# for box plots a range of values - normal values with different means
	box.data <- replicate(5, rnorm(100, mean=runif(1)*5))

	# set up the plot with 6 panels
	par(mfrow=c(3,2))

	#Plot x/y plots

	#Default
	plot(x,y)

	#Annotated
	plot(x,y, xlab="Scatter plot X values", ylab="Y values",
		main="Scatter plot", xlim=c(0,10), ylim=c(0,1500), col="blue")
	mtext("This is \n margin text",side=4, cex=0.7, line=1)
	legend("topright",legend="random values",fill="blue")

	#Plot barplots
	#Default
	barplot(bar.data)
	#Annotated
	barplot(bar.data, xlab="Bar plot X values", ylab="Y values",
		main="Bar plot",col=rainbow(5))
	mtext("This is \n margin text",side=4, cex=0.7, line=1)

	#Plot boxplots
	#Default
	boxplot(box.data)
	#Annotated
	boxplot(box.data, xlab="Box plot X categories", ylab="Y values",
		main="Box plot", ylim=c(0,15),col=rainbow(5))
	mtext("This is \n margin text",side=4, cex=0.7, line=1)
	legend("topright",legend=1:5,fill=rainbow(5),ncol=3)
}

#Make the plots
#Plot to screen

plotXYBarBox()

#plot to postscript
postscript(file="sixPanels.ps", paper="a4", width=21, height=27,horizontal=F)
plotXYBarBox()
dev.off()

#plot to jpeg
jpeg(file="sixPanels.jpg", height=800, width=560,res=150)
par(oma=c(1,1,1,1))
plotXYBarBox()
dev.off()

#plot to pdf
pdf(file="sixPanels.pdf",paper="a4", width=21, height=27)
par(oma=c(1,1,1,1))
plotXYBarBox()
dev.off()

###################END SCRIPT######################################

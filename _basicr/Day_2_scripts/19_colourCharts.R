#################################################################
## 18_colourCharts.R                                           ##
## Graphics example, a function which displays colour palettes ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###

## Define a function that generates the plot / two arguments 1) colour chart 2) plot title

paletteMap<-function(myCols, aTitle) {
	nData<-1  # a simple counter that records the index number of the current colour to plot
	x<-0	# coordinates of x
	y<-1    # coordinates of y

	plot(NULL, xlim=c(0,6), ylim=c(0,12),xaxt="n",yaxt="n",ylab="",xlab="",main=aTitle) # sets up an empty plotting chart

	while (nData <= nCols){  # For loop sets x / y plotting coordinates
		if(x<5){
			x<-x+1
		} else {
			x<-1
			y<-y+1
		}
		points(x,y,pch=15, cex=1, col=myCols[nData]) # actually plots the colour (as a big square)
		nData<-nData+1	# increments the colour index number
	}
}
### END FUNC ###

### START PROC ###
nCols<-50	# default number of colours / steps in gradients

#jpeg("HeatMaps.jpg", width=600, height=800, res=150)
par(mfrow=c(3,2))

#Rainbow palette colours
rainCols<-rainbow(nCols) 		# Generate colours
paletteMap(rainCols,"Rainbow Colours")	# Produce plot

#Heatmap palette colours
heatCols<-heat.colors(nCols)
paletteMap(heatCols, "Heat Colours")

#Terrain palette colours
terCols<-terrain.colors(nCols)
paletteMap(terCols, "Terrain Colours")

#Topological palette colours
topoCols<-topo.colors(nCols)
paletteMap(topoCols, "Topological Colours")

#Cyan/Magenta gradient colours
cmCols<-cm.colors(nCols)
paletteMap(cmCols, "Cyan Magenta Colours")

#dev.off()
## END PROC ###



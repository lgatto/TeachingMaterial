#################################################################
## 09_normalTests.R                                            ##
## Exercise in using classical tests                           ##
## Simple normal test                                          ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ##
#  Create data
y<-rnorm(1000) 			# normally distributed data
yy<-exp(rnorm(1000)) 		# exponential data

#  Define a Y scale range, and make plots
yScale<-c(min(y),max(yy)) 	# y axis scale
qqnorm(y, ylim=yScale) 		# plot normal series

par(new=T)			# Overlay a new plot
qqnorm(yy, ylim=yScale, col="blue") # plot exponential series
qqline(y, lty=2, col="black")       # add line of best fit (black for normal)
qqline(yy, lty=2, col="blue")	    # add line of best fit (blue for exponential)

shapiro.test(y)	# Normality test, hypothesis test is false, normal data is not significantly different to normal!!

shapiro.test(yy) # Normality test, hypothesis test is true, exponential data is significantly different to normal!!

### END ###


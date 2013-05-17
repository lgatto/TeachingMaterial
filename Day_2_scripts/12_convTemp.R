#################################################################
## 12_convTemp.R                                               ##
## Exercise user functions with arguments & error checking     ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###
convTemp <- function(t=0, unit="c"){
	# check parameters, stop if non numeric temperature
	if(!is.numeric(t)){
		stop("Non numeric temparture entered")
	}
	
	# check parameters, stop if unrecognised unit code
	if(!(unit %in% c("c","f","k"))){
		stop("Unrecognized temperature unit. \n Enter either (c)entigrade, (f)ahreneinheit or (k)elvin")
	}

	# convert based on input unit

	# If entered centigrade, report farenheit & kelvin
	if(unit=="c"){
		fTemp <- 9/5 * t + 32
		kTemp <- t + 273.15
		# produce a message for the output string
		output <- paste(t,"C is: \n",fTemp,"F \n",kTemp,"K \n")
		# display message on console
		cat(output)
	}

	# If entered farenheit, report centigrade & kelvin
	if(unit=="f"){
		cTemp <- 5/9 * (t-32)
		kTemp <- cTemp + 273.15
		output <- paste(t,"F is: \n",cTemp,"C \n",kTemp,"K \n")
		cat(output)
	}

	# If entered kelvin, report farenheit & centigrade
	if(unit=="k"){
		cTemp <- t - 273.15
		fTemp <- 9/5 * cTemp + 32
		output <- paste(t,"K is: \n",cTemp,"C \n",fTemp,"F \n")
		cat(output)
	}
}

### END ###

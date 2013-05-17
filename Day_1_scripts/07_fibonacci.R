#################################################################
## 07_fibonacci.R                                              ##
## Exercise using loops and conditional branching to generate  ##
## a Fibonacci number series.  Three solutions given, first is ##
## procedural, last two are custom functions                   ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###

# Solution 1 PROCEDURE
# computes first 10 elements using a basic for loop
#
# Define an object to hold results of Fib sequence as it is computed
fib <- rep(0, 10)

# Loop through each element of Fib sequence. i is the loop index counter
for (i in 1:10){
	# if i is greater than positiion 2, then we can start to work out
	# fib sequence at current position by adding together values from 2 former iterms
	# otherwise else simply fills 1 in to lower index positions.
	if ( i > 2 ){
		fib[i] <- fib[i-2] + fib[i-1]
	} else {
		fib[1] <- 1
		fib[2] <- 1
	}
}


##
## NOTE: above is designed to demonstrate loops and branches
# Can be done simpler

fib <- rep(1, 10)

for (i in 3:10) {
	fib[i] <- fib[i-2] + fib[i-1]
}



# Solution 1 as a FUNCTION
# computes first n elements of series, defaults to 10
#
# Simply wrap code in a function statement, substituting '10' for n
# n is the argument of the function. it defaults to 10
# for when no value specificed
fibonacci <- function(n=10){

	fib <- rep(0, n)
	
	for (i in 1:n){
		if ( i > 2 ){
			fib[i] <- fib[i-2]+fib[i-1]
		} else{
			fib[1] <- 1
			fib[2] <- 1
		}
	}

	# We explicitly return the result of the function call
	return(fib)
}

# Solution 2 as a FUNCTION
# computes fibonacci number at position n
# Taken from M Crawley, The R Book

fibonacciAtPosition <- function(n){

	# Three variables are required.  Next variable. Previous variable & a variable to hold the currentVal of fib
	# at each itteration of the while loop.  While loop itterates until n is less than 1

	nextVal <- 1 # nextVal is the upper value in sequence.  i.e. pos 1 at start (n=1)
	prevVal <- 0 # prevVal is the lower value in sequence   i.e. pos 0 at start (n=1)
	# After calculation of current Value, the above will shift once to right
	# Making prevVal equal to nextVal.

	#  Loop through each value of fib up to  n itterations
	while( n > 1 ){
		# currentVal at n = 1 is 1
		currentVal <- nextVal
		# nextVal at n = 1 is 1 + 0 ... 1
		nextVal <- nextVal + prevVal
		# prevVal at n = 1 is 1
		prevVal <- currentVal
		# number of itterations left to do is decreased by 1
		n <- n-1
	}
    return(nextVal)
}

#####
# for n = 5, the values after the loop are:
#
# n     nextVal    prevVal    currentVal
# 5       1          1            1
# 4       2          1            1
# 3       3          2            2
# 2       5          3            3
# 

### END ###

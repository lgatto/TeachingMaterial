#################################################################
## 06_loopExercise.R                                           ##
## Exercise in looping with FOR and While and BREAK            ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###

## For Loop ##
# f is the loop counter, holding values 1 to 26 as appropriate for
# current itteration

for (f in 1:26){

# break line added. if object letters[f] is n ... i.e. f is 14!
	if ( letters[f]=="n") break

# output letter indexed by f to console
	print (letters[f])
}




## While Loop ##
## This one counts backwards

# Set upper limit of index variable 'f' to 26 (letter z)
f <- 26

# Code will loop only when f does not equal 0
while( f!=0 ){

	if (letters[f]=="n") break	# code will break if letters object indexed by f is n i.e f is 14

# Display letter on console	
	print(letters[f])

# Subtract 1 from current value of f
	f <- f -1

} 

### END ###

#################################################################
## 04_Objects.R                                                ##
## Example of defining various object modes and types          ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###

#Integer
i <- sample(20)
i

#Characters
s <- letters[i]
s

#Logical
l<-i >=10
l

#Double precision / floats
d <- runif(20)
d

#Print out mode and type
#Integer
typeof(i)

#Logical
typeof(l)

#Character
typeof(s)

#Double
typeof(d)


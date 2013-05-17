#################################################################
## 05_phoneBook.R                                              ##
## Exercise in defining objects and indexing a dataframe       ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###

# character vector, length 10
firstName<-c("Adam","Eve","John","Mary","Peter","Paul","Luke",
	"Matthew","David","Sally")
length(firstName)

# character vector, length 10
secondName<-c("Tiny","Large","Small","Davis","Thumb","Daniels",
	"Edwards","Smith","Howkins","Dutch")
length(secondName)

# Double vector, length 10
telNumber<-c(111111,222222,333333,444444,555555,666666,777777,
	888888,121212,232323)
length(telNumber)

# Logical vector, length 10
notListed<-c(TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,
	TRUE,FALSE,TRUE,FALSE)
length(notListed)

# Dataframe from vector objects
# Note column names are not defined in function call, hence column names will be taken
# from original vector names
phoneBook<-data.frame(firstName, secondName, paste(firstName,secondName),
	telNumber, notListed, stringsAsFactors=FALSE)

# Name columns in dataframe
# Not really required, could have be done above e.g. ... data.frame(First_Name=firstNam, ...
names(phoneBook)<-c("First_Name","Second_Name","Full_Name","Tel_Number",
	"Not_Listed")

# Display phonebook object
phoneBook

### END ###

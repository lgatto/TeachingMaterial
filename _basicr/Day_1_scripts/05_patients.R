#################################################################
## 05_patients.R                                               ##
## Exercise in defining objects and indexing a dataframe       ##
##                                                             ##
## Ian Roberts.  ir210@cam.ac.uk                               ## 
#################################################################
##
##
### START ###

# character vector
firstName<-c("Adam","Eve","John","Mary","Peter","Paul","Joanna","Matthew","David","Sally")

# character vector
secondName<-c("Jones","Parker","Evans","Davis","Baker","Daniels","Edwards","Smith","Roberts","Wilson")

# Factor
sex<-c("Male","Female","Male","Female","Male","Male","Female","Male","Male","Female")

# Numeric vector
age<-c(50,21,35,45,28,31,42,33,57,62)

# Double vector
weight<-c(70.8,67.9,75.3,61.9,72.4,69.9,63.5,71.5,73.2,64.8)


# Logical vector
consent<-c(TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,TRUE)

# Dataframe from vector objects
# Note column names are not defined in function call, hence column names will be taken
# from original vector names
patients<-data.frame(firstName, secondName, paste(firstName,secondName),factor(sex), age, weight, consent, stringsAsFactors=FALSE)

# Name columns in dataframe
# Not really required, could have be done above e.g. ... data.frame(First_Name=firstNam, ...
names(patients)<-c("First_Name","Second_Name","Full_Name","Sex","Age","Weight","Consent")

# Display data frame
print(patients)

### END ###

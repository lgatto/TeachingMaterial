
## ----, echo=FALSE--------------------------------------------------------
opts_chunk$set(eval=FALSE)


## ----, eval = TRUE-------------------------------------------------------
x
x <- 1
y <- 1:10
z <- "foo"
x
x <- 2
x


## ------------------------------------------------------------------------
sum(1, 2, 3)
sqrt(2)
x <- 2
x^2


## ----, eval=TRUE---------------------------------------------------------
environment()
ls()
rm(x)
ls()


## ------------------------------------------------------------------------
getwd()


## ----, eval=FALSE--------------------------------------------------------
## install.packages("devtools")
## library("devtools")


## ------------------------------------------------------------------------
library("devtools")
install_github("lgatto/camweather")


## ----, eval=FALSE--------------------------------------------------------
## install.packages("camweather_0.1.2.zip", repos = NULL)    ## on windows
## install.packages("camweather_0.1.2.tar.gz", repos = NULL) ## elsewhere
## install.packages("camweather_0.1.2.tar.gz", repos = NULL, type = "src") ## mac


## ------------------------------------------------------------------------
allpkgs <- installed.packages()
nrow(allpkgs)


## ----, eval=FALSE--------------------------------------------------------
## head(allpkgs)


## ----, eval = TRUE-------------------------------------------------------
ls()
rm(list=ls())
ls()


## ----, eval=FALSE--------------------------------------------------------
## packageDescription("devtools")
## help(package = "devtools")


## ----, eval=TRUE---------------------------------------------------------
sessionInfo()
version


## ----, eval = TRUE-------------------------------------------------------
a <- 1L
typeof(a)

b <- 1
typeof(b)

typeof(a + b)


## ----, eval = TRUE-------------------------------------------------------
f <- function(x) { 10 }
system.time(f(Sys.sleep(3)))


## ----, eval = TRUE-------------------------------------------------------
f <- function(x) { force(x); 10 }
system.time(f(Sys.sleep(3)))


## ------------------------------------------------------------------------
vector()
matrix()
array()
list()
data.frame()


## ------------------------------------------------------------------------
x <- 1
length(x)
x <- c(1, 2, 4) ## combine values into a vector
length(x)
x[2] ## second element of x


## ------------------------------------------------------------------------
names(x)
names(x) <- c("a", "b", "c")
x["b"]


## ------------------------------------------------------------------------
x[2] <- x[2]*10
x


## ------------------------------------------------------------------------
y <- x[-2]
y


## ------------------------------------------------------------------------
typeof(c(1, 10, 1))
typeof(c(1L, 10L, 1L))
typeof(c(1L, 10, 1.0))
typeof(letters)
typeof(c(TRUE, FALSE)) ## TRUE and FALSE are reserved words


## ----, eval=TRUE---------------------------------------------------------
gender <- factor(c("male", "female", "male"))
gender
class(gender)
typeof(gender)


## ------------------------------------------------------------------------
character()
logical()
numeric()
double()
factor()


## ----, eval = TRUE-------------------------------------------------------
x <- c(1, 2, 3)
y <- c(3, 2, 1)
x + y
x^2
sqrt(x)


## ----, eval = TRUE-------------------------------------------------------
x <- c("a", "b", "c")
paste(x, 1:3, sep = ".")


## ----, eval  = TRUE------------------------------------------------------
x <- c(1, 2, 3, 4)
y <- c(1, 2)
x + y
z <- c(1, 2, 3)
x + z


## ----, eval = TRUE-------------------------------------------------------
x <- c(1, 2, 3, 4)
x[2:3] <- 10
x


## ----, eval = TRUE-------------------------------------------------------
x <- c("a", "b", "c")
paste(x, 1, sep = ".")


## ------------------------------------------------------------------------
seq(1, 10, 2)
seq(1, 3, length = 7)


## ------------------------------------------------------------------------
1:10


## ------------------------------------------------------------------------
summary(rnorm(100))
runif(10)


## ------------------------------------------------------------------------
sample(LETTERS[1:5])
sample(LETTERS[1:5], 10, replace = TRUE)


## ----, eval=TRUE---------------------------------------------------------
m <- c(1, 2, 3, 4, 5, 6)
dim(m) <- c(2, 3) ## always rows first
m
class(m) ## a matrix
mode(m)  ## of numerics


## ------------------------------------------------------------------------
m <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
m <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2, nrow = 3)


## ------------------------------------------------------------------------
m[1:2, 1:2]
m[, 2:3]
m[1, 1] <- 10
m
m[, -1]


## ------------------------------------------------------------------------
colnames(m) <- letters[1:3]
rownames(m) <- LETTERS[1:2]
## or provideDimnames(m, base = list(letters, LETTERS))
m
m["A", "c"]


## ----, tidy=FALSE--------------------------------------------------------
M <- matrix(c(4, 8, 5, 6, 4, 2, 1, 5, 7), nrow=3)
dimnames(M) <- list(year =
                    c(2005, 2006, 2007),
                    "mode of transport" =
                    c("plane", "bus", "boat"))


## ------------------------------------------------------------------------
m[-1, ] ## becomes a vector
m[-1, , drop = FALSE] ## remains a matrix
m[, -(1:2)] ## becomes a vector
m[, -(1:2), drop = FALSE] ## remains a matrix


## ------------------------------------------------------------------------
x <- 1:6
t(x) ## becomes a matrix
dim(t(x)) 
dim(t(t(x)))


## ------------------------------------------------------------------------
x <- 1:30
dim(x) <- c(5, 3, 2)
## or array(1:30, dim = c(5, 3, 2))
x


## ------------------------------------------------------------------------
l <- list(M = matrix(1:10, ncol = 5), V = 1:2, F = sum)
l
l$M + l[["V"]]
l$F(l$M) ## same as sum(1:10)


## ------------------------------------------------------------------------
class(l[2])
class(l[[2]])


## ------------------------------------------------------------------------
class(l[2:3])
class(l[[2:3]])


## ----, tidy=FALSE--------------------------------------------------------
age <- c(50, 21, 35, 45, 28,
         31, 42, 33, 57, 62)
weight <- c(70.8, 67.9, 75.3, 61.9, 72.4,
            69.9, 63.5, 71.5, 73.2, 64.8)
firstName <- c("Adam", "Eve", "John", "Mary",
               "Peter", "Paul", "Joanna", "Matthew",
               "David", "Sally")
secondName <- c("Jones", "Parker", "Evans",
                "Davis", "Baker", "Daniels",
                "Edwards", "Smith", "Roberts", "Wilson")
consent <- c(TRUE, TRUE, FALSE, TRUE, FALSE,
             FALSE, FALSE, TRUE, FALSE, TRUE)
gender <- c("Male", "Female", "Male", "Female",
            "Male", "Male", "Female", "Male",
            "Male", "Female")
patients <- data.frame(firstName, secondName,
                       gender = factor(gender),
                       age, weight, consent,
                       stringsAsFactors=FALSE)
rm(age, consent, firstName, gender, secondName, weight)


## ------------------------------------------------------------------------
patients[patients$consent, ]
patients[patients$gender == "Male", ]
patients[patients$gender == "Male" & patients$consent, ]


## ------------------------------------------------------------------------
avgw <- mean(patients$weight)
patients[patients$weight < avgw, ]


## ------------------------------------------------------------------------
patients[order(patients$age), ]


## ------------------------------------------------------------------------
sum(replicate(100, system.time(mtcars["Volvo 142E", "carb"])["elapsed"]))
sum(replicate(100, system.time(mtcars[32, 11])["elapsed"]))
sum(replicate(100, system.time(mtcars[[11]][32])["elapsed"]))
sum(replicate(100, system.time(mtcars$carb[32])["elapsed"]))


## ------------------------------------------------------------------------
m1 <- matrix(1:12, ncol = 3)
m2 <- matrix(1:9, ncol = 3)
m3 <- matrix(1:16, ncol = 4)
cbind(m1, m3)
rbind(m1, m2)


## ------------------------------------------------------------------------
cbind(m1, m2)
rbind(m1, m3)


## ----, eval=FALSE--------------------------------------------------------
## chrsms13 <- weatherdata("2013-12-25")
## chrsms12 <- weatherdata("2012-12-25")
## chrsms11 <- weatherdata("2011-12-25")
## chrsms <- rbind(chrsms11, chrsms12, chrsms13)
## nrow(chrsms11) + nrow(chrsms12) + nrow(chrsms13)
## nrow(chrsms)


## ------------------------------------------------------------------------
class(NA)
class(NaN)
class(NULL)


## ------------------------------------------------------------------------
l[[2]] <- NULL
l
length(NULL)
c(1, NULL)
list(1, NULL)


## ------------------------------------------------------------------------
class(as.integer(NA))
class(as.numeric(NA))
class(as.character(NA))


## ------------------------------------------------------------------------
sum(1, 2, NA)
sum(1, 2, NA, na.rm = TRUE)


## ------------------------------------------------------------------------
as.numeric("1")
as.character(1+2)
as.integer(1.9) ## see also floor and ceiling
as.numeric("1a")


## ----, eval = TRUE-------------------------------------------------------
e <- new.env()
e
e$a <- 1
e$a
ls() ## list content of global environment
ls(e) ## list content of e
a <- 10 ## a different variable a
e$a
e[["a"]]


## ----, eval = TRUE-------------------------------------------------------
lockEnvironment(e)
e$b <- 10
e$a <- 10
lockBinding("a", e)
e$a <- 100


## ----, eval=TRUE---------------------------------------------------------
class(x <- rnorm(100))
y <- rnorm(100)
class(1L)
class("123")
class(sum)
class(model <- lm(y ~ x))
str(model)


## ----, eval=TRUE---------------------------------------------------------
library("affydata")
data("Dilution")
class(Dilution)
str(Dilution)



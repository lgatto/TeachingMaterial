### R code from vignette source 'Sec-Intro.Rnw'

###################################################
### code chunk number 1: intro
###################################################
x <- 1 ## a variable
x
x = 2 ## overwrite the value x
x
y <- length(x) ## calling a function
y


###################################################
### code chunk number 2: env
###################################################
ls()
rm(y)
ls()


### R code from vignette source 'Sec-DataTypes.Rnw'

###################################################
### code chunk number 1: env
###################################################
library(Biobase)
data(sample.ExpressionSet)


###################################################
### code chunk number 2: vec0
###################################################
c(1,3,9,-1)


###################################################
### code chunk number 3: vecs
###################################################
mode(1)
typeof(1)
mode(1L)
typeof(1L)


###################################################
### code chunk number 4: vecs2
###################################################
mode("1")
typeof("1")
mode(TRUE)
typeof(FALSE)
## as we are talking about booleans...
TRUE & TRUE
TRUE | FALSE


###################################################
### code chunk number 5: isas
###################################################
x <- 1
typeof(x)
y <- as.integer(x)
typeof(y)
is.integer(y)


###################################################
### code chunk number 6: specialvals (eval = FALSE)
###################################################
## NULL; NA; NaN; Inf; -Inf
## is.null(); is.na(); is.infinite()


###################################################
### code chunk number 7: class
###################################################
class(x)
class("a character")


###################################################
### code chunk number 8: createvecs0
###################################################
vector(mode = "character", length = 3)
vector(mode = "numeric", length = 4)
numeric(4)


###################################################
### code chunk number 9: createvecs
###################################################
x <- c(1, 4, 7, 10) ## concatenate
x
y <- 1:5 ## integer sequence 
y
z <- seq(from = 1, to = 10, by = 2)
z


###################################################
### code chunk number 10: args
###################################################
z1 <- seq(from = 1, to = 10, by = 2)
z2 <- seq(1, 10, 2)
z1 == z2
all(z1 == z2)
identical(z1, z2)


###################################################
### code chunk number 11: vec
###################################################
x <- 1:5; y <- 5:1
x
y
x + y
x^2


###################################################
### code chunk number 12: mat
###################################################
m <- matrix(1:12, nrow = 4, ncol = 3)
m
dim(m)


###################################################
### code chunk number 13: mat2
###################################################
matrix(1:11, 4, 3)
matrix(1:12, 3, 3)


###################################################
### code chunk number 14: mat
###################################################
x <- 1:12
class(x)
dim(x)
dim(x) <- c(4, 3)
x
class(x)


###################################################
### code chunk number 15: arrays
###################################################
array(1:16, dim = c(2, 4, 2))


###################################################
### code chunk number 16: list
###################################################
ll <- list(a = 1:3, c = length)
ll
ll[[1]]
ll$c(ll)


###################################################
### code chunk number 17: dfr
###################################################
dfr <- data.frame(type = c(
                    rep("case", 2), 
                    rep("ctrl", 2)),
                  time = rnorm(4))
dfr


###################################################
### code chunk number 18: dfr2
###################################################
dfr[1,]
dfr[1, "time"]
dfr$time


###################################################
### code chunk number 19: env
###################################################
e <- new.env()
e[["a"]] <- 1:3
assign("b", "CSAMA", envir = e)
ls(e)
e$a
get("b", e)


###################################################
### code chunk number 20: names
###################################################
x <- c(a = 1, b = 2)
x
names(x)


###################################################
### code chunk number 21: Sec-DataTypes.Rnw:232-238
###################################################
M <- matrix(c(4, 8, 5, 6, 4, 2, 1, 5, 7), nrow=3)
dimnames(M) <- list(year = 
                    c(2005, 2006, 2007),
                    "mode of transport" =                     
                    c("plane", "bus", "boat"))
M


###################################################
### code chunk number 22: Sec-DataTypes.Rnw:246-247
###################################################
sample.ExpressionSet$type


###################################################
### code chunk number 23: eset
###################################################
library(Biobase)
data(sample.ExpressionSet)
sample.ExpressionSet


###################################################
### code chunk number 24: eset2
###################################################
class(sample.ExpressionSet)
slotNames(sample.ExpressionSet)


###################################################
### code chunk number 25: eset2b (eval = FALSE)
###################################################
## class?ExpressionSet


###################################################
### code chunk number 26: eset3
###################################################
exprs(sample.ExpressionSet)[1:4, 1:3]
dim(exprs(sample.ExpressionSet)) ## or dim(sample.ExpressionSet)


###################################################
### code chunk number 27: eset3
###################################################
phenoData(sample.ExpressionSet)


###################################################
### code chunk number 28: eset3
###################################################
featureData(sample.ExpressionSet)


###################################################
### code chunk number 29: pdata
###################################################
head(pData(sample.ExpressionSet))


### R code from vignette source 'Sec-Manip.Rnw'

###################################################
### code chunk number 1: env
###################################################
library(Biobase)
data(sample.ExpressionSet)


###################################################
### code chunk number 2: sub1
###################################################
x <- 1:10
x[3:7]
x[9:11]
x[0:1]
x[c(1, 7, 2, NA)]


###################################################
### code chunk number 3: sub2
###################################################
x[2] <- 20
x[4:5] <- x[4:5] * 100
x[1:6]


###################################################
### code chunk number 4: sub3
###################################################
x <- 1:10
x[-c(3:7)]


###################################################
### code chunk number 5: sub3
###################################################
x[c(TRUE, TRUE, rep(FALSE, 8))]
x > 5
x[x > 5]
x[c(TRUE, FALSE)] ## recycled


###################################################
### code chunk number 6: sub4
###################################################
x <- c(a = 1, b = 2, c = 2)
x[c("a", "c")]
x[c("a", "d")]


###################################################
### code chunk number 7: submat
###################################################
M <- matrix(1:12, 3)
M[2,3] <- 0
M


###################################################
### code chunk number 8: submat
###################################################
M < 9
M[M < 9] <- -1
M


###################################################
### code chunk number 9: sublist
###################################################
ll <- list(a = 1:3, b = "CSAMA", c = length)
ll[1] ## still a list
ll[[1]] ## first element of the list


###################################################
### code chunk number 10: subexprs
###################################################
sample.ExpressionSet[1:10, 1:2]


### R code from vignette source 'Sec-Useful.Rnw'

###################################################
### code chunk number 1: read.csv
###################################################
read.table("./Data/data.csv", sep = ",",
           header = TRUE, row.names = 1)


###################################################
### code chunk number 2: read.csv
###################################################
read.csv("./Data/data.csv", row.names = 1)


###################################################
### code chunk number 3: read.csv
###################################################
x <- read.csv("./Data/data.csv", row.names = 1)
save(x, file = "./Data/data.rda")
rm(x)
load("./Data/data.rda")
x[1:3, ]


###################################################
### code chunk number 4: str1
###################################################
paste("abc", "def", sep = "-")
paste0("abc", "def")


###################################################
### code chunk number 5: str2
###################################################
month.name[1:4]
grep("Feb", month.name)
grep("Feb", month.name, value = TRUE)
grepl("Feb", month.name)


###################################################
### code chunk number 6: str3
###################################################
month.name[1]
length(month.name[1])
nchar(month.name[1])


###################################################
### code chunk number 7: str4
###################################################
strsplit("abc-def", "-")


###################################################
### code chunk number 8: comp
###################################################
set.seed(1)
x <- sample(letters[1:10], 6)
y <- sample(letters[1:10], 6)
x
y


###################################################
### code chunk number 9: comp2
###################################################
intersect(x, y)
setdiff(x, y)
union(x, y)


###################################################
### code chunk number 10: comp3
###################################################
x %in% y
x == y
match(x, y)


###################################################
### code chunk number 11: gen
###################################################
seq(1,7,3)
rep(1:2, 2)
rep(1:2, each = 2)


###################################################
### code chunk number 12: gen2
###################################################
runif(5)
rnorm(5)


###################################################
### code chunk number 13: aboutdata
###################################################
table(sample(letters, 100, replace = TRUE))
summary(rnorm(100))
head(x)
tail(x)


###################################################
### code chunk number 14: head
###################################################
M <- matrix(rnorm(1000), ncol=4)
head(M)


### R code from vignette source 'Sec-Plotting.Rnw'

###################################################
### code chunk number 1: env
###################################################
library(Biobase)
data(sample.ExpressionSet)


###################################################
### code chunk number 2: plotcode
###################################################
plot(exprs(sample.ExpressionSet[, 1]), 
     exprs(sample.ExpressionSet[, 2]), 
     log = "xy", 
     xlab = sampleNames(sample.ExpressionSet)[1],
     ylab = sampleNames(sample.ExpressionSet)[2])
abline(0, 1)
grid()


###################################################
### code chunk number 3: plotfig1
###################################################
plot(exprs(sample.ExpressionSet[, 1]), 
     exprs(sample.ExpressionSet[, 2]), 
     log = "xy", 
     xlab = sampleNames(sample.ExpressionSet)[1],
     ylab = sampleNames(sample.ExpressionSet)[2])
abline(0, 1)
grid()


###################################################
### code chunk number 4: plotcode
###################################################
pairs(log2(exprs(sample.ExpressionSet)[, 1:4]),
      pch = 19,
      col = "#0000FF20")


###################################################
### code chunk number 5: plotfig2
###################################################
pairs(log2(exprs(sample.ExpressionSet)[, 1:4]),
      pch = 19,
      col = "#0000FF20")


###################################################
### code chunk number 6: plotcode
###################################################
boxplot(log2(exprs(sample.ExpressionSet)))


###################################################
### code chunk number 7: plotfig3
###################################################
boxplot(log2(exprs(sample.ExpressionSet)))


###################################################
### code chunk number 8: plotcode
###################################################
smoothScatter(log2(exprs(sample.ExpressionSet)[, 1:2]))


###################################################
### code chunk number 9: plotfig4
###################################################
smoothScatter(log2(exprs(sample.ExpressionSet)[, 1:2]))


### R code from vignette source 'Sec-Programming.Rnw'

###################################################
### code chunk number 1: flow (eval = FALSE)
###################################################
## for (var in seq) expr
## while (cond) expr
## repeat expr
## break


###################################################
### code chunk number 2: for
###################################################
for (i in 1:4) { ## bad
  print(i^2)
}
(1:4)^2 ## good


###################################################
### code chunk number 3: apply
###################################################
M <- matrix(1:9, ncol = 3)
M
apply(M, 1, max)
apply(M, 2, max)
2


###################################################
### code chunk number 4: slapply
###################################################
sapply(month.name[1:2], paste0, "_2012")
lapply(month.name[1:2], paste0, "_2012")


###################################################
### code chunk number 5: replicate
###################################################
mean(rnorm(100))
replicate(3, mean(rnorm(100)))
replicate(2, rnorm(3))


###################################################
### code chunk number 6: flow (eval = FALSE)
###################################################
## if (cond) expr1 else expr2
## ifelse(cond, expr1, expr2)
## switch


###################################################
### code chunk number 7: ifelse
###################################################
x <- 2
if (x > 0) { ## bad
  log2(x)
} else {
  log2(-x)
}
log2(abs(x)) ## better 


###################################################
### code chunk number 8: fun (eval = FALSE)
###################################################
## myFun <- function(param1, param2, ...) {
##   ## function body
##   ## acting on copies of the params
##   ans <- param1 + param2  
##   return(ans) 
## }


###################################################
### code chunk number 9: f1
###################################################
x <- 1
f <- function(x) { x <- x + 10; x }
f(x)
x


###################################################
### code chunk number 10: f2
###################################################
x <- 1
f <- function() { x <- x + 10; x }
f()
x


###################################################
### code chunk number 11: anonymfun
###################################################
M <- matrix(rnorm(50), ncol = 5)
M[sample(50, 10)] <- NA
sum(is.na(M))
apply(M, 1, function(x) sum(is.na(x)))
apply(M, 2, function(x) sum(is.na(x)))


### R code from vignette source 'Sec-Packages.Rnw'

###################################################
### code chunk number 1: bioclite (eval = FALSE)
###################################################
## source("http://www.bioconductor.org/biocLite.R")
## ## or, if you have already done so in the past
## library("BiocInstaller")
## biocLite("packageName")  


###################################################
### code chunk number 2: pckhelp (eval = FALSE)
###################################################
## help(package = "Biobase")


###################################################
### code chunk number 3: pckvig (eval = FALSE)
###################################################
## vignette(package = "Biobase")
## vignette("Bioconductor", package = "Biobase")


###################################################
### code chunk number 4: pckdemo (eval = FALSE)
###################################################
## demo("lattice", package = "lattice")


###################################################
### code chunk number 5: pkgversion
###################################################
packageDescription("Biobase")


### R code from vignette source 'Sec-RBioc.Rnw'


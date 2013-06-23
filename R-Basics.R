
## @knitr intro
x <- 1 ## a variable
x
x = 2 ## overwrite the value x
x
y <- length(x) ## calling a function
y


## @knitr env
ls()
rm(y)
ls()



## @knitr env0, echo=FALSE, message = FALSE
library(Biobase)
data(sample.ExpressionSet)


## @knitr vec0
c(1,3,9,-1)


## @knitr vecs
mode(1)
typeof(1)
mode(1L)
typeof(1L)


## @knitr vecs2
mode("1")
typeof("1")
mode(TRUE)
typeof(FALSE)
## as we are talking about booleans...
TRUE & TRUE
TRUE | FALSE


## @knitr isas
x <- 1
typeof(x)
y <- as.integer(x)
typeof(y)
is.integer(y)


## @knitr specialvals, eval=FALSE
## NULL; NA; NaN; Inf; -Inf
## is.null(); is.na(); is.infinite()


## @knitr class
class(x)
class("a character")


## @knitr createvecs0
vector(mode = "character", length = 3)
vector(mode = "numeric", length = 4)
numeric(4)


## @knitr createvecs
x <- c(1, 4, 7, 10) ## concatenate
x
y <- 1:5 ## integer sequence 
y
z <- seq(from = 1, to = 10, by = 2)
z


## @knitr args
z1 <- seq(from = 1, to = 10, by = 2)
z2 <- seq(1, 10, 2)
z1 == z2
all(z1 == z2)
identical(z1, z2)


## @knitr vec
x <- 1:5; y <- 5:1
x
y
x + y
x^2


## @knitr mat
m <- matrix(1:12, nrow = 4, ncol = 3)
m
dim(m)


## @knitr mat2
matrix(1:11, 4, 3)
matrix(1:12, 3, 3)


## @knitr mat3
x <- 1:12
class(x)
dim(x)
dim(x) <- c(4, 3)
x
class(x)


## @knitr arrays
array(1:16, dim = c(2, 4, 2))


## @knitr list
ll <- list(a = 1:3, c = length)
ll
ll[[1]]
ll$c(ll)


## @knitr dfr, tidy = FALSE
dfr <- data.frame(type = c(
                    rep("case", 2), 
                    rep("ctrl", 2)),
                  time = rnorm(4))
dfr


## @knitr dfr2
dfr[1,]
dfr[1, "time"]
dfr$time


## @knitr env
e <- new.env()
e[["a"]] <- 1:3
assign("b", "CSAMA", envir = e)
ls(e)
e$a
get("b", e)


## @knitr names
x <- c(a = 1, b = 2)
x
names(x)


## @knitr matdimnames, tidy = FALSE
M <- matrix(c(4, 8, 5, 6, 4, 2, 1, 5, 7), nrow=3)
dimnames(M) <- list(year = 
                    c(2005, 2006, 2007),
                    "mode of transport" =                     
                    c("plane", "bus", "boat"))
M


## @knitr 
sample.ExpressionSet$type


## @knitr eset
library(Biobase)
data(sample.ExpressionSet)
sample.ExpressionSet


## @knitr eset2
class(sample.ExpressionSet)
slotNames(sample.ExpressionSet)


## @knitr eset2b, eval=FALSE, tidy=FALSE
## ?ExpressionSet


## @knitr eset3
exprs(sample.ExpressionSet)[1:4, 1:3]
dim(sample.ExpressionSet)


## @knitr eset4
phenoData(sample.ExpressionSet)


## @knitr eset5
featureData(sample.ExpressionSet)


## @knitr pdata
head(pData(sample.ExpressionSet))



## @knitr env, echo=FALSE, message = FALSE
library(Biobase)
data(sample.ExpressionSet)


## @knitr sub1
x <- 1:10
x[3:7]
x[9:11]
x[0:1]
x[c(1, 7, 2, NA)]


## @knitr sub2
x[2] <- 20
x[4:5] <- x[4:5] * 100
x[1:6]


## @knitr sub3
x <- 1:10
x[-c(3:7)]


## @knitr sub4
x[c(TRUE, TRUE, rep(FALSE, 8))]
x > 5
x[x > 5]
x[c(TRUE, FALSE)] ## recycled


## @knitr sub5
x <- c(a = 1, b = 2, c = 2)
x[c("a", "c")]
x[c("a", "d")]


## @knitr submat0
M <- matrix(1:12, 3)
M[2,3] <- 0
M


## @knitr submat1
M < 9
M[M < 9] <- -1
M


## @knitr sublist
ll <- list(a = 1:3, b = "CSAMA", c = length)
ll[1] ## still a list
ll[[1]] ## first element of the list


## @knitr subexprs
sample.ExpressionSet[1:10, 1:2]



## @knitr read.csv0, tidy = FALSE
read.table("./Data/data.csv", sep = ",",
           header = TRUE, row.names = 1)


## @knitr read.csv1
read.csv("./Data/data.csv", row.names = 1)


## @knitr read.csv2
x <- read.csv("./Data/data.csv", row.names = 1)
save(x, file = "./Data/data.rda")
rm(x)
load("./Data/data.rda")
x[1:3, ]


## @knitr str1
paste("abc", "def", sep = "-")
paste0("abc", "def")


## @knitr str2
month.name[1:4]
grep("Feb", month.name)
grep("Feb", month.name, value = TRUE)
grepl("Feb", month.name)


## @knitr str3
month.name[1]
length(month.name[1])
nchar(month.name[1])


## @knitr str4
strsplit("abc-def", "-")


## @knitr comp
set.seed(1)
x <- sample(letters[1:10], 6)
y <- sample(letters[1:10], 6)
x
y


## @knitr comp2
intersect(x, y)
setdiff(x, y)
union(x, y)


## @knitr comp3
x %in% y
x == y
match(x, y)


## @knitr gen
seq(1,7,3)
rep(1:2, 2)
rep(1:2, each = 2)


## @knitr gen2
runif(5)
rnorm(5)


## @knitr aboutdata
table(sample(letters, 100, replace = TRUE))
summary(rnorm(100))
head(x)
tail(x)


## @knitr head
M <- matrix(rnorm(1000), ncol=4)
head(M)



## @knitr env, echo=FALSE, massage = FALSE
suppressPackageStartupMessages(library(Biobase))
data(sample.ExpressionSet)


## @knitr plotcode, eval = FALSE, tidy = FALSE
## plot(exprs(sample.ExpressionSet[, 1]),
##      exprs(sample.ExpressionSet[, 2]),
##      log = "xy",
##      xlab = sampleNames(sample.ExpressionSet)[1],
##      ylab = sampleNames(sample.ExpressionSet)[2])
## abline(0, 1)
## grid()


## @knitr label=plotfig1, echo=FALSE, fig.width=5, fig.height=4, tidy=FALSE, warning = FALSE
plot(exprs(sample.ExpressionSet[, 1]), 
     exprs(sample.ExpressionSet[, 2]), 
     log = "xy", 
     xlab = sampleNames(sample.ExpressionSet)[1],
     ylab = sampleNames(sample.ExpressionSet)[2])
abline(0, 1)
grid()


## @knitr plotcode0, eval = FALSE, tidy = FALSE
## pairs(log2(exprs(sample.ExpressionSet)[, 1:4]),
##       pch = 19,
##       col = "#0000FF20")


## @knitr label=plotfig2,echo=FALSE,fig.width=3.5,fig.height=3.5,tidy=FALSE, warning = FALSE
pairs(log2(exprs(sample.ExpressionSet)[, 1:3]),
      pch = 19,
      col = "#0000FF20")


## @knitr plotcode1, eval = FALSE
## boxplot(log2(exprs(sample.ExpressionSet)))


## @knitr label=plotfig3,echo=FALSE,fig.width=5,fig.height=4,tidy=FALSE
boxplot(log2(exprs(sample.ExpressionSet)))


## @knitr label=plotfig4,echo=TRUE,fig.width=5,fig.height=4,tidy=FALSE, warning = FALSE, message = FALSE
smoothScatter(log2(exprs(sample.ExpressionSet)[, 1:2]))



## @knitr flowctrl, eval = FALSE, tidy = FALSE
## for (var in seq) expr
## 
## while (cond) expr
## 
## repeat expr
## break
## 


## @knitr for
for (i in 1:4) { ## bad
  print(i^2)
}
(1:4)^2 ## good


## @knitr apply
M <- matrix(1:9, ncol = 3)
M
apply(M, 1, max)
apply(M, 2, max)
2


## @knitr slapply
sapply(month.name[1:2], paste0, "_2012")
lapply(month.name[1:2], paste0, "_2012")


## @knitr replicate
mean(rnorm(100))
replicate(3, mean(rnorm(100)))
replicate(2, rnorm(3))


## @knitr flow, eval = FALSE, tidy = FALSE
## if (cond) expr1 else expr2
## 
## ifelse(cond, expr1, expr2)
## 
## switch


## @knitr ifelse
x <- 2
if (x > 0) { ## bad
  log2(x)
} else {
  log2(-x)
}
log2(abs(x)) ## better 


## @knitr fun, eval=FALSE, tidy = FALSE
## myFun <- function(param1, param2, ...) {
##   ## function body
##   ## acting on copies of the params
##   ans <- param1 + param2
##   return(ans)
## }


## @knitr f1
x <- 1
f <- function(x) { x <- x + 10; x }
f(x)
x


## @knitr f2
x <- 1
f <- function() { x <- x + 10; x }
f()
x


## @knitr anonymfun
M <- matrix(rnorm(50), ncol = 5)
M[sample(50, 10)] <- NA
sum(is.na(M))
apply(M, 1, function(x) sum(is.na(x)))
apply(M, 2, function(x) sum(is.na(x)))



## @knitr bioclite, echo=TRUE, eval=FALSE
## source("http://www.bioconductor.org/biocLite.R")
## ## or, if you have already done so in the past
## library("BiocInstaller")
## biocLite("packageName")


## @knitr pckhelp, eval = FALSE
## help(package = "Biobase")


## @knitr pckvig0, eval = FALSE, tidy = FALSE
## vignette(package = "Biobase")


## @knitr pckvig1, eval = FALSE, tidy = FALSE
## vignette("Bioconductor", package = "Biobase")


## @knitr pckdemo, eval = FALSE
## demo("lattice", package = "lattice")


## @knitr pkgversion
packageDescription("Biobase")




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




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


## @knitr eset2b, eval=FALSE
## class?ExpressionSet


## @knitr eset3
exprs(sample.ExpressionSet)[1:4, 1:3]
dim(sample.ExpressionSet)


## @knitr eset4
phenoData(sample.ExpressionSet)


## @knitr eset5
featureData(sample.ExpressionSet)


## @knitr pdata
head(pData(sample.ExpressionSet))



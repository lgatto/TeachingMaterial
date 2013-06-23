
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



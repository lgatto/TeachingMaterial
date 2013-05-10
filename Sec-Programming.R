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



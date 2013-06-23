
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


## @knitr slapply
sapply(month.name[1:2], paste0, "_2012")
lapply(month.name[1:2], paste0, "_2012")


## @knitr apply
M <- matrix(1:9, ncol = 3)
M
apply(M, 1, sum) ## better rowSums
apply(M, 2, sum) ## better colSums
2


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



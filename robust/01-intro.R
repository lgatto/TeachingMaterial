
library("formatR")
tidy_source(text = "a=1+1;a  # print the value
                    matrix ( rnorm(10),5)",
            arrow = TRUE)


head(cars)
head(cars[, 1])
head(cars[, 1, drop = FALSE])


e <- new.env()
e$x <- 1
f <- function(myenv) myenv$x <- 2
f(e)
e$x


library("Biobase")
getClass("eSet")
getClass("AssayData")
new("ExpressionSet")


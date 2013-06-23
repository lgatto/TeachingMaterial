## @knitr env, echo=FALSE, message = FALSE
options(width = 50)


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


## @knitr str4b
strsplit(c("abc-def", "ghi-jkl"), "-")


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


## @knitr aboutdata, size="scriptsize"
table(sample(letters, 100, replace = TRUE))
summary(rnorm(100))
head(x)
tail(x)


## @knitr head
M <- matrix(rnorm(1000), ncol=4)
head(M)



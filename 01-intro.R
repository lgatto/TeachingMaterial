
library("formatR")
tidy_source(text = "a=1+1;a  # print the value
                    matrix ( rnorm(10),5)",
            arrow = TRUE)


head(cars)
head(cars[, 1])
head(cars[, 1, drop = FALSE])


df1 <- data.frame(x = 1:3, y = LETTERS[1:3])
sapply(df1, class)
df2 <- data.frame(x = 1:3, y = Sys.time() + 1:3)
sapply(df2, class)


lapply(df1, class)
lapply(df2, class)


vapply(df1, class, "1")
vapply(df2, class, "1")


x <- 1
f <- function(x) {
    x <- 2
    x
}
x
f(x)
x


e <- new.env()
e$a <- 1
e$b <- LETTERS[1:5]
e$c <- TRUE
e$d <- mean


e$a <- e$b
e$a <- LETTERS[1:5]


e <- new.env()
parent.env(e)


environment()
parent.env(globalenv())
parent.env(parent.env(globalenv()))


globalenv()
emptyenv()
baseenv()


search()
as.environment("package:stats")


ls() ## default is R_GlobalEnv
ls(envir = e)
ls(pos = 1)


search()


e1 <- new.env()
e2 <- new.env()
e1$a <- 1:10
e2$a <- e1$a


e3 <- e1
e3
e1
identical(e1, e3)


e <- new.env()
e$a <- 1
e$b <- 2  ## add
e$a <- 10 ## modify


lockEnvironment(e)
e$k <- 1
e$a <- 100


lockBinding("a", e)
e$a <- 10
e$b <- 10

lockEnvironment(e, bindings = TRUE)
e$b <- 1


e <- new.env()
e$foo <- 1
bar <- 2
where("foo")
where("bar")
where("foo", env = e)
where("bar", env = e)


search()
mean <- function(x) cat("The mean is", sum(x)/length(x), "\n")
mean(1:10)
base::mean(1:10)
rm(mean)
mean(1:10)


library("fortunes")
fortune(174)


rm(list = ls())
x
f1 <- function() x <<- 1
f1()
x


f2 <- function() x <<- 2
f2()
x


f3 <- function() x <- 10
f3()
x


f4 <- function(x) x <-10
f4(x)
x


modify <- function(x) {
	x$a <- 2
	invisible(TRUE)
}


x_l <- list(a = 1)
modify(x_l)
x_l$a


x_e <- new.env()
x_e$a <- 1
modify(x_e)
x_e$a


e <- new.env()
e$a <- 1
e
parent.env(e)

parent.env(e) <- emptyenv()
parent.env(e)
e


e <- new.env(parent.env = empty.env())


x <- 1
e1 <- new.env()
get("x", envir = e1)


get("x", envir = e1, inherits = FALSE)


e2 <- new.env(parent = emptyenv())
get("x", envir = e2)


get("x", envir = e1, inherits = FALSE)


e <- new.env()
e$x <- 1
f <- function(myenv) myenv$x <- 2
f(e)
e$x


library("Biobase")
getClass("eSet")
getClass("AssayData")
new("ExpressionSet")


.pRolocEnv <- new.env(parent=emptyenv(), hash=TRUE)

stockcol <- c("#E41A1C", "#377EB8", "#238B45", "#FF7F00", "#FFD700", "#333333",
              "#00CED1", "#A65628", "#F781BF", "#984EA3", "#9ACD32", "#B0C4DE",
              "#00008A", "#8B795E", "#FDAE6B", "#66C2A5", "#276419", "#CD8C95",
              "#6A51A3", "#EEAD0E", "#0000FF", "#9ACD32", "#CD6090", "#CD5B45",
              "#8E0152", "#808000", "#67000D", "#3F007D", "#6BAED6", "#FC9272")

assign("stockcol", stockcol, envir = .pRolocEnv)

getStockcol <- function() get("stockcol", envir = .pRolocEnv)

setStockcol <- function(cols) {
    if (is.null(cols)) {
        assign("stockcol", stockcol, envir = .pRolocEnv)
    } else {
		assign("stockcol", cols, envir = .pRolocEnv)
	}
}


...
if (missing(col))
  col <- getStockcol()
...


setStockcol <- function(cols) {
	prevcols <- getStockcol()
    if (is.null(cols)) {
        assign("stockcol", stockcol, envir = .pRolocEnv)
    } else {
		assign("stockcol", cols, envir = .pRolocEnv)
	}
	invisible(prevcols)
}


library("dplyr")
surveys <- read.csv("http://datacarpentry.github.io/dc_zurich/data/portal_data_joined.csv")
head(surveys)

surveys %>%
  filter(weight < 5) %>%
  select(species_id, sex, weight)

surveys %>%
  mutate(weight_kg = weight / 1000) %>%
  filter(!is.na(weight)) %>%
  head

surveys %>%
  group_by(sex) %>%
  tally()

surveys %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE))

surveys %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE),
            min_weight = min(weight, na.rm = TRUE)) %>%
  filter(!is.nan(mean_weight))


quote(1:10)
quote(paste(letters, LETTERS, sep = "-"))


eval(quote(1 + 1))
eval(quote(1:10))

x <- 10
eval(quote(x + 1))

e <- new.env()
e$x <- 1
eval(quote(x + 1), env = e)

eval(quote(x), list(x = 30))

dfr <- data.frame(x = 1:10, y = LETTERS[1:10])
eval(quote(sum(x)), dfr)


x <- 10
substitute(sqrt(x))

e <- new.env()
e$x <- 1
substitute(sqrt(x), env = e)


parse(text = "1:10")
parse(file = "lineprof-example.R")


x <- 123
deparse(substitute(x))


foo <- "bar"
as.name(foo)
string <- "1:10"
parse(text=string)
eval(parse(text=string))


varName1 <- "varName2"
assign(varName1, "123")
varName1
get(varName1)
varName2


test <- function(x) {
    y <- deparse(substitute(x))
    print(y)
    print(x)
}
var <- c("one","two","three")
test(var)


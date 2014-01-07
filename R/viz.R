
## ------------------------------------------------------------------------
plot(1:10, 1:10)
rect(2, 2, 8, 8, col = "black")
rect(3, 3, 7, 7, col = "white")
abline(0, 1, col = "red")


## ------------------------------------------------------------------------
par(mfrow = c(2, 2))
plot(1:10, type = "p", main = "points (default)")
plot(1:10, type = "l", main = "lines")
plot(1:10, 10:1, type = "b", main = "both (points and lines)")
plot(1:10, type = "h", main = "histogram" )


## ------------------------------------------------------------------------
par(mfrow = c(1, 2))
barplot(c(1, 2, 3, 4))
s <- sample(letters, 1e3, replace = TRUE)
barplot(table(s))


## ------------------------------------------------------------------------
par(mfrow = c(1, 2))
l <- lapply(1:10, rnorm)
boxplot(l)
m <- matrix(rnorm(1e3), ncol = 10)
boxplot(m, names = LETTERS[1:10])


## ------------------------------------------------------------------------
par(mfrow = c(1, 2))
x <- rnorm(1e4)
hist(x)
hist(x, breaks = 50, freq = FALSE)
lines(density(x), col = "red")


## ------------------------------------------------------------------------
pie(c(1, 2, 3, 4))


## ------------------------------------------------------------------------
curve(x^2, 0, 10)


## ------------------------------------------------------------------------
m <- matrix(rnorm(30), ncol = 3)
dimnames(m) <- list(genes = paste("Gene", 1:10),
                    sample = LETTERS[1:3])
matplot(t(m), type = "b")


## ------------------------------------------------------------------------
heatmap(m, col = cm.colors(256))


## ------------------------------------------------------------------------
plot(1:10)
points(1:3, 3:1, pch = 19, col = "red")
lines(c(10, 1), c(1, 10))
lines(c(9, 1), c(1, 9), lty = "dotted", lwd = 3)
rect(8, 8, 10, 10, col = "black")
arrows(5, 1, 5, 10)
abline(v = 2, col = "blue")
abline(h = 2, col = "steelblue")
grid()


## ------------------------------------------------------------------------

(m <- matrix(c(1, 1, 2, 3, 4, 4, 4, 4), nrow = 2))
layout(m)
plot(1:10, main = 1)
curve(sin(x), from = -2, to = 2, main = 2)
curve(cos(x), from = -2, to = 2, main = 3)
image(t(m), main = 4)


## ----, eval=FALSE--------------------------------------------------------
## pdf("fig.pdf")
## plot(1:10)
## dev.off()


## ------------------------------------------------------------------------
library("scales")
cl <- col2hcl("steelblue", alpha = .1)
x <- rnorm(3e3)
y <- rnorm(3e3)
plot(x, y, pch = 19, col = cl, cex = 2)


## ------------------------------------------------------------------------
smoothScatter(x, y)


## ------------------------------------------------------------------------
brblramp <- colorRampPalette(c("brown", "steelblue"))
par(mfrow = c(2, 2))
plot(1:10, pch = 19, col = rainbow(10), cex = 3)
plot(1:10, pch = 19, col = cm.colors(10), cex = 3)
plot(1:10, pch = 19, col = heat.colors(10), cex = 3)
plot(1:10, pch = 19, col = brblramp(10), cex = 3)


## ------------------------------------------------------------------------
library("RColorBrewer")
display.brewer.all()


## ------------------------------------------------------------------------
library("camweather")
d <- nounits(weatherdata("2013-06-01"))
library("ggplot2")
p <- ggplot(data = d, aes(x = Time, y = Temp))
p + geom_point() ## add a layer of points 
p <- ggplot(data = d, aes(x = Time, y = Temp, colour = WindDr))
p + geom_point(size = 3) 
p <- ggplot(data = d, aes(x = Time, y = Temp, colour = WindDr, size = WindSp))
p + geom_point() 
p + geom_point() + facet_wrap(~WindDr)


## ------------------------------------------------------------------------
library("lattice")
xyplot(Temp ~ as.POSIXct(Time), data = d, col = d$WindDr, pch = 19)
xyplot(Temp ~ Press | WindDr, data = d, pch = 19)
splom(d[, c("Temp", "Press", "WindSp", "Humid")])



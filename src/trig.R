## Example script file, to generate a little plot.
## Demonstrate trigonometric functions.
## Sept 2007

x <- seq(from=0, to=2*pi, length=100)
y <- sin(x)
z <- cos(2*x)
z ## will not appear when source'd
print(y[1:10]) ## should use print()
plot(x, y, type='l')
lines(x, z, type='l', col='red')

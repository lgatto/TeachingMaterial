pdf(file='mfrow_eg.pdf', width=6,
    height=4)
par(mfrow=c(2,3))
par(mar=c(3.5, 3.5, 1.5, 0.5),
    mgp=c(2.5, 1, 0))
x <- seq(from=0, to=2*pi, len=100)
plot(x, sin(x), main="sin (x)",
     type='l')
plot(x, sin(2*x), main="sin (2x)",
     type='l')
plot(x, sin(3*x), main="sin (3x)",
     type='l')
plot(x, cos(x), main="cos (x)",
     type='l')
plot(x, cos(2*x), main="cos (2x)",
     type='l')
plot(x, cos(3*x), main="cos (3x)",
     type='l')
dev.off()

pdf("moredesigns.pdf", width = 10, height = 9)

library("ggplot2")

xs <- c(1:4, 6:9, 11:14, 16:19)
cex <- 6

plot(0, type = "n", xlim = range(xs), ylim = c(1,4),
     xlab = "Experimental units", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n")

points(xs, rep(4, 16), pch = 22, cex = cex, col = "grey")

points(xs, rep(3, 16), pch = 22, cex = cex, col = "black",
       bg = rep(c("white", "grey"), each = 8))

points(xs, rep(2, 16), pch = 22, cex = cex, col = "black",
       bg = c("white", "white", "grey", "grey",
              "white", "white", "grey", "white",
              "white", "white", "grey", "grey",
              "grey", "white", "grey", "grey"))

points(xs, rep(1, 16), pch = 22, cex = cex, col = "black",
       bg = c("grey", "white", "white", "grey",
              "white", "grey", "white", "grey",
              "white", "grey", "grey", "white",
              "grey", "white", "grey", "white"))

text(3, 3, "Biaised design", pos = 3, offset = 2)
text(3, 2, "Randomised design", pos = 3, offset = 2)
text(3, 1, "Block randomised design", pos = 3, offset = 2)
dev.off()




pdf("interactions.pdf", width = 10, height = 6)
strain <- rep(c("A", "B"), each = 10)
mns <- c(rep(5, 5),
         rep(7, 5),
         c(5, 6, 5, 6, 9),
         c(7, 8, 9, 10, 3))
food <- rep(c(rep("low", 5), rep("high", 5)), 2)
facet <- rep(1:5, 4)
x <- data.frame(strain, growth = mns, food, facet)
ggplot(aes(x = strain, y = growth, colour = food), data = x) +
    geom_point() +
    facet_grid(~facet)
dev.off()


pdf("power.pdf", width = 10, height = 7)
deltas <- c(2, 2, 2, 1, 2)
sds <- c(1, 1, 2, 1.5, 4)
sig <- c(0.01, 0.05, 0.05, 0.05, 0.05)
pw <- c(seq(0.5, 0.95, 0.05), 0.99)
res <- matrix(NA, 11, length(deltas))
rownames(res) <- pw
for (i in 1:length(deltas)) {
    res[, i] <- sapply(pw,
                       function(p) power.t.test(delta = deltas[i],
                                                sd = sds[i],
                                                sig.level = sig[i],
                                                power = p)$n)
}
matplot(res, type = "b",
        col = 1:5,
        lty = 1,
        xaxt = "n",
        xlab = "Power",
        ylab = "Sample size (per group)")
axis(side = 1, at = 1:length(pw), labels = pw)
legend("topleft",
       paste0("Delta: ", deltas, " sd: ", sds, " signif: ", sig),
       lty = 1, bty = "n",
       col = 1:6)
grid()
abline(h = 3)
dev.off()

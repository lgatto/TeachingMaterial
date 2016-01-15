pdf("moredesigns.pdf", width = 10, height = 9)

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

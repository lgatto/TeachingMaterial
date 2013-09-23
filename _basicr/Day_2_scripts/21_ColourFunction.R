colourRange <- function(pallette) {
        par(mar = c(1,1,1,1))
        plot(NULL, xlim = c(0,1), ylim = c(0,1), axes = F, xlab = "", ylab = "")
        rht = 1/length(pallette)
        for(i in 1:length(pallette)){
                rect(0, (i-1) * rht, 1, i*rht, col = pallette[i], border = NA)
        }
}

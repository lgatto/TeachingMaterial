fdrdiagramme <- function(ntests = 1000,
                         preal = 0.1,
                         power = 0.8,
                         sig.level = 0.05,
                         plot = TRUE) {

    text2 <- function(x, ...) {
        xoffset <- 0.6
        yoffset <- 0.1
        rect(x[1] - xoffset, x[2] - yoffset,
             x[1] + xoffset, x[2] + yoffset,
             col = "white")
        text(x[1], x[2], ...)

    }
    
    arrow2 <- function(x, y, text = "") {
        text((x[1] + y[1])/2, (x[2] + y[2])/2, text, pos = 2, offset = 0.5)
        segments(x[1], x[2], y[1], y[2])
    }


    eff <- preal * ntests
    noeff <- (1 - preal) * ntests
    tp <- eff * power
    fn <- eff * (1 - power)
    tn <- noeff * (1 - sig.level)
    fp <- noeff * sig.level
        fdr <- round(fp/(tp+fp), 3)
    
    if (plot) {
        plot(0, type = "n", xlab = "", ylab = "",
             xaxt = "n", yaxt = "n",
             bty = "n", xlim = c(-1, 5))
         
        p1 <- c(0, 0)
        p2 <- c(2, 0.5)
        p3 <- c(2, -0.5)
        p4 <- c(4, 0.75)
        p5 <- c(4, 0.25)
        p6 <- c(4, -0.25)
        p7 <- c(4, -0.75)

        arrow2(p1, p2, paste0("P(real)=", preal))
        arrow2(p1, p3)
        arrow2(p2, p4, paste0("power=", power))
        arrow2(p2, p5)
        arrow2(p3, p6, paste0("sig=", sig.level))
        arrow2(p3, p7)
        
        text2(p1, paste(ntests, "tests"))
        text2(p2, paste(eff, "real effect"))
        text2(p3, paste(noeff, "no effect"))

        text2(p4, paste(tp, "tests positive\nTP"))
        text2(p5, paste(fn, "tests negative\nFN"))
        
        text2(p6, paste(tn, "tests negative\nTN"))
        text2(p7, paste(fp, "tests positive\nFP"))

        title(main = paste0("False discovery rate: ",
                            fp, "/(", tp, "+", fp, ") = ", fdr))
    }

    invisible(fdr)
    
}

library("genefilter")
library("multtest")

mkdata <- function(n, nk = 15) {
  m <- 6
  x <- abs(rnorm(n) + rpois(n, 2))
  y <- replicate(m - 1, abs(x + rnorm(n)))
  ans <- cbind(x, y)
  colnames(ans) <- 1:m
  ## tweak data
  k <- sample(nrow(ans), nk)
  x <- rnorm(nk, 5)
  y <- rnorm(nk, 1)
  for (l in 1:3)
    ans[k, l] <- abs(x + rnorm(nk, 0, 0.2)) 
  for (l in 4:6)
    ans[k, l] <- abs(y + rnorm(nk, 0, 0.2))  
  return(ans)
}

mkname <- function(n) {
  dups <- TRUE
  while (dups) {
    nms <- replicate(n, paste0(paste(sample(LETTERS, 2), collapse = ""),
                               paste(sample(9, 8, replace = TRUE),
                                     collapse = "")))
    dups <- any(duplicated(nms))
  }
  return(sort(nms))
}

mkpw <- function(n, m = 100) 
  paste("PathWay", sample(m, n, replace = TRUE), sep = "_")



set.seed(1)


es <- mkdata(1000)

pd <- data.frame(sample = LETTERS[1:ncol(es)],
                 group = rep(c("ctrl", "cond"), each = 3),
                 replicate = rep(1:3, 2))

tt <- rowttests(es, pd$group)
adj <- mt.rawp2adjp(tt$p.value)
bh <- adj$adjp[order(adj$index), "BH"]

fd <- data.frame(genes = mkname(nrow(es)),
                 pathway = mkpw(nrow(es), 25),
                 pv = tt$p.value,
                 bh = bh)

ma1 <- list(expression = es,
            featuremeta = fd,
            samplemeta = pd)

write.csv(es, file = "MAdata1.csv")
write.csv(fd, file = "fmeta1.csv")
write.csv(pd, file = "smeta1.csv")
save(ma1, file = "MA1.rda")



makeall <- function(n = 1000, nk=15) {
  require("genefilter")
  require("multtest")
  es <- mkdata(n, nk)
  pd <- data.frame(sample = LETTERS[1:ncol(es)],
                   group = rep(c("ctrl", "cond"), each = 3),
                   replicate = rep(1:3, 2))

  tt <- rowttests(es, pd$group)
  adj <- mt.rawp2adjp(tt$p.value)
  bh <- adj$adjp[order(adj$index), "BH"]

  fd <- data.frame(genes = mkname(nrow(es)),
                   pathway = mkpw(nrow(es), 25),
                   pv = tt$p.value,
                   bh = bh)

  ma <- list(expression = es,
             featuremeta = fd,
             samplemeta = pd)
  return(ma)
}

ma2 <- makeall(nk = 5)
ma3 <- makeall(nk = 0)

es2 <- ma2$expression
es3 <- ma3$expression
fd2 <- ma2$featuremeta
fd3 <- ma3$featuremeta
sd2 <- ma2$samplemeta
sd3 <- ma3$samplemeta

write.csv(es2, file = "MAdata2.csv")
write.csv(fd2, file = "fmeta2.csv")
write.csv(sd2, file = "smeta2.csv")

write.csv(es3, file = "MAdata3.csv")
write.csv(fd3, file = "fmeta3.csv")
write.csv(sd3, file = "smeta3.csv")




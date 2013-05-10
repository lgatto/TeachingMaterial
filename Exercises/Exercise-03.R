expdata <- matrix(rnorm(200), nrow = 50, ncol = 4)

dimnames(expdata) <-
  list(features = paste0("gene", 1:nrow(expdata)),
       samples = paste0("sample", 1:ncol(expdata)))


smdata <- data.frame(feature = colnames(expdata),
                     group = c("ctrl", "ctrl",
                       "cond1", "cond1"),
                     replicate = rep(1:2, each = 2))

fmdata <- data.frame(feature = rownames(expdata),                     
                     description = paste("Important gene",
                       rownames(expdata)))

marray <- list(
  expression = expdata,
  featuremeta = fmdata,
  samplemeta = smdata))


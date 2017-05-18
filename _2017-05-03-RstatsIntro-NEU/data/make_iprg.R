iprg <- read.csv("iPRG_example_runsummary_0.csv")
iprg$TechReplicate <-
    factor(sub("\\.raw", "", sub(".+sample.+_", "", sub("-", "_", iprg$Run))))
save(iprg, file = "iprg.rda")
write.csv(iprg, file = "iPRG_example_runsummary.csv", row.names = FALSE)

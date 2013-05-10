con <- file("version1.dat", "r")
while (length(dat <- scan(con,n = 5,quiet = TRUE)) > 0) {
  print(mean(dat))
}
close(con)

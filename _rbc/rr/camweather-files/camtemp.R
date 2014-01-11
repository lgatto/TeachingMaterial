camtemp <- function(date) {
  ## Show weather for DATE given in YYYY-MMMM-DD form.

  matches <- grep(date, as.character(camweatherraw[,"Date"]))

  temps <- camweatherraw[matches,"Temp"] / 10
  times <- as.POSIXct(camweatherraw[matches,"Date"])

  plot(times, temps, type='l', main=date,
       ylim=c(-10, 30), ylab='Temperature (C)')
}

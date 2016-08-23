print.oCPD <-
function(x, ...) {
  cat("Changepoints:\n")
  if(is.null(x$time)) print(x$changes[-1])
  else                print(x$time[x$changes[-1]])
}

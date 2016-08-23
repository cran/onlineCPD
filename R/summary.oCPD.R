summary.oCPD <-
function(object, ...) { 
  cat("  An oCPD object:\n\n")
  
  print(object)
  
  cat("\n-",dim(object$data)[2], "- variate data\n\n")
  cat("- Colnames(data) is \'",colnames(object$data),"\'")
  
  if(is.null(object$time)) cat("\n\n- Time is NULL \n")
  else                     cat("\n\n- Time is not NULL \n")
}

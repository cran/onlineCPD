studentpdf <-
function(x, mu, var, nu) {
  c <- exp(lgamma(nu/2 + 0.5) - lgamma(nu/2)) / sqrt(nu * pi * var)
  
  return( c * (1 + (1/(nu * var)) * t(x - t(mu))^2)^(-(nu + 1)/2) )
}

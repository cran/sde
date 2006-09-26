checkBS <- function(theta){
  if(theta[2]<=0) stop("variance must be positive")
}

rcBS <- function(n=1, t, x0, theta){
  checkBS(theta)
  ml <- log(x0) + (theta[1]-theta[2]^2/2)*t
  sl <- sqrt(t)*theta[2]
  rlnorm(n, meanlog = ml, sdlog = sl)
}

dcBS <- function(x, t, x0, theta, log = FALSE){
  checkBS(theta)
  ml <- log(x0) + (theta[1]-theta[2]^2/2)*t
  sl <- sqrt(t)*theta[2]
  dlnorm(x, meanlog = ml, sdlog = sl, log=log)
}

pcBS <- function(x, t, x0, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkBS(theta)
  ml <- log(x0) + (theta[1]-theta[2]^2/2)*t
  sl <- sqrt(t)*theta[2]
  plnorm(x, meanlog = ml, sdlog = sl,
	lower.tail = lower.tail, log.p = log.p)
}

qcBS <- function(p, t, x0, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkBS(theta)
  ml <- log(x0) + (theta[1]-theta[2]^2/2)*t
  sl <- sqrt(t)*theta[2]
  qlnorm(p, meanlog = ml, sdlog = sl,
	lower.tail = lower.tail, log.p = log.p) 
}




rsBS <- function(n=1, theta){
  checkBS(theta)
  return(NA)
  rnorm(n, mean=theta[2], sd=theta[3]/(2*theta[1]))
}

dsBS <- function(x, theta, log = FALSE){
  checkBS(theta)
  return(NA)
  dnorm(x, mean=theta[2], sd=theta[3]/(2*theta[1]), log=log)
}

psBS <- function(x, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkBS(theta)
  return(NA)
  pnorm(x, mean=theta[2], sd=theta[3]/(2*theta[1]),
	lower.tail = lower.tail, log.p = log.p)
}

qsBS <- function(p, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkBS(theta)
  return(NA)
  qnorm(p, mean=theta[2], sd=theta[3]/(2*theta[1]),
	lower.tail = lower.tail, log.p = log.p) 
}

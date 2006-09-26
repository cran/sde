checkOU <- function(theta){
  if(theta[1]<=0) stop("the process is not stationary")
  if(theta[3]<=0) stop("variance must be positive")
}

rcOU <- function(n=1, t, x0, theta){
  checkOU(theta)
  rnorm(n, mean=theta[2]+(x0-theta[2])*exp(-theta[1]*t), 
    sd=theta[3]*sqrt((1-exp(-2*theta[1]*t))/(2*theta[1])))
}

dcOU <- function(x, t, x0, theta, log = FALSE){
  checkOU(theta)
  dnorm(x, mean=theta[2]+(x0-theta[2])*exp(-theta[1]*t), 
    sd=theta[3]*sqrt((1-exp(-2*theta[1]*t))/(2*theta[1])), log=log)
}

pcOU <- function(x, t, x0, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkOU(theta)
  pnorm(x, mean=theta[2]+(x0-theta[2])*exp(-theta[1]*t), 
    sd=theta[3]*sqrt((1-exp(-2*theta[1]*t))/(2*theta[1])),
	lower.tail = lower.tail, log.p = log.p)
}

qcOU <- function(p, t, x0, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkOU(theta)
  qnorm(p, mean=theta[2]+(x0-theta[2])*exp(-theta[1]*t), 
    sd=theta[3]*sqrt((1-exp(-2*theta[1]*t))/(2*theta[1])),
	lower.tail = lower.tail, log.p = log.p) 
}


rsOU <- function(n=1, theta){
  checkOU(theta)
  rnorm(n, mean=theta[2], sd=theta[3]/(2*theta[1]))
}

dsOU <- function(x, theta, log = FALSE){
  checkOU(theta)
  dnorm(x, mean=theta[2], sd=theta[3]/(2*theta[1]), log=log)
}

psOU <- function(x, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkOU(theta)
  pnorm(x, mean=theta[2], sd=theta[3]/(2*theta[1]),
	lower.tail = lower.tail, log.p = log.p)
}

qsOU <- function(p, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkOU(theta)
  qnorm(p, mean=theta[2], sd=theta[3]/(2*theta[1]),
	lower.tail = lower.tail, log.p = log.p) 
}

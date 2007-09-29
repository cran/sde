require(sde)
# ex3.10.R
# CIR-Model
# theta = -1
# alpha = 10
# sigma = 1
# x0 = 10
set.seed(123); 
d <- expression(10 - x)
s <- expression(sqrt(x)) 
x0 <- 10
sde.sim(X0=x0,drift=d, sigma=s,N=1000,delta=0.1) -> X

# estimator for alpha
(sum(X)^2)/(2*(length(X)*sum(X^2)-sum(X)^2))

# estimator for theta
(-length(X)*sum(X))/(2*(length(X)*sum(X^2)-sum(X)^2))

# ex1.09.R
set.seed(123)
n <- 100
Z <- rnorm(n)
h <- seq(1e-7, 1e-2,length=30)
W <- sum(Z*sapply(1:n, function(x) phi(x,0.5,T)))
for(i in h)
 Wh <- sum(Z*sapply(1:n, function(x) phi(x,0.5+i,T)))
inc.ratio <- abs(Wh-W)/h
plot(h,inc.ratio,type="l",log="y",
    ylab=expression(abs(B(0.5+h)-B(0.5))/h))
max(inc.ratio,na.rm=T)

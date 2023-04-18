source("R/opt.R")
set.seed(19630103)
p <- 1000
m <- p-1
num_groups <- 4
group_size <- p / num_groups
err.var <- 1
signal <- c(1,-1,2,-2)
beta0 <- c(rep(signal[1],group_size),rep(signal[2],group_size),rep(signal[3],group_size),rep(signal[4],group_size))

Y <- beta0 + sqrt(err.var)*rnorm(p)
tmp <- sort(rep(1:p,2))
idx <- tmp[1:(2*m)]
jdx <- tmp[2:(2*m+1)]
val <- rep(c(1,-1),m)

# diff <- 2*abs(Y-mean(Y))
# diff[1] <- diff[1] / 2; diff[p] <- diff[p] / 2
# lambda_max <- max(diff)
#lambda <- exp(seq(-log(1e3),log(1e3),length.out=10))
lambda <- c(1e3)

out <- filtering(Y,val,idx,jdx,lambda,m,varying_rho=FALSE)
############ complete subgraphs with random erroneous edges #############

source("R/opt.R")
cat("loading data....... \n")
set.seed(19630103)
n = 200
num.groups <- 200
num.vars.per.group <- 11
num.active.groups <- 4
cor <- .7
err.var <- 1


p <- num.groups * num.vars.per.group
beta.true <- rep(0,p)
X <- matrix(rnorm(n*p),n,p)
num.active.vars <- num.active.groups*num.vars.per.group
for (i in 1:num.groups)
{
	if (i <= num.active.groups)
	{
		beta.true[(1+(i-1)*num.vars.per.group):(i*num.vars.per.group)] <- floor((i+1)/2)*(-1)^{i+1}	
	}
	for (j in 2:num.vars.per.group)
	{
		X[,j+(i-1)*num.vars.per.group] = (sqrt(cor)*X[,1+(i-1)*num.vars.per.group]
									   + sqrt(1-cor)*X[,j+(i-1)*num.vars.per.group])
	}
}

X <- sweep(X,2,sqrt(n-1)*apply(X,2,sd),'/') # normalize columns of X
Y = X %*% beta.true + sqrt(err.var)*rnorm(n)

 ##### number of erroneous edges in the graph #####
# num.err.edges <- num.active.groups*num.vars.per.group*(num.vars.per.group-1)/2
num.err.edges <- (num.vars.per.group-1)*num.active.vars
m <- num.groups*num.vars.per.group*(num.vars.per.group-1)/2 + num.err.edges
jdx <- rep(0,2*m)

a <- sample(1:num.active.vars, replace=TRUE, num.err.edges)
b <- sample((1+num.active.vars):p,num.err.edges)
jdx[1:(2*num.err.edges)] <- as.vector(rbind(a,b))
c <- NULL
for (i in 1:(num.vars.per.group-1)){
	for (j in (i+1):num.vars.per.group){
		c <- c(c,c(i,j))
	}
}
for (i in 1:num.groups){
	jdx[(1+2*num.err.edges+(i-1)*num.vars.per.group*(num.vars.per.group-1)):(2*num.err.edges+i*num.vars.per.group*(num.vars.per.group-1))] <- (c + (i-1)*num.vars.per.group)
}
rm(a); rm(b); rm(c)
idx <- sapply(2:(2*m+1),function(x) floor(x/2))
val <- rep(c(1,-1),m)
lambda <- exp(seq(-log(1e10),log(1e10),length.out=100))


ptm <- proc.time()
out1 <- linreg_path_new(Y,X,val,idx,jdx,lambda,p,m)
time.diff <- proc.time() - ptm 
summary1 <- c(time.diff,mean(out1$admm_iter_num),sum(out1$chol_num))
save(summary1,file="variantADMM_100grid.RData")

ptm <- proc.time()
out2 <- linreg_path_new(Y,X,val,idx,jdx,lambda,p,m, variantADMM=FALSE)
time.diff <- proc.time() - ptm
summary2 <- c(time.diff,mean(out2$admm_iter_num),sum(out2$chol_num)) 
save(summary2,file="ADMM_100grid.RData")



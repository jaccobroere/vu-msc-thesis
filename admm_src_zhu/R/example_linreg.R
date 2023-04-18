source("R/opt.R")
cat("loading data....... \n")
set.seed(19630103)
n = 300
num.groups <- 20
num.vars.per.group <- 5
p <- num.groups * num.vars.per.group
num.active.groups <- 3
beta.true <- rep(0,p)
X <- matrix(rnorm(n*p),n,p)
cor <- .5
err.var <- .1
num.of.err.edges <- 2
num.of.active.vars <- num.active.groups*num.vars.per.group
for (i in 1:num.groups)
{
	if (i <= num.active.groups)
	{
		beta.true[(1+(i-1)*num.vars.per.group):(i*num.vars.per.group)] <- 1	
	}
	for (j in 2:num.vars.per.group)
	{
		X[,j+(i-1)*num.vars.per.group] = (sqrt(cor)*X[,1+(i-1)*num.vars.per.group]
									   + sqrt(1-cor)*X[,j+(i-1)*num.vars.per.group])
	}
}

X <- sweep(X,2,sqrt(n-1)*apply(X,2,sd),'/')
Y = X %*% beta.true + sqrt(err.var)*rnorm(n)
# beta_ols <- solve(t(X)%*%X,t(X)%*%Y)

I <- J <- rep(0,num.of.err.edges*num.active.groups+num.groups*(num.vars.per.group-1))
for (i in 1:num.groups)
{
	if (i <= num.active.groups)
	{
		I.tmp = sort(sample(((1+(i-1)*num.vars.per.group):(i*num.vars.per.group)),
		replace=TRUE, num.of.err.edges))
		J.tmp = sample((1+num.of.active.vars):p,num.of.err.edges)
		I.tmp <- c(rep(1+(i-1)*num.vars.per.group,num.vars.per.group-1),I.tmp)
		J.tmp <- c((2+(i-1)*num.vars.per.group):(i*num.vars.per.group),J.tmp)
	
		I[(1+(i-1)*(num.of.err.edges+num.vars.per.group-1)):
		(i*(num.of.err.edges+num.vars.per.group-1))] <- I.tmp
		J[(1+(i-1)*(num.of.err.edges+num.vars.per.group-1)):
		(i*(num.of.err.edges+num.vars.per.group-1))] <- J.tmp
	} else 
	{
		I.tmp <- rep(1+(i-1)*num.vars.per.group,num.vars.per.group-1)
		J.tmp <- (2+(i-1)*num.vars.per.group):(i*num.vars.per.group)
		I[(1+num.of.err.edges*num.active.groups+(i-1)*(num.vars.per.group-1)):
		(num.of.err.edges*num.active.groups+i*(num.vars.per.group-1))] <- I.tmp
		J[(1+num.of.err.edges*num.active.groups+(i-1)*(num.vars.per.group-1)):
		(num.of.err.edges*num.active.groups+i*(num.vars.per.group-1))] <- J.tmp
	}
}

m <- length(I)
jdx <- rep(0,2*m)
jdx[seq(1,(2*m-1),2)] <- I
jdx[seq(2,2*m,2)] <- J
rdx <- seq(1,2*m+1,2)
val <- rep(c(1,-1),m)
lambda <- exp(seq(-log(1e2),log(1e2),length.out=50))*log(p)/n
beta_path <- linreg_path(Y,X,val,jdx,rdx,lambda,p,m)



gen_data <- function(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var){
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
	X <- sweep(X,2,apply(X,2,mean),'-')
	X <- sweep(X,2,sqrt(n-1)*apply(X,2,sd),'/') # normalize columns of X
	Y = X %*% beta.true + sqrt(err.var)*rnorm(n)
	
	 ##### number of erroneous edges in the graph #####
	# num.err.edges <- num.active.groups*num.vars.per.group*(num.vars.per.group-1)/2
	num.err.edges <- (num.vars.per.group-1)*num.active.vars
	m <- num.groups*num.vars.per.group*(num.vars.per.group-1)/2 + num.err.edges
	jdx <- rep(0,2*m)

	# a <- sample(1:num.active.vars, replace=TRUE, num.err.edges)
	# b <- sample((1+num.active.vars):p,num.err.edges)
	# jdx[1:(2*num.err.edges)] <- as.vector(rbind(a,b))
	
	for (i in 1:num.active.vars) 
	{
		a <- rep(i,num.vars.per.group-1)
		b <- sample((1+num.active.vars):p,num.vars.per.group-1)
		jdx[(1+2*(i-1)*(num.vars.per.group-1)):(2*i*(num.vars.per.group-1))] <- as.vector(rbind(a,b))
	}
	
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
	return (list(Y=Y,X=X,val=val,idx=idx,jdx=jdx,p=p,m=m,beta.true=beta.true))
}
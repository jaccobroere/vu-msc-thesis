require(genlasso)
source("R/opt.R")
source("R/gen_data.R")

########### compare four methods on 20 different values of lambda #########
timing_comp_single_lam <- function(num.groups, gamma)
{
  cor = .7 
  num.of.grid = 20
  minlam = 1e-4
  n = 200
  num.vars.per.group <- 11
  num.active.groups <- 4
  err.var <- .1
  p <- num.vars.per.group*num.groups
  set.seed(19630103)
  data <- gen_data(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var)
  D = sparseMatrix(i=data$idx,j=data$jdx,x=data$val)
  out.tmp = fusedlasso(data$Y,data$X,D=D,minlam=minlam,gamma = gamma)
  lambda_graph <- exp(seq(log(max(out.tmp$lambda)),log(minlam),length.out=num.of.grid))
  dir <- "../code_fgfl_aaai14/data/"
  write.table(lambda_graph,file=paste(dir,"lambdas_single_p",p,"_gamma_",gamma,".txt",sep=""),row.name=FALSE,col.names=FALSE)
  obj_val <- function(beta,lambda)
  {
    tmp = norm(data$Y - data$X %*% beta,"F")
    fun.val = .5*tmp*tmp + lambda*(norm(D %*% beta,"o")+gamma*norm(matrix(beta),"o"))
    fun.val
  }
  
  ########## get approximate true optimal solution ###########
  out.tmp <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,gamma=gamma,data[[6]],data[[7]],eps_abs=1e-14,eps_rel=1e-13)
  
  
  ########### set termination tolerance so that sub-optimalities are similar #############
  eps_augADMM = 1e-6
  eps_stanADMM = 1e-5
  
  fun.val.genlasso <- rep(0,length(lambda_graph))
  fun.val.augADMM <- rep(0,length(lambda_graph))
  fun.val.stanADMM <- rep(0,length(lambda_graph))
  fun.val.true <- rep(0,length(lambda_graph))
  
  time.genlasso <- rep(0,length(lambda_graph))
  time.augADMM  <- rep(0,length(lambda_graph))
  time.stanADMM <- rep(0,length(lambda_graph))
  
  for (i in seq_along(lambda_graph))
  {
    ########## run three methods #############
    ptm <- proc.time()
    out.genlasso = fusedlasso(data$Y,data$X,D=D,minlam=lambda_graph[i], gamma = gamma, maxsteps=1e6)
    time.diff.genlasso <- proc.time() - ptm
    
    ptm <- proc.time()
    out.augADMM <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph[i],gamma=gamma,data[[6]],data[[7]],eps_abs=.1*eps_augADMM,eps_rel=eps_augADMM)
    time.diff.augADMM <- proc.time() - ptm
    
    ptm <- proc.time()
    out.stanADMM <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph[i],gamma=gamma,data[[6]],data[[7]],eps_abs=.1*eps_stanADMM,eps_rel=eps_stanADMM,standard_ADMM=TRUE)
    time.diff.stanADMM <- proc.time() - ptm
    
    ################### store timing #############
    time.genlasso[i] = time.diff.genlasso[3]
    time.augADMM[i] = time.diff.augADMM[3]
    time.stanADMM[i] = time.diff.stanADMM[3]
    
    ########### store objective function ###########
    beta.genlasso = coef(out.genlasso, lambda_graph[i])$beta
    fun.val.genlasso[i] <- obj_val(beta.genlasso,lambda_graph[i])
    fun.val.augADMM[i] <- obj_val(matrix(out.augADMM$beta_path[,1]),lambda_graph[i])
    fun.val.stanADMM[i] <- obj_val(matrix(out.stanADMM$beta_path[,1]),lambda_graph[i])
    fun.val.true[i] <- min(obj_val(matrix(out.tmp$beta_path[,i]),lambda_graph[i]),fun.val.genlasso[i],fun.val.augADMM[i],fun.val.stanADMM[i])
  }
  
  funVals <- list(fun.val.true=fun.val.true,fun.val.genlasso=fun.val.genlasso,fun.val.augADMM=fun.val.augADMM,fun.val.stanADMM=fun.val.stanADMM)
  save(funVals,file=paste("tmp/funVals_single_p",p,"_gamma_",gamma,".RData",sep=""))
  runtime = list(time.genlasso=time.genlasso,time.augADMM=time.augADMM,time.stanADMM=time.stanADMM)
  save(runtime,file=paste("tmp/timing_single",p,"_gamma_",gamma,".RData",sep="")) 
}


args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)!=2){
  stop("please supply both num.groups and gamma.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

timing_comp_single_lam(num.groups,gamma)


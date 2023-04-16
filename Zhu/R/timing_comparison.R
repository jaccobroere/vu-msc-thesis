require(genlasso)
source("R/opt.R")
source("R/gen_data.R")

timing_comp_path = function(num.groups,gamma)
{
  n <- 200
  num.vars.per.group <- 11
  num.active.groups <- 4
  err.var <- .1
  cor <- .7
  p <- num.groups * num.vars.per.group
  num.of.grid <- 1e2
  set.seed(19630103)
  data <- gen_data(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var)
  D = sparseMatrix(i=data$idx,j=data$jdx,x=data$val)
  num.of.grid <- 1e2
  lambda_graph <- exp(seq(log(1e4),log(1e-4),length.out=num.of.grid))
  obj_val <- function(beta)
  {
    fun.val <- rep(0,dim(beta)[2])
    for (i in 1:dim(beta)[2]) {
      tmp = norm(data$Y - data$X %*% beta[,i],"F")
      fun.val[i] = .5*tmp*tmp + lambda_graph[i]*(norm(D %*% beta[,i],"o")+gamma*norm(matrix(beta[,i]),"o"))
    }
    fun.val
  }
  #degree <- apply(D,2,function(x) sum(abs(x)))
  
  ###############  get estimate of best lambda_graph ##########
  # lambda_graph <-  exp(seq(log(1e3),-log(1e3),length.out=2e2))
  # out <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],graph_only=TRUE)
  # opt.score = rep(0,length(lambda_graph))
  # for (i in 1:length(lambda_graph)){
  #   opt.score[i] <- norm(matrix(data$beta.true - out$beta_path[,,i]), "F")
  # }
  # lambda_star <- lambda_graph[which.min(opt.score)]
  # lambda_graph <- exp(seq(log(1e3), log(lambda_star/2),length.out=num.of.grid))
  
  
  ########### set termination tolerance so that suboptimalities are similar #############
  
  eps_augADMM = 1e-6
  eps_stanADMM = 1e-5
  
  if (num.groups == 10){
    eps_augADMM = 1e-8
    eps_stanADMM = 1e-5
  }
  
  if (num.groups == 20) # checked
  {
    eps_augADMM = 1e-5
    eps_stanADMM = 1e-4
  }
  
  if (num.groups == 40) # checked
  {
    eps_augADMM = 1e-6
    eps_stanADMM = 1e-4
  }
  
  if (num.groups == 60) # checked 
  { 
    eps_augADMM = 1e-5
    eps_stanADMM = 1e-4
  }
  
  if (num.groups == 80){
    eps_augADMM = 1e-6
    eps_stanADMM = 1e-4
  }
  
  
  if (num.groups == 100){
    eps_augADMM = 1e-6
    eps_stanADMM = 1e-4
  }
  
  if (num.groups == 120){
    eps_augADMM = 1e-6
    eps_stanADMM = 1e-5
  }
  
  
  if (num.groups == 140){
    eps_augADMM = 1e-6
    eps_stanADMM = 1e-5
  }
  
  
  if (num.groups == 160){
    eps_augADMM = 1e-6
    eps_stanADMM = 1e-5
  }
  
  if (num.groups == 180){
    eps_augADMM = 1e-6
    eps_stanADMM = 1e-5
  }
  
  if (num.groups == 200){
    eps_augADMM = 1e-6
    eps_stanADMM = 1e-5
  }
  
  if (num.groups == 260){
    eps_augADMM = 1e-6
    eps_stanADMM = 1e-5
  }
  
  if (num.groups == 300){
    eps_augADMM = 1e-6
    eps_stanADMM = 1e-5
  }
  
  
  
  
  
  plot_rel_diff <- function(val1,val2,caption="suboptimality"){
    rel_diff <- (val1 - val2) / (abs(val2) + 1e-10)
    plot(rel_diff,xlab="lambda",ylab="relative differences",main=caption)
  }
  
  ############ run augADMM for 100 lambdas #############
  ptm <- proc.time()
  #out <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],eps_abs=.1*eps_augADMM,eps_rel=eps_augADMM,graph_only=TRUE)
  out.augADMM100 <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,gamma=gamma,data[[6]],data[[7]],eps_abs=.1*eps_augADMM,eps_rel=eps_augADMM)
  time.diff.augADMM <- proc.time() - ptm
  fun.val.augADMM <- obj_val(out.augADMM100$beta_path)
  cat("time running augADMM:",time.diff.augADMM[3],"\n",sep="")
  beta=out.augADMM100$beta_path
  save(beta,file=paste("tmp/beta_p",p,".RData",sep=""))
  
  
  ########## run standard ADMM for 100 lambdas ###############
  ptm <- proc.time()
  #out <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],eps_abs=.1*eps_stanADMM,eps_rel=eps_stanADMM,graph_only=TRUE,standard_ADMM=TRUE)
  out.stanADMM100 <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,gamma=gamma,data[[6]],data[[7]],eps_abs=.1*eps_stanADMM,eps_rel=eps_stanADMM,standard_ADMM=TRUE)
  time.diff.stanADMM <- proc.time() - ptm
  fun.val.stanADMM <- obj_val(out.stanADMM100$beta_path)
  cat("time running stanADMM:",time.diff.stanADMM[3],"\n",sep="")
  
  ########## run genlasso #############
  ptm <- proc.time()
  out.genlasso = fusedlasso(data$Y,data$X,D=D,gamma=gamma, minlam=min(lambda_graph),maxsteps=1e6)
  time.diff.genlasso <- proc.time() - ptm
  cat("time running genlasso:",time.diff.genlasso[3],"\n",sep="")
  beta.genlasso = coef(out.genlasso, lambda_graph)$beta
  fun.val.genlasso <- obj_val(beta.genlasso)
  
  ######## save lambdas #########
  dir <- "../code_fgfl_aaai14/data/"
  write.table(lambda_graph,file=paste(dir,"lambdas100_p",p,".txt",sep=""),row.name=FALSE,col.names=FALSE)
  
  
  ####### plot relative differences of augADMM and stanADMM #########
  pdf(paste("tmp/sub_optimality_p",p,"_gamma_",gamma,".pdf",sep=""))
  par(mfrow=c(3,1))
  plot_rel_diff(fun.val.augADMM,fun.val.genlasso,caption="augADMM vs genlasso")
  plot_rel_diff(fun.val.stanADMM,fun.val.genlasso,caption="stanADMM vs genlasso")
  plot_rel_diff(fun.val.augADMM,fun.val.stanADMM,caption="augADMM vs stanADMM")
  dev.off()
  
  ####### report timing and suboptimality ###########
  cat("genlasso:", time.diff.genlasso[3],"\n","augADMM100:",time.diff.augADMM[3],"\n","stanADMM100:",time.diff.stanADMM[3],"\n",sep="")
  
  
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

timing_comp_path(num.groups,gamma)


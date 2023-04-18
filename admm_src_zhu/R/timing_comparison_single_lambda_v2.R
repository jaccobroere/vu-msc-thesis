require(genlasso)
source("R/opt.R")
source("R/gen_data.R")
require(Matrix)

############ compare three methods at 20 different values of lambda with ################
############  precision equal to that of genlasso #####################
timing_comp_single_lam <- function(num.groups, gamma, subopt=c(1e-2,1e-4,1e-6,1e-8))
{
  EPS = min(subopt)
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
  
  
  ######### run genlasso just to get the lambda_max and construct lambda_graph #########
  D = sparseMatrix(i=data$idx,j=data$jdx,x=data$val)
  out.tmp = fusedlasso(data$Y,data$X,D=D,minlam=minlam,gamma = gamma)
  lambda_graph <- exp(seq(log(max(out.tmp$lambda)),log(minlam),length.out=num.of.grid))
  dir <- "../code_fgfl_aaai14/data/"
  write.table(lambda_graph,file=paste(dir,"lambdas_single_p",p,"_gamma_",gamma,".txt",sep=""),row.name=FALSE,col.names=FALSE)
  write.table(subopt,file=paste(dir,"subopt",".txt",sep=""),row.name=FALSE,col.names=FALSE)
  
  obj_val <- function(beta,lambda)
  {
    tmp = norm(data$Y - data$X %*% beta,"F")
    fun.val = .5*tmp*tmp + lambda*(norm(D %*% beta,"o")+gamma*norm(matrix(beta),"o"))
    fun.val
  }
  rel_diff <- function(val1,val2){ 
    tmp <- (val1 - val2) / (abs(val2) + 1e-10) 
    tmp[tmp < .0] = .0
    tmp
  }
  
  
  max_iter_augADMM <- rep(2e3,length(lambda_graph))
  max_iter_stanADMM <- rep(2e3,length(lambda_graph))
  fun.val.star <- rep(0,length(lambda_graph))
  
  ########## get iteration numbers that achieves EPS suboptimality ###########
  iter.augADMM = iter.stanADMM = matrix(0,length(lambda_graph),length(subopt))
  for (i in seq_along(lambda_graph)){
    cat("lambda is: ", lambda_graph[i],"\n")
    ###### get optimal value at lambda_i #######
    out.tmp <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph[i],gamma=gamma,data[[6]],data[[7]],eps_abs=1e-16,eps_rel=1e-15)
    fun.val.star[i] <- out.tmp$fun_val
    repeat
    {
      out.report <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph[i],gamma=gamma,data[[6]],data[[7]],reporting=TRUE,max_num_iter=max_iter_augADMM[i])
      funPath.augADMM = rel_diff(out.report$fun_path,fun.val.star[i])
      if (funPath.augADMM[max_iter_augADMM[i]] < EPS) {
        for (k in seq_along(subopt)) iter.augADMM[i,k] = min(which(funPath.augADMM < subopt[k]))
        break
      }
      max_iter_augADMM[i] = 2*max_iter_augADMM[i]
      cat("maxiter augADMM: ",max_iter_augADMM[i],"\t")
    }
    
    repeat
    {
      out.report <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph[i],gamma=gamma,data[[6]],data[[7]],reporting=TRUE,max_num_iter=max_iter_stanADMM[i],standard_ADMM=TRUE)
      funPath.stanADMM = rel_diff(out.report$fun_path,fun.val.star[i])
      if (funPath.stanADMM[max_iter_stanADMM[i]] < EPS) {
        for (k in seq_along(subopt)) iter.stanADMM[i,k] = min(which(funPath.stanADMM < subopt[k]))
        break
      }
      max_iter_stanADMM[i] = 2*max_iter_stanADMM[i]
      cat("maxiter stanADMM: ",max_iter_stanADMM[i],"\t")
    }
    cat("\n")
  }
  
  ########### report timing for different lambdas' and different subopts #########
  time.augADMM  <- matrix(0,length(lambda_graph),length(subopt))
  time.stanADMM <- matrix(0,length(lambda_graph),length(subopt))
  
  for (i in seq_along(lambda_graph))
  {
    for (k in seq_along(subopt)){
      ########## run three methods #############
      ptm <- proc.time()
      out.augADMM <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph[i],gamma=gamma,data[[6]],data[[7]],max_num_iter=iter.augADMM[i,k])
      t <- proc.time() - ptm; time.augADMM[i,k] = t[3]
      
      ptm <- proc.time()
      out.stanADMM <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph[i],gamma=gamma,data[[6]],data[[7]],standard_ADMM=TRUE,max_num_iter=iter.stanADMM[i,k])
      t <- proc.time() - ptm; time.stanADMM[i,k] = t[3]
    }
  }
  
  ######### save timing and optimal values##########
  runtime = list(augADMM=time.augADMM,stanADMM=time.stanADMM)
  save(runtime,file=paste("tmp/timing_20lambda_",p,"_gamma_",gamma,".RData",sep=""))
  write.table(fun.val.star,file=paste(dir,"opt_funVals_20lambdas_p",p,"_gamma_",gamma,".txt",sep=""),row.name=FALSE,col.names=FALSE)
  iteration = list(augADMM=iter.augADMM,stanADMM=iter.stanADMM)
  save(iteration,file=paste("tmp/iteration_20lambda_",p,"_gamma_",gamma,".RData",sep=""))
}

######### then run fGFL with lambda_graph and optimal function values as input ###########
######### to get its runtimes at each lambda and each specified suboptimality. ##########


args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)<2){
  stop("please supply num.groups and gamma")
}else if (length(args)==2){
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
  timing_comp_single_lam(num.groups,gamma)
} else if (length(args)==3){
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
  timing_comp_single_lam(num.groups,gamma,subopt)
} else {
  stop("too many input arguments.")
}


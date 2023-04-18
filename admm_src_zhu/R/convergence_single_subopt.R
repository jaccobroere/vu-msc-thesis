source("R/opt.R")
source("R/gen_data.R")
require(Matrix)
require(irlba)
check_convergence <- function(num.groups,lambda_graph,gamma, max_iter=1e3){
  rel_diff <- function(val1,val2){ 
    tmp <- (val1 - val2) / (abs(val2) + 1e-10) 
    tmp[tmp < .0] = .0
    tmp
  }
  
  lambda_sparse = gamma * lambda_graph
  n <- 200
  num.vars.per.group <- 11
  num.active.groups <- 4
  err.var <- 20
  cor <- .7
  p <- num.groups*num.vars.per.group
  set.seed(19630103)
  data <- gen_data(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var)
  A = sparseMatrix(i=data$idx,j=data$jdx,x=data$val)
  #max_sig = max(2 * apply(A,2,function(x) sum(abs(x))) + 1)
  max_sig = (irlba(A,nu=1,nv=1)$d)^2 + 1
  ############## get fGFL reslults #############
  fun_val_fgfl <- as.matrix(read.table(paste("tmp/funVal_p",p,"_lambda_graph_",lambda_graph,"_gamma_",gamma,".txt",sep="")))
  
  ############ get approximate optimal value ##############
  out.tmp <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,gamma,data[[6]],data[[7]],standard_ADMM=TRUE,reporting=TRUE,max_num_iter=1e4)
  fun.val.star = min(c(out.tmp$fun_path,fun_val_fgfl))
  
  funPath.fGFL = rel_diff(fun_val_fgfl,fun.val.star)
  if (min(funPath.fGFL) > 1e-10) warning("fGFL does not converge within 1e-10 of optimal value. please tighten the stopping criterion of fGFL! \n")
  
  funPath.augADMM = funPath.stanADMM = 1
  repeat
  {
    cat("max_iter: ", max_iter, ".\n")
    ############ running Augmented ADMM #############
    out.report <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,gamma,data[[6]],data[[7]],reporting=TRUE,max_num_iter=max_iter)
    out.report2 <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,gamma,data[[6]],data[[7]],reporting=TRUE,max_num_iter=max_iter,diag_mat=max_sig)
    out.report3 <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,gamma,data[[6]],data[[7]],reporting=TRUE,max_num_iter=max_iter,diag_mat=5*max_sig)
    out.report4 <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,gamma,data[[6]],data[[7]],reporting=TRUE,max_num_iter=max_iter,diag_mat=10*max_sig)
    ############# running standard ADMM ##############
    out.report.stand <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,gamma,data[[6]],data[[7]],standard_ADMM=TRUE,reporting=TRUE,max_num_iter=max_iter)
    
    ########## calculate relative differences ########### 
    funPath.augADMM = rel_diff(out.report$fun_path,fun.val.star)
    funPath.augADMM2 = rel_diff(out.report2$fun_path,fun.val.star)
    funPath.augADMM3 = rel_diff(out.report3$fun_path,fun.val.star)
    funPath.augADMM4 = rel_diff(out.report4$fun_path,fun.val.star)
    funPath.stanADMM = rel_diff(out.report.stand$fun_path,fun.val.star)
    
    if (max(min(funPath.augADMM),min(funPath.stanADMM),min(funPath.augADMM2),min(funPath.augADMM3),min(funPath.augADMM4)) < 1e-10) break
    max_iter = 2*max_iter
  }
  
  
  ############ save timing and relative differences ##########
  out <- list(funPath.augADMM=funPath.augADMM,funPath.augADMM2=funPath.augADMM2,funPath.augADMM3=funPath.augADMM3,funPath.augADMM4=funPath.augADMM4,funPath.stanADMM=funPath.stanADMM,funPath.fGFL=funPath.fGFL)
  output_file = paste("tmp/convergence_funPaths_p",p,"_lam1_",lambda_graph,"_lam2_",lambda_sparse,".RData",sep="")
  save(out,file=output_file)
}

args=(commandArgs(TRUE))

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}
if(length(args) == 3){
  check_convergence(num.groups,lambda_graph,gamma)
} else if (length(args == 4)){
  check_convergence(num.groups,lambda_graph,gamma, max_iter=max_iter)
} else {
  stop("please supply num.groups, lambda_graph, gamma, and max_iter!")
}

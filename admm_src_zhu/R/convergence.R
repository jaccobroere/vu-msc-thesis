source("R/opt.R")
source("R/gen_data.R")

check_convergence <- function(num.groups,lambda_graph,lambda_sparse, max_iter=1e3){
  rel_diff <- function(val1,val2){ 
    tmp <- (val1 - val2) / (abs(val2) + 1e-10) 
    tmp[tmp < .0] = .0
    tmp
  }
  
  n <- 200
  num.vars.per.group <- 11
  num.active.groups <- 4
  err.var <- .1
  cor <- .7
  p <- num.groups*num.vars.per.group
  set.seed(19630103)
  data <- gen_data(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var)
  
  ############## get fGFL reslults #############
  fun_val_fgfl <- as.matrix(read.table(paste("tmp/funVal_p",p,".txt",sep="")))
  #time.fGFL = as.numeric(read.table(paste("tmp/convergence_timing_p",p,".txt",sep="")))
  
  ############ get approximate optimal value ##############
  out.tmp <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],standard_ADMM=TRUE,reporting=TRUE,max_num_iter=1e4)
  fun.val.star = min(c(out.tmp$fun_path,fun_val_fgfl))
  
  funPath.fGFL = rel_diff(fun_val_fgfl,fun.val.star)
  if (min(funPath.fGFL) > 1e-10) warning("fGFL does not converge within 1e-10 of optimal value. please tighten the stopping criterion of fGFL! \n")
  
  funPath.augADMM = funPath.stanADMM = 1
  repeat
  {
    cat("max_iter: ", max_iter, ".\n")
    ############ running Augmented ADMM #############
    out.report <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],reporting=TRUE,max_num_iter=max_iter)
    
    ############# running standard ADMM ##############
    out.report.stand <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],standard_ADMM=TRUE,reporting=TRUE,max_num_iter=max_iter)
    
    
    ########## calculate relative differences ########### 
    funPath.augADMM = rel_diff(out.report$fun_path,fun.val.star)
    funPath.stanADMM = rel_diff(out.report.stand$fun_path,fun.val.star)
    
    if (max(min(funPath.augADMM),min(funPath.stanADMM)) < 1e-10) break
    max_iter = 2*max_iter
  }
  
  ########## get runtime at which the sub-optimality reaches 1e-4,1e-6,1e-8,1e-10. 
  precision = c(1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10)
  iter.augADMM = iter.stanADMM = iter.fGFL = time.augADMM = time.stanADMM = rep(0,length(precision))
  for (k in seq_along(precision)){
    eps = precision[k]
    iter.augADMM[k] = min(which(funPath.augADMM < eps))
    iter.stanADMM[k] = min(which(funPath.stanADMM < eps))
    iter.fGFL[k] = min(which(funPath.fGFL < eps))
    
    ptm <- proc.time()
    linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],max_num_iter=iter.augADMM[k])
    t = proc.time() - ptm; time.augADMM[k] = t[3];
    
    ptm <- proc.time()
    out.report.stand <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],standard_ADMM=TRUE,max_num_iter=iter.stanADMM[k])
    t = proc.time() - ptm; time.stanADMM[k] <- t[3];
  }
  cat("iter augADMM: ", iter.augADMM,"\n")
  cat("iter stanADMM: ", iter.stanADMM,"\n")
  cat("iter fGFL: ", iter.fGFL,"\n")
  
  ########### save runtime at four different termination error ############
  output_file = paste("tmp/convergence_timing_seq_p",p,"_lam1_",lambda_graph,"_lam2_",lambda_sparse,".RData",sep="")
  tmp.time = list(time.stanADMM=time.stanADMM,time.augADMM=time.augADMM)
  save(tmp.time,file=output_file)
  
  
  ############ save timing and relative differences ##########
  out <- list(funPath.augADMM=funPath.augADMM,funPath.stanADMM=funPath.stanADMM,funPath.fGFL=funPath.fGFL)
  output_file = paste("tmp/convergence_funPaths_p",p,"_lam1_",lambda_graph,"_lam2_",lambda_sparse,".RData",sep="")
  save(out,file=output_file)
  #output_file = paste("tmp/convergence_timing_p",p,"_lam1_",lambda_graph,"_lam2_",lambda_sparse,"_iter_",max_iter,".RData",sep="")
  #out.timing <- list(time.augADMM=time.augADMM,time.stanADMM=time.stanADMM,time.fGFL=time.fGFL)
  #save(out.timing,file=output_file)
  
  
  ############ save fGFL iterations #############
  output_file = paste("../code_fgfl_aaai14/data/fGFL_iter_seq_p",p,"_lam1_",lambda_graph,"_lam2_",lambda_sparse,".txt",sep="")
  write.table(iter.fGFL,file=output_file,row.name=FALSE,col.names=FALSE)
}

args=(commandArgs(TRUE))

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}
if(length(args) == 3){
  check_convergence(num.groups,lambda_graph,lambda_sparse)
} else if (length(args == 4)){
  check_convergence(num.groups,lambda_graph,lambda_sparse, max_iter=max_iter)
} else {
  stop("please supply num.groups, lambda_graph, lambda_sparse, and max_iter!")
}

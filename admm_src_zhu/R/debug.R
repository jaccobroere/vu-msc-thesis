source("R/opt.R")
source("R/gen_data.R")
set.seed(19630103)
n <- 200
num.vars.per.group <- 11
num.active.groups <- 4
err.var <- .1
num.groups <- 200
cor <- .7
p <- num.groups * num.vars.per.group
data <- gen_data(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var)

lambda_graph_single <- c(0.021544)
lambda_sparse_single <- c(0.001000)

out_single_long <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
								lambda_graph_single,lambda_sparse_single,
								data[[6]],data[[7]],varying_rho=FALSE,max_num_iter=2e4,
								eps_abs=.0, eps_rel=.0)

opt_val <- out_single_long$fun_val[1]


out.acc.specialADMM <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
								lambda_graph_single,lambda_sparse_single,
								data[[6]],data[[7]],max_num_iter=2e3)
								
								#### choose another set of lambda's 
								lambda_graph_single <- 0.001000
								lambda_sparse_single <- 0.464159

								out_single_long <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
																lambda_graph_single,lambda_sparse_single,
																data[[6]],data[[7]],max_num_iter=2e4,
																eps_abs=.0, eps_rel=.0)

								opt_val <- out_single_long$fun_val[1]


								out.acc.specialADMM <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
																lambda_graph_single,lambda_sparse_single,
																data[[6]],data[[7]],max_num_iter=2e3,reporting=TRUE)
								fun.path.acc.specialADMM <- log10(out.acc.specialADMM$fun_path - opt_val)
								save(fun.path.acc.specialADMM,file="fun.acc.specialADMM.RData")


								out.acc.ADMM <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
																lambda_graph_single,lambda_sparse_single,
																data[[6]],data[[7]],max_num_iter=2e3,variantADMM=FALSE,reporting=TRUE)
								fun.path.acc.ADMM <- log10(out.acc.ADMM$fun_path  - opt_val)
								save(fun.path.acc.ADMM,file="fun.acc.ADMM.RData")

								out.nonacc.specialADMM <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
																lambda_graph_single,lambda_sparse_single,
																data[[6]],data[[7]],max_num_iter=2e3,varying_rho=FALSE,reporting=TRUE)
								fun.path.nonacc.specialADMM <- log10(out.nonacc.specialADMM$fun_path  - opt_val)
								save(fun.path.nonacc.specialADMM,file="fun.nonacc.specialADMM.RData")


								out.nonacc.ADMM <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
																lambda_graph_single,lambda_sparse_single,
																data[[6]],data[[7]],max_num_iter=2e3,varying_rho=FALSE,variantADMM=FALSE,reporting=TRUE)
								fun.path.nonacc.ADMM <- log10(out.nonacc.ADMM$fun_path - opt_val)
								save(fun.path.nonacc.ADMM,file="fun.nonacc.ADMM.RData")


								plot(fun.path.acc.specialADMM,type="l",col="blue",xlab="ADMM iterations",ylab="function value")
								lines(fun.path.acc.ADMM,type="l",col="red")
								lines(fun.path.nonacc.specialADMM,type="l",col="blue",lty=2)
								lines(fun.path.nonacc.ADMM,type="l",col="red",lty=2)
								legend(800,1.5, c("accelerated special ADMM","accelerated standard ADMM","special ADMM","standard ADMM"), 
										col = c("blue","red","blue","red"), lty=c(1,1,2,2),bty="n",cex=.6)
		
		
								out <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
										                    lambda_graph[1],lambda_sparse[1],
										                    data[[6]],data[[7]])
######### delete from timing_comparison.R
#lambda_graph <- exp(seq(log(max(a2$lambda)),log(min(a2$lambda)),length.out=num.of.grid))
#lambda_graph <- seq(min(a2$lambda),max(a2$lambda),length.out=1e3)
#lambda_graph <- a2$lambda
#save(beta.out, file="beta.genlasso.RData")



######### run augADMM again for 1000 lambdas ###########
num.of.grid <- 1e3
#lambda_graph <- exp(seq(log(1e3), log(lambda_star/2),length.out=num.of.grid))
#lambda_graph <- seq(max(a2$lambda),min(a2$lambda),length.out=num.of.grid)
lambda_graph <- exp(seq(log(max(a2$lambda)),log(min(a2$lambda)),length.out=num.of.grid))
ptm <- proc.time()
#out <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],eps_abs=.1*eps_augADMM,eps_rel=eps_augADMM,graph_only=TRUE)
out <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,gamma=gamma,data[[6]],data[[7]],eps_abs=.1*eps_augADMM,eps_rel=eps_augADMM)
time.diff3 <- proc.time() - ptm

############ compare suboptimality against genlasso #######
beta.out = coef(a2, lambda_graph)$beta
fun.val.genlasso <- obj_val(beta.out)

fun.val.aADMM <- obj_val(out$beta_path)
#save(beta, file="beta.augADMM.RData")

rel_diff_augADMM = (fun.val.genlasso - fun.val.aADMM)/abs(fun.val.aADMM)

#plot(rel_diff_augADMM,ylab="relative objective function differences",xlab="lambda")


######### run standard ADMM again for 1000 lambdas ###########
ptm <- proc.time()
#out <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],eps_abs=.1*eps_stanADMM,eps_rel=eps_stanADMM,graph_only=TRUE,standard_ADMM=TRUE)
out <- linreg_path_v2(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,gamma=gamma,data[[6]],data[[7]],eps_abs=.1*eps_stanADMM,eps_rel=eps_stanADMM,standard_ADMM=TRUE)
time.diff5 <- proc.time() - ptm

fun.val.stanADMM <- obj_val(out$beta_path)
#save(beta, file="beta.stanADMM.RData")

rel_diff_stanADMM = (fun.val.genlasso - fun.val.stanADMM)/abs(fun.val.stanADMM)
#plot(rel_diff_stanADMM,ylab="relative objective function differences",xlab="lambda")
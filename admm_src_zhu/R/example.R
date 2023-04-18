require(xtable)
source("R/opt.R")
source("R/gen_data.R")
n <- 200
num.vars.per.group <- 11
num.active.groups <- 4
err.var <- .1
out.table <- matrix(0,24,3)
# row_names <- rep(c("accelerated special ADMM","accelerated ADMM","special ADMM","ADMM"),3)
# col_names <- rep(c("time(s)","#iter","#chol"),2)
index <- 0
lambda_graph <-  exp(seq(-log(1e3),log(1e3),length.out=10))
lambda_sparse <-  exp(seq(-log(1e3),log(1e3),length.out=10))

for (cor in c(.5,.7,.9)){
	for (num.groups in c(200,300)){
		set.seed(19630103)
		data <- gen_data(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var)
    	for (variant in c(TRUE,FALSE)) {
        	for (nonvary_rho in c(TRUE,FALSE)) {
                index <- index + 1
                ptm <- proc.time()
                out <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
                lambda_graph,lambda_sparse,data[[6]],data[[7]],variantADMM=variant,nonadaptive_varying_rho=nonvary_rho)
                time.diff <- proc.time() - ptm
                out.table[index,] <- round(c(time.diff[3],sum(out$chol_num),mean(out$admm_iter_num)),digit=3)
                #names(summary) <- c("time(s)","#iter","#chol")
                #filename = paste("variant",toString(variant),"groups",toString(num.groups),"cor",toString(cor),".RDtata",sep="")
                #save(summary,file=filename)
            }
        }
    }
}
save(out.table,file="example.RData")
out.latex.src <- xtable(out.table)
save(out.latex.src,file="example_table.RData")

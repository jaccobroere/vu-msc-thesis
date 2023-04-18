### get results for constant rho ###
require(xtable)
source("R/opt.R")
source("R/gen_data.R")
set.seed(19630103)
n <- 200
num.vars.per.group <- 11
num.active.groups <- 4
err.var <- .1
# row_names <- rep(c("accelerated special ADMM","accelerated ADMM","special ADMM","ADMM"),3)
# col_names <- rep(c("time(s)","#iter","#chol"),2)
i <- j <- k <- l <- 0
lambda_graph <-  exp(seq(-log(1e3),log(1e3),length.out=10))
lambda_sparse <-  exp(seq(-log(1e3),log(1e3),length.out=10))
rho_seq <- exp(seq(-log(1e3),log(1e3),length.out=20))
cor <- .7
num.groups <- 200
variant <- FALSE
set.seed(19630103)
data <- gen_data(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var)

out.table <- matrix(0,length(rho_seq),3)

for (l in 1:length(rho_seq)) {
    rho <- rho_seq[l]
    ptm <- proc.time()
    out <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
    lambda_graph,lambda_sparse,data[[6]],data[[7]],rho=rho,variantADMM=variant,constant_rho=TRUE)
    time.diff <- proc.time() - ptm
    out.table[l,] <- round(c(time.diff[3],sum(out$chol_num),mean(out$admm_iter_num)),digit=1)
    #names(summary) <- c("time(s)","#iter","#chol")
    #filename = paste("variant",toString(variant),"groups",toString(num.groups),"cor",toString(cor),".RDtata",sep="")
    #save(summary,file=filename)
}


save(out.table,file="example_constant_rho.RData")
out.latex.src <- xtable(out.table)
save(out.latex.src,file="example_constant_rho_table.RData")

require(xtable)
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

lambda_graph <-  exp(seq(-log(1e3),log(1e3),length.out=10))
lambda_sparse <-  exp(seq(-log(1e3),log(1e3),length.out=10))
#lambda_graph <- c(1e-2)
#lambda_sparse <- c(1)
data <- gen_data(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var)
ptm <- proc.time()
out <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
                    lambda_graph,lambda_sparse,
                    data[[6]],data[[7]],linearized_ADMM=TRUE)
time.diff <- proc.time() - ptm

error <- matrix(0,length(lambda_graph),length(lambda_sparse))
for (i in 1:length(lambda_graph))
{
    for (j in 1:length(lambda_sparse))
    {
        beta.current <- out$beta_path[,j,i]
        error[i,j] <- sqrt(sum((beta.current-data[[8]])^2) / p)
    }
}

index <- which.min(error)
j.opt <- ceiling(index / length(lambda_graph))
i.opt <- index - (j.opt-1)*length(lambda_graph)

lambda_graph_single <- lambda_graph[i.opt]
lambda_sparse_single <- lambda_sparse[j.opt]

out_single <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
lambda_graph_single,lambda_sparse_single,
data[[6]],data[[7]])
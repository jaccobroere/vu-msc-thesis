#sample size 
source("R/opt.R")
source("R/gen_data.R")
set.seed(19630103)
n <- 1000
num.vars.per.group <- 11
num.active.groups <- 4
err.var <- 1

        for (num.groups in c(200,300)){
			p <- num.groups * num.vars.per.group
                for (cor in c(.5,.7,.9)){
                        data <- gen_data(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var)
						lambda_graph <-  exp(seq(-log(1e3),log(1e3),length.out=10))
						lambda_sparse <-  exp(seq(-log(1e3),log(1e3),length.out=10))
                        ptm <- proc.time()
                        out <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
							lambda_graph,lambda_sparse,data[[6]],data[[7]])
                        time.diff <- proc.time() - ptm
                        summary <- c(time.diff,mean(out$admm_iter_num),sum(out$chol_num))
                        filename = paste("special_ADMM_groups",toString(num.groups),"cor",toString(cor),"accelerated.RDtata",sep="")
                        save(summary,file=filename)
                }
        }

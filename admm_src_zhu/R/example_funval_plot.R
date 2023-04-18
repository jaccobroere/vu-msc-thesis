source("R/opt.R")
source("R/gen_data.R")
set.seed(19630103)
n <- 200
num.vars.per.group <- 11
num.active.groups <- 4
err.var <- .1
for (num.groups in c(200,300)){
    for (cor in c(.5,.7,.9)) {
        
        p <- num.groups * num.vars.per.group
        lambda_graph <-  exp(seq(-log(1e3),log(1e3),length.out=10))
        lambda_sparse <-  exp(seq(-log(1e3),log(1e3),length.out=10))
        #lambda_graph <- c(1e-2)
        #lambda_sparse <- c(1)
        data <- gen_data(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var)
        out <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
        lambda_graph,lambda_sparse,
        data[[6]],data[[7]])
        
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
        
        
        out_single_long <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
        lambda_graph_single,lambda_sparse_single,
        data[[6]],data[[7]],max_num_iter=2e4,
        eps_abs=.0, eps_rel=.0)
        
        opt_val <- out_single_long$fun_val[1]
        
        
        out.acc.specialADMM <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
        lambda_graph_single,lambda_sparse_single,
        data[[6]],data[[7]],max_num_iter=2e3,reporting=TRUE)
        fun.path.acc.specialADMM <- log10(out.acc.specialADMM$fun_path - opt_val)
        #save(fun.path.acc.specialADMM,file="fun.acc.specialADMM.RData")
        
        
        out.acc.ADMM <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
        lambda_graph_single,lambda_sparse_single,
        data[[6]],data[[7]],max_num_iter=2e3,variantADMM=FALSE,reporting=TRUE)
        fun.path.acc.ADMM <- log10(out.acc.ADMM$fun_path  - opt_val)
        #save(fun.path.acc.ADMM,file="fun.acc.ADMM.RData")
        
        out.nonacc.specialADMM <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
        lambda_graph_single,lambda_sparse_single,
        data[[6]],data[[7]],max_num_iter=2e3,nonadaptive_varying_rho=TRUE,reporting=TRUE)
        fun.path.nonacc.specialADMM <- log10(out.nonacc.specialADMM$fun_path  - opt_val)
        #save(fun.path.nonacc.specialADMM,file="fun.nonacc.specialADMM.RData")
        
        
        out.nonacc.ADMM <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],
        lambda_graph_single,lambda_sparse_single,
        data[[6]],data[[7]],max_num_iter=2e3,nonadaptive_varying_rho=TRUE,variantADMM=FALSE,reporting=TRUE)
        fun.path.nonacc.ADMM <- log10(out.nonacc.ADMM$fun_path - opt_val)
        #save(fun.path.nonacc.ADMM,file="fun.nonacc.ADMM.RData")
        
        
        pdf(paste("fun_val","cor",toString(cor*10),"num_groups",toString(num.groups),".pdf",sep=""))
        # postscript(paste("fun_val","cor",toString(cor*10),"num_groups",toString(num.groups),".eps",sep=""))
        plot(fun.path.acc.specialADMM,type="l",col="blue",xlab="",ylab="",ylim=c(-10,2))
        mtext(expression(log(p^k - p^N)),2,2)
        mtext(paste("iteration",expression(k)),1,2)
        
        lines(fun.path.acc.ADMM,type="l",col="red")
        lines(fun.path.nonacc.specialADMM,type="l",col="blue",lty=2)
        lines(fun.path.nonacc.ADMM,type="l",col="red",lty=2)
        legend(800,1.5, c("accelerated augmented ADMM","accelerated standard ADMM","augmented ADMM","standard ADMM"),
        col = c("blue","red","blue","red"), lty=c(1,1,2,2),bty="n",cex=.8)
        dev.off()
        
    }
}

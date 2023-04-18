source("R/opt.R")
source("R/gen_data.R")
n <- 200
num.vars.per.group <- 11
num.active.groups <- 4
err.var <- .1
cor <- .7
num.groups <- 300
p <- num.groups*num.vars.per.group
set.seed(19630103)
data <- gen_data(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var)
#D = sparseMatrix(i=data$idx,j=data$jdx,x=data$val)
lambda_graph <-  c(.05)
lambda_sparse <- c(.5)

############ running Augmented ADMM #############
out.report <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],reporting=TRUE,max_num_iter=2000)

############# running standard ADMM ##############
out.report.stand <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],standard_ADMM=TRUE,reporting=TRUE,max_num_iter=2000)


############ get approximate optimal value ##############
#out.tmp <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],eps_abs=1e-20,eps_rel=1e-19)
out.tmp <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],standard_ADMM=TRUE,reporting=TRUE,max_num_iter=3000)
fun.val.star = min(out.tmp$fun_path)

options(digits=20)
num.of.iter = 1700
fun_val_fgfl <- as.matrix(read.table(paste("tmp/funVal_p",p,".txt",sep=""))) - fun.val.star
plot(1:num.of.iter,log10(out.report$fun_path[1:num.of.iter]-fun.val.star),type="l",col="blue",ylim=c(-12,4),ylab="",xlab="",main=paste("p=",p,sep=""))
lines(1:num.of.iter,log10(out.report.stand$fun_path[1:num.of.iter]-fun.val.star),type="l",col="orange")
lines(log10(fun_val_fgfl[1:num.of.iter]),type="l",col="red")
legend(1000,3,c("Augmented ADMM", "Standard ADMM", "genlasso"), col = c("blue","orange","red"), lty=c(1,1,1),bty="n",cex=.8)
#legend(1400,3,c("Augmented ADMM", "genlasso"), col = c("blue","red"), lty=c(1,1),bty="n",cex=.8)
mtext(expression(log(p^k - p^N)),2,2)
mtext(paste("iteration",expression(k)),1,2)



source("R/opt.R")
source("R/gen_data.R")
n <- 200
num.vars.per.group <- 11
num.active.groups <- 4
err.var <- .1
cor <- .5
num.groups <- 200
p <- num.groups*num.vars.per.group
set.seed(19630103)
data <- gen_data(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var)
D = sparseMatrix(i=data$idx,j=data$jdx,x=data$val)
lambda_graph <-  c(.1)
lambda_sparse <- c(.1)

obj_val <- function(beta)
{
  tmp = norm(data$Y - data$X %*% beta,"F")
  fun.val = .5*tmp*tmp + lambda_sparse[1] * norm(beta,"o") + lambda_graph[1] * norm(D %*% beta,"o")
  fun.val
}


############ running Augmented ADMM #############
ptm <- proc.time()
out <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],eps_abs=1e-17,eps_rel=1e-16)
time.diff2 <- proc.time() - ptm
beta = matrix(out$beta_path[,,1])
fun.val = obj_val(beta)
out.report <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],reporting=TRUE,max_num_iter=2000)

############# running standard ADMM ##############
ptm <- proc.time()
out.stand <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],eps_abs=1e-17,eps_rel=1e-16,standard_ADMM=TRUE)
time.diff3 <- proc.time() - ptm
beta = matrix(out$beta_path[,,1])
fun.val.stand = obj_val(beta)
out.report.stand <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],standard_ADMM=TRUE,reporting=TRUE,max_num_iter=2000)


############ get approximate optimal value ##############
out.tmp <- linreg_path(data[[1]],data[[2]],data[[3]],data[[4]],data[[5]],lambda_graph,lambda_sparse,data[[6]],data[[7]],eps_abs=1e-20,eps_rel=1e-19)
beta = matrix(out.tmp$beta_path[,,1])
fun.val.star = obj_val(beta)

############## save data for fGFL method ################
nonzero_idx <- function(x) { which(x != 0) }
edge_set = apply(D,1,nonzero_idx)
dir <- "../code_fgfl_aaai14/data/"
write.table(edge_set[1,],file=paste(dir,"in_node.txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)
write.table(edge_set[2,],file=paste(dir,"out_node.txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)
write.table(data$X,file=paste(dir,"data_matrix.txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)
write.table(data$Y,file=paste(dir,"response.txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)
write.table(data$idx,file=paste(dir,"D_idx.txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)
write.table(data$jdx,file=paste(dir,"D_jdx.txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)
write.table(data$val,file=paste(dir,"D_val.txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)

options(digits=20)
fun_val_fgfl <- as.matrix(read.table("tmp/funVal.txt")) - fun.val.star
plot(1:2000,log10(out.report$fun_path[1:2000]-fun.val.star),type="l",col="blue",ylim=c(-11,4),ylab="",xlab="",main=paste("p=",p,sep=""))
lines(1:1500,log10(out.report.stand$fun_path[1:1500]-fun.val.star),type="l",col="orange")
lines(log10(fun_val_fgfl),type="l",col="red")
legend(1100,3,c("Augmented ADMM", "Standard ADMM", "genlasso"), col = c("blue","orange","red"), lty=c(1,1,1),bty="n",cex=.8)
#legend(1400,3,c("Augmented ADMM", "genlasso"), col = c("blue","red"), lty=c(1,1),bty="n",cex=.8)
mtext(expression(log(p^k - p^N)),2,2)
mtext(paste("iteration",expression(k)),1,2)



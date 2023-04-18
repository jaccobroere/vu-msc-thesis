source("R/gen_data.R")
require(Matrix)

dir <- "../code_fgfl_aaai14/data/"
TMP = seq(from=20,to=300,by=20)
for (num.groups in TMP){
  cor = .7
  n = 200
  num.vars.per.group <- 11
  num.active.groups <- 4
  err.var <- 20
  p <- num.vars.per.group*num.groups
  set.seed(19630103)
  data <- gen_data(num.groups,num.vars.per.group,n,num.active.groups,cor,err.var)
  D = sparseMatrix(i=data$idx,j=data$jdx,x=data$val)
  nonzero_idx <- function(x) { which(x != 0) }
  edge_set = apply(D,1,nonzero_idx)
  write.table(edge_set[1,],file=paste(dir,"in_node_p",p,".txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)
  write.table(edge_set[2,],file=paste(dir,"out_node_p",p,".txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)
  write.table(data$X,file=paste(dir,"data_matrix_p",p,".txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)
  write.table(data$Y,file=paste(dir,"response_p",p,".txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)
  write.table(data$idx,file=paste(dir,"D_idx_p",p,".txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)
  write.table(data$jdx,file=paste(dir,"D_jdx_p",p,".txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)
  write.table(data$val,file=paste(dir,"D_val_p",p,".txt",sep=""),sep=",",row.name=FALSE,col.names=FALSE)
}
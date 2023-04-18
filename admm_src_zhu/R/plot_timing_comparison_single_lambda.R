require(ggplot2)
require(gridExtra)
require(reshape2)
rel_diff <- function(val1,val2){ 
  tmp <- (val1 - val2) / (abs(val2) + 1e-10) 
  tmp[tmp < .0] = .0
  tmp
}
gamma <- 10
num.of.grid <- 20
p = 3300
dir <- "../code_fgfl_aaai14/data/"
lambda <- read.table(paste(dir,"lambdas_single_p",p,"_gamma_",gamma,".txt",sep=""),col.names="lambda")

######## plot relative difference ##########

load(paste("tmp/funVals_single_p",p,"_gamma_",gamma,".RData",sep=""))
funVals_fGFL = as.numeric(read.table(paste("tmp/fGFL_funVals_p",p,"_gamma_",gamma,".txt",sep=""),sep=","))
#pdf(paste("tmp/reldiff_funval_single_p",p,"_gamma_",gamma,".pdf",sep=""))
rel_augADMM = rel_diff(funVals$fun.val.augADMM,funVals$fun.val.true)+1e-20
rel_stanADMM = rel_diff(funVals$fun.val.stanADMM,funVals$fun.val.true)+1e-20
rel_genlasso = rel_diff(funVals$fun.val.genlasso,funVals$fun.val.true)+1e-20
rel_fGFL = rel_diff(funVals_fGFL,funVals$fun.val.true)+1e-20
rel_funval <- data.frame(id=1:num.of.grid,augADMM=rel_augADMM,stanADMM=rel_stanADMM,fGFL=rel_fGFL,genlasso=rel_genlasso)
mrel <- melt(rel_funval,id.var="id")
mrel[,"lambda"] <- do.call("rbind", replicate(4, lambda, simplify = FALSE))
names(mrel) = c("id", "method","suboptimality","lambda")
p.tmp = ggplot(mrel, aes(lambda,suboptimality,group=method,colour=method,shape=method)) + scale_y_log10() + scale_x_log10()
p.tmp = p.tmp + geom_point() + geom_line() + labs(title=bquote(Problem~dimension==.(p)))
p.tmp = p.tmp + theme(legend.position = "right") #theme(legend.position = c(.75, 0.85),legend.justification = "left")
#p.tmp = p.tmp + theme(legend.text = element_text(face = "bold", size = 15),legend.title = element_text(face = "bold", size = 15))
p.tmp = p.tmp + xlab(expression(lambda)) + ylab("Relative suboptimality")
#p.tmp = p.tmp + theme(axis.text = element_text(size = 15))
#p.tmp = p.tmp + theme(axis.title = element_text(size = 15))
#p.tmp + theme(legend.key.size=unit(.3, "in"))
p.tmp = p.tmp + scale_shape_discrete(name  ="", breaks=c("augADMM", "stanADMM","fGFL","genlasso"), labels=c("augADMM", "stanADMM","fGFL","genlasso"))
p.tmp = p.tmp + scale_colour_discrete(name  ="", breaks=c("augADMM", "stanADMM","fGFL","genlasso"), labels=c("augADMM", "stanADMM","fGFL","genlasso"))
p.tmp
ggsave(p.tmp, file=paste("../plots/single_lam_subopt_p",p,".pdf",sep=""),width=6, height=4.5)


######### plot timing ############
load(paste("tmp/timing_single",p,"_gamma_",gamma,".RData",sep=""))
timing_fGFL = as.numeric(read.table(paste("tmp/fGFL_timing_p",p,"_gamma_",gamma,".txt",sep=""),sep=","))
#pdf(paste("tmp/timing_single",p,"_gamma_",gamma,".pdf",sep=""))
timing <- data.frame(grid=1:num.of.grid,augADMM=runtime$time.augADMM,stanADMM=runtime$time.stanADMM,fGFL=timing_fGFL,genlasso=runtime$time.genlasso)
mtiming <- melt(timing,id.var="grid")
mtiming[,"lambda"] <- do.call("rbind", replicate(4, lambda, simplify = FALSE))
names(mtiming) = c("id","method","time","lambda")
p.tmp = ggplot(mtiming, aes(lambda,time,group=method,colour=method,shape=method)) + scale_y_log10() + scale_x_log10()
p.tmp = p.tmp + geom_point() + geom_line() + labs(title=bquote(Problem~dimension==.(p)))
p.tmp = p.tmp + theme(legend.justification = "right") #legend.position = c(.75, 0.85)
#p.tmp = p.tmp + theme(legend.text = element_text(face = "bold", size = 15),legend.title = element_text(face = "bold", size = 15))
p.tmp = p.tmp + xlab(expression(lambda)) + ylab("Time (sec)")
#p.tmp = p.tmp + theme(axis.text = element_text(size = 15))
#p.tmp = p.tmp + theme(axis.title = element_text(size = 20))
# p.tmp + theme(legend.key.size=unit(.5, "in"))
p.tmp = p.tmp + scale_shape_discrete(name  ="", breaks=c("augADMM", "stanADMM","fGFL","genlasso"), labels=c("augADMM", "stanADMM","fGFL","genlasso"))
p.tmp = p.tmp + scale_colour_discrete(name  ="", breaks=c("augADMM", "stanADMM","fGFL","genlasso"), labels=c("augADMM", "stanADMM","fGFL","genlasso"))
p.tmp
ggsave(p.tmp, file=paste("../plots/single_lam_timing_p",p,".pdf",sep=""),width=6, height=4.5)


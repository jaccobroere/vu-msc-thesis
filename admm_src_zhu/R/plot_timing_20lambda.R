require(reshape2)
require(ggplot2)
subopt=c(1e-2,1e-4,1e-6,1e-8)
EPS = min(subopt)
num.groups = 200
gamma = 5
p <- 11*num.groups
load(paste("tmp/timing_20lambda_",p,"_gamma_",gamma,".RData",sep=""))
load(paste("tmp/iteration_20lambda_",p,"_gamma_",gamma,".RData",sep=""))
dir <- "../code_fgfl_aaai14/data/"
lambda <- as.matrix(read.table(paste(dir,"lambdas_single_p",p,"_gamma_",gamma,".txt",sep="")))

timing.fGFL = as.matrix(read.table(paste("tmp/fGFL_20lambda_timing_p",p,"_gamma_",gamma,".txt",sep=""),sep=","),ncol=length(subopt))
timing.fGFL = as.vector(timing.fGFL)
timing.augADMM = as.vector(runtime$augADMM)
timing.stanADMM = as.vector(runtime$stanADMM)
suboptimality = as.factor(c(rep(1e-2,20),rep(1e-4,20),rep(1e-6,20),rep(1e-8,20)))
lambda_rep = rep(lambda,4)

timing = data.frame(augADMM = timing.augADMM, stanADMM = timing.stanADMM, fGFL= timing.fGFL, suboptimality = suboptimality,lambda=lambda_rep)
mtiming <- melt(timing,id.var=c("suboptimality","lambda"))
fn = function(x)  rep(x,20)
group.var = as.factor(sapply(1:12,function(x) rep(x,20)))
mtiming[,'group'] = group.var
names(mtiming) = c("suboptimality","lambda","method","runtime","group")
p.tmp = ggplot(mtiming, aes(lambda,runtime,group=group)) + scale_y_log10() + scale_x_log10()
p.tmp = p.tmp + geom_point(aes(shape=suboptimality)) 
p.tmp = p.tmp + geom_line(aes(colour=method))
p.tmp = p.tmp +  scale_colour_discrete(name  ="Algorithm",breaks=c("augADMM", "stanADMM","fGFL"),labels=c("augADMM", "stanADMM","fGFL"))
p.tmp= p.tmp+scale_shape_discrete(name  = "Suboptimality",breaks=c(levels(suboptimality)),labels=c(bquote(10^{-8}),bquote(10^{-6}),bquote(10^{-4}),bquote(10^{-2})))
p.tmp = p.tmp + labs(title=bquote(p==.(p)~~~nu==.(gamma)))
p.tmp = p.tmp + xlab(expression(lambda)) + ylab("Time (sec)")
p.tmp
ggsave(p.tmp, file=paste("../plots/single_20lam_timing_p",p,"_gamma_",gamma,".pdf",sep=""),width=6, height=4.5)

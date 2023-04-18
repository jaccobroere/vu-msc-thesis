require(ggplot2)
require(gridExtra)
require(reshape2)
num.groups = 200
scale = 1
gamma = 10
lambda = .01
lambda_graph = lambda * scale 
lambda_sparse = lambda_graph * gamma


p <- 11 * num.groups

############# ggplots #############
file_name = paste("tmp/convergence_funPaths_p",p,"_lam1_",lambda_graph,"_lam2_",lambda_sparse,".RData",sep="")
load(file_name)
max_iter = length(out$funPath.augADMM)
rel_augADMM = out$funPath.augADMM
rel_augADMM2 = out$funPath.augADMM2
rel_augADMM3 = out$funPath.augADMM3
rel_augADMM4 = out$funPath.augADMM4
rel_stanADMM = out$funPath.stanADMM
rel_fGFL = out$funPath.fGFL
rel_funval <- data.frame(id=1:max_iter,augADMM=rel_augADMM,augADMM2=rel_augADMM2,augADMM3=rel_augADMM3,augADMM4=rel_augADMM4,stanADMM=rel_stanADMM,fGFL=c(rel_fGFL,rep(NA,max_iter - length(rel_fGFL))))
mrel <- melt(rel_funval,id.var="id",na.rm=TRUE)
names(mrel) = c("id", "method","suboptimality")
mrel <- mrel[which(mrel$suboptimality > 1e-12),]
p.tmp = ggplot(mrel, aes(id,suboptimality,group=method,colour=method)) + scale_y_log10() + scale_x_log10()
tmp = toString(lambda_graph)
p.tmp = p.tmp + geom_line() + labs(title=bquote(lambda==.(lambda_graph)~~~nu==.(gamma)))
p.tmp = p.tmp + xlab("Iteration") + ylab("Suboptimality")
p.tmp = p.tmp + scale_shape_discrete(name  ="Algorithm", breaks=c("augADMM", "augADMM2", "augADMM3", "augADMM4","stanADMM","fGFL"), labels=c("augADMM", "augADMMx1", "augADMMx5", "augADMMx10","stanADMM","fGFL"))
p.tmp = p.tmp + scale_colour_discrete(name  ="Algorithm", breaks=c("augADMM", "augADMM2", "augADMM3", "augADMM4","stanADMM","fGFL"), labels=c("augADMM", "augADMMx1", "augADMMx5", "augADMMx10","stanADMM","fGFL"))
#p.tmp = p.tmp + scale_shape_discrete(name  ="Algorithm")
#p.tmp = p.tmp + scale_colour_discrete(name  ="Algorithm")
p.tmp = p.tmp + theme(legend.position = c(.15, .3))
p.tmp
ggsave(p.tmp, file=paste("../plots/subopt_plot_lam_",scale,"_gamma_",gamma,".pdf",sep=""),width=6, height=4.5)


# sub_opt <-  c(1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10)
# load(paste("tmp/convergence_timing_seq_p",p,"_lam1_",lambda_graph,"_lam2_",lambda_sparse,".RData",sep=""))
# time.fGFL = as.numeric(read.table(paste("tmp/convergence_timing_p",p,"_lam1_",lambda_graph,"_lam2_",lambda_sparse,".txt",sep=""),sep=","))
# time.augADMM = tmp.time$time.augADMM
# time.stanADMM = tmp.time$time.stanADMM
# df.time <- data.frame(id=1:length(sub_opt), augADMM=time.augADMM,stanADMM=time.stanADMM,fGFL=time.fGFL)
# mdf.time = melt(df.time,id.var="id")
# mdf.time[,"suboptimality"] <- do.call("rbind", replicate(3, as.matrix(sub_opt), simplify = FALSE))
# names(mdf.time) = c("id", "method","time","suboptimality")
# p.tmp = ggplot(mdf.time, aes(time,suboptimality,group=method,colour=method,shape=method)) + scale_y_log10() + scale_x_log10()
# p.tmp = p.tmp + geom_point() + geom_path() + labs(title=bquote(lambda==.(lambda_graph)~~~nu==.(gamma)))
# p.tmp = p.tmp + xlab("Time (sec)") + ylab("Relative suboptimality")
# p.tmp = p.tmp + scale_shape_discrete(name  ="", breaks=c("augADMM", "stanADMM","fGFL"), labels=c("augADMM", "stanADMM","fGFL"))
# p.tmp = p.tmp + scale_colour_discrete(name  ="", breaks=c("augADMM", "stanADMM","fGFL"), labels=c("augADMM", "stanADMM","fGFL"))
# ggsave(p.tmp, file=paste("../plots/time_subopt_lam_",scale,".pdf",sep=""),width=6, height=4.5)
# 

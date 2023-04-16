num.groups = seq(20,300,by=20)
p = 11 * num.groups

########## results for gamma = 10 ##########
genlasso = c(108.606,312.825,527.217,532.808,1244.996,1536.324,1953.478,1952.781,2084.38,3044.548,3001.943,3610.79,4324.766,5596.162,5431.095)
augADMM = c(10.707,28.526,25.985,26.042,68.213,101.338,91.105,90.828,116.756,151.693,151.669,204.032,201.658,208.903,198.443)
stanADMM = c(4.028,10.708,21.696,21.838,57.99,135.99,160.638,159.136,225.552,349.122,348.911,480.713,640.551,719.875,737.876)

df = data.frame(p = p, augADMM = augADMM, stanADMM = stanADMM, genlasso = genlasso)
mdf = melt(df,id.var="p")
names(mdf) = c("dimension", "method","time")
p.tmp = ggplot(mdf, aes(dimension,time,group=method,colour=method,shape=method))   # + scale_y_log10() + scale_x_log10()
p.tmp = p.tmp + geom_point() + geom_line() + labs(title=expression(paste(nu," = 10")))
p.tmp = p.tmp + xlab("Problem dimension") + ylab("Time (sec)")
p.tmp = p.tmp + scale_shape_discrete(name  ="Algorithm", breaks=c("augADMM", "stanADMM","genlasso"), labels=c("augADMM", "stanADMM","genlasso"))
p.tmp = p.tmp + scale_colour_discrete(name  ="Algorithm", breaks=c("augADMM", "stanADMM","genlasso"), labels=c("augADMM", "stanADMM","genlasso"))
p.tmp = p.tmp + theme(legend.position = c(.15, .8))
p.tmp
ggsave(p.tmp, file="../plots/plot_timing_path_gamma10.pdf",width=6, height=4.5)

######### results for gamma = 5 #########
genlasso = c(127.751,365.016,665.907,977.036,1471.148,1701.768,2081.172,2563.679,3034.94,4062.789,3711.112,3999.143,4424.267,5278.318,5416.328)
augADMM = c(12.512,26.555,21.844,53.025,69.634,98.968,84.718,123.756,101.353,147.727,181.486,169.315,180.005,184.74,187.869)
stanADMM = c(4.453,10.181,20.56,36.176,57.697,125.116,151.806,219.785,222.169,333.723,397.433,453.486,539.081,592.789,699.594)

df = data.frame(p = p, augADMM = augADMM, stanADMM = stanADMM, genlasso = genlasso)
mdf = melt(df,id.var="p")
names(mdf) = c("dimension", "method","time")
p.tmp = ggplot(mdf, aes(dimension,time,group=method,colour=method,shape=method))   # + scale_y_log10() + scale_x_log10()
p.tmp = p.tmp + geom_point() + geom_line() + labs(title=expression(paste(nu," = 5")))
p.tmp = p.tmp + xlab("Problem dimension") + ylab("Time (sec)")
p.tmp = p.tmp + scale_shape_discrete(name  ="Algorithm", breaks=c("augADMM", "stanADMM","genlasso"), labels=c("augADMM", "stanADMM","genlasso"))
p.tmp = p.tmp + scale_colour_discrete(name  ="Algorithm", breaks=c("augADMM", "stanADMM","genlasso"), labels=c("augADMM", "stanADMM","genlasso"))
p.tmp = p.tmp + theme(legend.position = c(.15, .8))
p.tmp
ggsave(p.tmp, file="../plots/plot_timing_path_gamma5.pdf",width=6, height=4.5)

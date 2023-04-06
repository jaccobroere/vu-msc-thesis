library(splash)
library(genlasso)
library(data.table)

setwd("/Users/jacco/Documents/repos/vu-msc-thesis/")

y <- t(fread("out/sim_y.csv", header = T, skip = 0))

model <- splash(y)

model$AB

#!/bin/bash
# script running three methods at a single value of tuning parameter
# it reports the function paths and vectors of runtime and sub-optimality
######## used to produce Figure 2 ###########
cd ../code_fgfl_aaai14/
mcc -R -singleCompThread -R -nodesktop -R -nodisplay -R -nosplash -R -logfile -R output.log -m comparison_convergence.m -a ./GFL
./run_comparison_convergence.sh /usr/local/matlab_r2015a  200 .1 10 2e3 &
wait
cd ../
make clean
make
#R CMD BATCH --no-save --no-restore '--args num.groups=20 lambda_graph=.05 gamma=10 max_iter=1e3' R/convergence.R convergence_out &
R CMD BATCH --no-save --no-restore '--args num.groups=200 lambda_graph=0.1 gamma=10 max_iter=1e3' R/convergence_single_subopt.R convergence_single_out &
#R CMD BATCH --no-save --no-restore '--args num.groups=200 lambda_graph=0.01 gamma=5 max_iter=1e3' R/convergence_single_subopt.R convergence_single_out &
#R CMD BATCH --no-save --no-restore '--args num.groups=200 lambda_graph=1 gamma=5 max_iter=1e3' R/convergence_single_subopt.R convergence_single_out &
#R CMD BATCH --no-save --no-restore '--args num.groups=200 lambda_graph=.1 gamma=10 max_iter=1e3' R/convergence.R convergence_out &
#R CMD BATCH --no-save --no-restore '--args num.groups=200 lambda_graph=.05 gamma=10 max_iter=1500' R/convergence.R convergence_out &
#R CMD BATCH --no-save --no-restore '--args num.groups=200 lambda_graph=.01 gamma=10 max_iter=10000' R/convergence.R convergence_out &
#R CMD BATCH --no-save --no-restore '--args num.groups=200 lambda_graph=.001 gamma=10 max_iter=5e3' R/convergence.R convergence_out &
wait
echo "finish running three methods!"

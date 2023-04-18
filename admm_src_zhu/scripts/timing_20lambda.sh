#!/bin/bash
######## used to produce Figure 1 ###########
cd ../
make clean
make
R CMD BATCH --no-save --no-restore '--args num.groups=100 gamma=10' R/timing_comparison_single_lambda_v2.R timing_single_out &
wait
cd code_fgfl_aaai14/
mcc -R -singleCompThread -R -nodesktop -R -nodisplay -R -nosplash -R -logfile -R output.log -m timing_comparison_20lambda.m -a ./GFL
./run_timing_comparison_20lambda.sh /usr/local/matlab_r2015a  100 10 &
wait
echo "finished running three methods!"

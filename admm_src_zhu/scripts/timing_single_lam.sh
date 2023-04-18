#!/bin/bash
# script running three methods to check their runtime and
# sub-optimalities
cd ..
make clean
make
R CMD BATCH --no-save --no-restore '--args num.groups=300 gamma=10' R/timing_comparison_single_lambda.R timing_single_out &
wait
cd code_fgfl_aaai14/
matlab -nodesktop -nodisplay -singleCompThread < timing_comparison_single_lambda.m &> out &

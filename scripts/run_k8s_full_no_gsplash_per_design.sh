#!/bin/bash

cd $PROJ_DIR

# Define your list of sim_design_ids
sim_design_ids=("designB_T500_p9" "designB_T1000_p9" "designB_T2000_p9")

# Loop over the array
for sim_design_id in "${sim_design_ids[@]}"; do
    # Run the MC simulation job without GSPLASH (i.e. false in 2nd argument)
    bash scripts/run_k8s_mc_simulation.sh $sim_design_id "false" 
done

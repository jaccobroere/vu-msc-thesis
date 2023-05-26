#!/bin/bash

cd $PROJ_DIR

# Define your list of sim_design_ids
sim_design_ids=("designC_T500_p36" "designC_T1000_p36" "designB_T500_p36" "designB_T1000_p36")

# Loop over the array
for sim_design_id in "${sim_design_ids[@]}"; do
    # Run the MC simulation job with all models and full cross-validation
    bash scripts/run_k8s_mc_simulation.sh $sim_design_id "fullcv" 
done
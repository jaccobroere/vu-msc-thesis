#!/bin/bash

cd $PROJ_DIR

# Define your list of sim_design_ids
sim_design_ids=("designB_T500_p49" "designB_T1000_p49" "designB_T2000_p49")

# Loop over the array
for sim_design_id in "${sim_design_ids[@]}"; do
    # Run the determine lambda job, and is skipped if the grid file already exists
    bash scripts/run_k8s_determine_lambda.sh $sim_design_id

    # Run the MC simulation job
    bash scripts/run_k8s_mc_simulation.sh $sim_design_id
done

#!/bin/bash

cd $PROJ_DIR

# Read CLI argument for the design ID
sim_design_id=$1

# Run the determine lambda job, and is skipped if the grid file already exists
bash scripts/run_k8s_determine_lambda.sh $sim_design_id

# Run the MC simulation job
bash scripts/run_k8s_mc_simulation.sh $sim_design_id
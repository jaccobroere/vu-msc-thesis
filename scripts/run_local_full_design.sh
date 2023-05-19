#!/bin/bash
# NOTE: Before running this scripts, make sure that the PV and PVC has been setup
# Use kubectl -f scripts/k8s/setup_pv.yml to set these up if they are not running yet.

cd $PROJ_DIR

# Insert the design ID into the k8s job YML files
sim_design_id=$1
# replace_string='s/REPLACEME/'$sim_design_id'/g'
# sim_design_id_dashes=${sim_design_id//_/-}
# replace_string_dashes='s/MEREPLACE/'${sim_design_id_dashes,,}'/g'

# Start the script
echo "STARTING: Running MC simulation for design $sim_design_id"

# Run the detlam job
for ((i=1; i<=5; i++))
do
    # Task to repeat
    bash scripts/determine_lambda_preliminary.sh $sim_design_id
done

# Run the simulation job
for ((i=1; i<=20; i++))
do
    # Task to repeat
    bash scripts/model_fit_montecarlo.sh $sim_design_id
done


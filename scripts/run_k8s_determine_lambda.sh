#!/bin/bash
# NOTE: Before running this scripts, make sure that the PV and PVC has been setup
# Use kubectl -f scripts/k8s/setup_pv.yml to set these up if they are not running yet.

cd $PROJ_DIR

# Insert the design ID into the k8s job YML files
sim_design_id=$1
grid_dir="$PROJ_DIR/out/simulation/lambdas/${sim_design_id}"

if [ -f "$grid_dir/grid_gfsplash_a05.csv" ]
then
    nlines_grid=$(wc -l < $grid_dir/grid_gfsplash_a05.csv)
    if [ "$nlines_grid" -gt 20 ]
    then
        echo "SKIPPING: Lambda grid search already run for design $sim_design_id"
        exit 0
    fi
fi

# If it has not been run yet, then remove the files from the directory and continue
echo "STARTING: Running lambda grid search for design $sim_design_id"
# Remove the files from the directory
rm -r $grid_dir/*

# Insert the design ID into the k8s job YML files
replace_string='s/REPLACEME/'$sim_design_id'/g'
replace_string_dashes='s/REPLACEME/'${sim_design_id//_/-}'/g'
sed -E $replace_string scripts/k8s/determine_lambda_TEMPLATE.yml > scripts/k8s/determine_lambda_REPLACED.yml
sed -E $replace_string_dashes scripts/k8s/determine_lambda_REPLACED.yml

# Update docker image scripts
# docker build -t jaccusaurelius/vu-msc-thesis:kube .

# Clear running pods and jobs
kubectl apply -f scripts/k8s/setup_pv.yml

# Run the determine_lambda job
kubectl apply -f scripts/k8s/determine_lambda_REPLACED.yml
kubectl wait --for=condition=complete --timeout=5h job/detlam

# Delete the determine_lambda job
kubectl delete -f scripts/k8s/determine_lambda_REPLACED.yml
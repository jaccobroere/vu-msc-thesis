#!/bin/bash
# NOTE: Before running this scripts, make sure that the PV and PVC has been setup
# Use kubectl -f scripts/k8s/setup_pv.yml to set these up if they are not running yet.

cd $PROJ_DIR

# Insert the design ID into the k8s job YML files
sim_design_id=$1
SOURCE="data-access-nfs:/app/out/simulation/lambdas/${sim_design_id}/grid_gfsplash_a05.csv"
DESTINATION="$PROJ_DIR/out/simulation/lookup/${sim_design_id}/grid_gfsplash_a05.csv"
# Attempt to copy the file
kubectl cp $SOURCE $DESTINATION 2>/dev/null
# Check the exit status of the previous command
if [ $? -eq 0 ]; then
  echo "STARTING: Running lambda grid search for design $sim_design_id"
else
  echo "SKIPPING: Lambda grid search already run for design $sim_design_id"
  exit 0
fi

# Remove the files from the directory
kubectl exec -it data-access-nfs -- bash -c "rm -r /app/out/simulation/lambdas/${sim_design_id}/*"

# Insert the design ID into the k8s job YML files
replace_string='s/REPLACEME/'$sim_design_id'/g'
sim_design_id_dashes=${sim_design_id//_/-}
replace_string_dashes='s/MEREPLACE/'${sim_design_id_dashes,,}'/g'
sed -E $replace_string scripts/k8s/determine_lambda_TEMPLATE.yml > scripts/k8s/determine_lambda_REPLACED.yml.tmp
sed -E $replace_string_dashes scripts/k8s/determine_lambda_REPLACED.yml.tmp > scripts/k8s/determine_lambda_REPLACED.yml 
rm scripts/k8s/determine_lambda_REPLACED.yml.tmp

# Update docker image scripts
# docker build -t jaccusaurelius/vu-msc-thesis:kube .

# Clear running pods and jobs
# kubectl apply -f scripts/k8s/setup_pv.yml

# Run the determine_lambda job
kubectl apply -f scripts/k8s/determine_lambda_REPLACED.yml
kubectl wait --for=condition=complete --timeout=5h "job/detlam-${sim_design_id_dashes,,}"

# Delete the determine_lambda job
kubectl delete -f scripts/k8s/determine_lambda_REPLACED.yml
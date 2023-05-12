#!/bin/bash
# NOTE: Before running this scripts, make sure that the PV and PVC has been setup
# Use kubectl -f scripts/k8s/setup_pv.yml to set these up if they are not running yet.

cd $PROJ_DIR

# Insert the design ID into the k8s job YML files
sim_design_id=$1
replace_string='s/REPLACEME/'$sim_design_id'/g'
sed -E $replace_string scripts/k8s/determine_lambda_TEMPLATE.yml > scripts/k8s/determine_lambda_REPLACED.yml

# Update docker image scripts
docker build -t jaccusaurelius/vu-msc-thesis:kube .

# Clear running pods and jobs
# kubectl delete pods --all
# kubectl delete jobs --all

# # Setup pod to access the data in the PV
# kubectl apply -f scripts/k8s/data_access.yml
# kubectl wait --for=condition=ready --timeout=1m pod/data-access
kubectl apply -f scripts/k8s/setup_pv.yml

# Run the determine_lambda job
kubectl apply -f scripts/k8s/determine_lambda_REPLACED.yml
kubectl wait --for=condition=complete --timeout=5h job/detlam

# Delete the determine_lambda job
kubectl delete -f scripts/k8s/determine_lambda_REPLACED.yml
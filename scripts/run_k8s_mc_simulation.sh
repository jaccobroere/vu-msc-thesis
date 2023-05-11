#!/bin/bash
# NOTE: Before running this scripts, make sure that the PV and PVC has been setup
# Use kubectl -f scripts/k8s/setup_pv.yml to set these up if they are not running yet.

cd $PROJ_DIR

# Insert the design ID into the k8s job YML files
sim_design_id=$1
replace_string='s/REPLACEME/'$sim_design_id'/g'
sed -E $replace_string scripts/k8s/mc_simulation_TEMPLATE.yml > scripts/k8s/mc_simulation_REPLACED.yml

# Update docker image
docker build --quiet -t jaccusaurelius/vu-msc-thesis:kube .

# Clear running pods, jobs and delete PVCs
kubectl delete pods --all
kubectl delete jobs --all

# Setup pod to access the data in the PV
kubectl apply -f scripts/k8s/data_access.yml
kubectl wait --for=condition=ready --timeout=1m pod/data-access

# Run the simulation job
kubectl apply -f scripts/k8s/mc_simulation_REPLACED.yml
kubectl wait --for=condition=complete --timeout=12h job/modelfit

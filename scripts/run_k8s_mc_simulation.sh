#!/bin/bash

cd $PROJ_DIR

# Update docker image
docker build --quiet -t jaccusaurelius/vu-msc-thesis:kube .
# Clear running pods, jobs and delete PVCs and PVs
bash scripts/clear_k8s.sh
# Setup the PV and PVC
kubectl apply -f scripts/k8s/setup_pv.yml
# Setup pod to access the data in the PV
kubectl apply -f scripts/k8s/data_access.yml
kubectl wait --for=condition=ready --timeout=30s pod/data_access
# Run the simulation job
kubectl apply -f scripts/k8s/mc_simulation.yml

# kubectl wait --for=condition=complete job/mc-sim

#!/bin/bash

cd $PROJ_DIR

# Update docker image
docker build --quiet -t jaccusaurelius/vu-msc-thesis:kube .
# Clear running pods, jobs and delete PVCs
kubectl delete pods --all
kubectl delete jobs --all
kubectl delete pvc --all
# Setup pod to access the data in the PV
kubectl apply -f scripts/k8s/data_access.yml
kubectl wait --for=condition=ready --timeout=30s pod/data-access
# Run the simulation job
kubectl apply -f scripts/k8s/mc_simulation.yml

kubectl wait --for=condition=complete --timeout=3h job/modelfit

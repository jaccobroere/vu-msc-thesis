#!/bin/bash

cd $PROJ_DIR

sim_design_id=$1
replace_string='s/REPLACEME/'$sim_design_id'/g'
sed -E $replace_string scripts/k8s/determine_lambda.yml > scripts/k8s/determine_lambda_replaced.yml

# Update docker image scripts
docker build -t jaccusaurelius/vu-msc-thesis:kube .

kubectl delete pods --all
kubectl delete jobs --all

# bash scripts/clear_k8s.sh

# kubectl apply -f scripts/k8s/setup_pv.yml

kubectl apply -f scripts/k8s/data_access.yml

kubectl wait --for=condition=ready --timeout=1m pod/data-access

kubectl apply -f scripts/k8s/determine_lambda_replaced.yml

kubectl wait --for=condition=complete --timeout=30m job/detlam

kubectl exec -it data-access -- bash -c "python3 /app/vu-msc-thesis/src/compute/save_best_lambda.py $sim_design_id" 
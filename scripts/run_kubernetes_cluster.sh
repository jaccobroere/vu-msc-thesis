#!/bin/bash

cd $PROJ_DIR

# Update docker image scripts
docker build -t --quiet jaccusaurelius/vu-msc-thesis:kube .

bash scripts/clear_k8s.sh

kubectl apply -f scripts/k8s/setup_pv.yml

kubectl apply -f scripts/k8s/dataaccess.yml

kubectl wait --for=condition=ready --timeout=30s pod/dataaccess

kubectl apply -f scripts/k8s/determine_lambda.yml

kubectl wait --for=condition=complete --timeout=30m job/detlam

kubectl exec -it dataaccess -- bash -c "python3 src/compute/save_best_lambda.py"

kubectl wait --timeout=10s

kubectl cp dataaccess:app/vu-msc-thesis/out/simulation k8s_export

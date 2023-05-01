#!/bin/bash

cd $PROJ_DIR

# Update docker image scripts
docker build --quiet -t jaccusaurelius/vu-msc-thesis:kube .

bash scripts/clear_k8s.sh

kubectl apply -f scripts/k8s/setup_pv.yml

kubectl apply -f scripts/k8s/data_access.yml

kubectl wait --for=condition=ready --timeout=30s pod/data-access

kubectl apply -f scripts/k8s/determine_lambda.yml

kubectl wait --for=condition=complete --timeout=30m job/detlam

kubectl exec -it data-access -- bash -c "python3 src/compute/save_best_lambda.py"

# kubectl wait --timeout=10s # Probably not needed anymore

# kubectl cp data-access:/app/vu-msc-thesis/out/simulation k8s_export # Probably not needed anymore because of host mount of PV

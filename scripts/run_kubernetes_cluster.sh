#!/bin/bash

cd $PROJ_DIR

bash scripts/clear_k8s.sh

kubectl apply -f scripts/k8s/setup_pv.yml

kubectl apply -f scripts/k8s/dataaccess.yml

kubectl wait --for=condition=ready --timeout=30s pod/dataaccess

kubectl apply -f scripts/k8s/determine_lambda.yml

kubectl wait --for=condition=complete --timeout=30m job/detlam

kubectl cp dataaccess:app/vu-msc-thesis/out/simulation k8s_export

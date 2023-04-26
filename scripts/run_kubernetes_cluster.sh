cd $PROJ_DIR

kubectl apply -f scripts/k8s/setup_pv.yml

kubectl apply -f scripts/k8s/determine_lambda.yml

kubectl wait --for=condition=complete --timeout=30m job/detlam

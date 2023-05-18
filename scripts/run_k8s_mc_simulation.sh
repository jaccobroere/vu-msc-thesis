#!/bin/bash
# NOTE: Before running this scripts, make sure that the PV and PVC has been setup
# Use kubectl -f scripts/k8s/setup_pv.yml to set these up if they are not running yet.

cd $PROJ_DIR
# Insert the design ID into the k8s job YML files
sim_design_id=$1
rungsplash=${2:-true}
replace_string='s/REPLACEME/'$sim_design_id'/g'
sim_design_id_dashes=${sim_design_id//_/-}
replace_string_dashes='s/MEREPLACE/'${sim_design_id_dashes,,}'/g'

if [ $rungsplash = "true" ]; then
    sed -E $replace_string scripts/k8s/mc_simulation_TEMPLATE.yml > scripts/k8s/mc_simulation_REPLACED.yml.tmp
    sed -E $replace_string_dashes scripts/k8s/mc_simulation_REPLACED.yml.tmp > scripts/k8s/mc_simulation_REPLACED.yml 
    rm scripts/k8s/mc_simulation_REPLACED.yml.tmp
else
    sed -E $replace_string scripts/k8s/mc_simulation_no_gsplash_TEMPLATE.yml > scripts/k8s/mc_simulation_REPLACED.yml.tmp
    sed -E $replace_string_dashes scripts/k8s/mc_simulation_REPLACED.yml.tmp > scripts/k8s/mc_simulation_REPLACED.yml 
    rm scripts/k8s/mc_simulation_REPLACED.yml.tmp
fi

# Start the script
echo "STARTING: Running MC simulation for design $sim_design_id"

# Update docker image
# docker build --quiet -t jaccusaurelius/vu-msc-thesis:kube .

# Run the simulation job
kubectl apply -f scripts/k8s/mc_simulation_REPLACED.yml
kubectl wait --for=condition=complete --timeout=24h "job/modelfit-${sim_design_id_dashes,,}"

# Delete the modelfit job
kubectl delete -f scripts/k8s/mc_simulation_REPLACED.yml

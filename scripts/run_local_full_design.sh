#!/bin/bash
# NOTE: Before running this scripts, make sure that the PV and PVC has been setup
# Use kubectl -f scripts/k8s/setup_pv.yml to set these up if they are not running yet.

cd $PROJ_DIR

# Insert the design ID into the k8s job YML files
# sim_design_id=$1
sim_design_ids=(
    "designB_T500_p25"
    "designB_T500_p36"
    "designB_T1000_p25"
    "designB_T1000_p36"
    "designB_T2000_p25"
    "designB_T2000_p36"
    "designC_T500_p25"
    "designC_T500_p36"
    "designC_T1000_p25"
    "designC_T1000_p36"
    "designC_T2000_p25"
    "designC_T2000_p36"
    "designD_T500_p25"
    "designD_T500_p36"
    "designD_T1000_p25"
    "designD_T1000_p36"
    "designD_T2000_p25"
    "designD_T2000_p36"
)

# rundetlam=${2:-"false"}
# replace_string='s/REPLACEME/'$sim_design_id'/g'
# sim_design_id_dashes=${sim_design_id//_/-}
# replace_string_dashes='s/MEREPLACE/'${sim_design_id_dashes,,}'/g'

# # Run the detlam job
# if [ $rundetlam = "true" ]; then
#     for ((i=1; i<=5; i++))
#     do
#         # Task to repeat
#         bash scripts/determine_lambda_preliminary.sh $sim_design_id
#     done
# fi

# Run the simulation job
for sim_design_id in "${sim_design_ids[@]}"; do
    echo "STARTING: Running MC simulation for design $sim_design_id"
    for ((i=1; i<=250; i++))
    do
        # Task to repeat
        bash scripts/model_fit_montecarlo.sh $sim_design_id
    done
done


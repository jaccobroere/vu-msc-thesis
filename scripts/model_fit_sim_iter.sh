#!/bin/bash

# Progress bar
total_steps=4
current_step=0

print_progress_bar() {
  local current_step=$1
  local total_steps=$2
  local width=$3
  local percentage=$((current_step * 100 / total_steps))
  local completed_chars=$((width * current_step / total_steps))
  local pending_chars=$((width - completed_chars))
  local progress_bar=""

  for ((i = 0; i < completed_chars; i++)); do
    progress_bar+="#"
  done

  for ((i = 0; i < pending_chars; i++)); do
    progress_bar+="-"
  done

  printf "\r[%s] %d%%" "$progress_bar" "$percentage"
  echo ""
} 

step_create_directories () {
  # If folder structure is not present yet, create it
  mkdir -p out/simulation/fit/${sim_design_id}/${uuidtag}
}

step_sim() {
    # Simulation step
    echo "Running $path ..."
    julia --project=$JULIA_DIR $path $p $T $h_A $h_B ${sim_design_id} $uuidtag
    echo "$path completed."
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

# Transform data
step_transform () {
    # Run Julia script for step 1
    echo "Running transform_bootstrap_graph.jl ..."
    julia --project=$JULIA_DIR src/compute/transform_bootstrap_graph.jl ${sim_design_id} $uuidtag
    echo "transform_bootstrap_graph.jl completed."
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

# Calculate performance for each lambda value once
step_detlam () {
    # rm -rf out/simulation/lambdas/${sim_design_id}
    echo "Running determine_lambda.R ..."
    Rscript src/compute/model_fitting.R ${sim_design_id} $uuidtag > /dev/null 2>&1
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

# DESIGN A
# p=25
# T=500
# h_A=3
# h_B=3
# path=src/simulation/simulation_designA.jl
# prefix=designA

# DESIGN B
p=25 # m^2 
T=500
h_A=3 # Placeholder
h_B=3 # Placeholder
path=src/simulation/simulation_designB.jl
prefix=designB
sim_design_id=${prefix}_T${T}_p${p}

uuidtag=$(uuidgen)
step_create_directories && step_sim && step_transform && step_detlam
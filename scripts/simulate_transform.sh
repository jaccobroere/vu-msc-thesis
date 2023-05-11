#!/bin/bash

# Progress bar
total_steps=3
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
    echo "Creating directories ..."
    mkdir -p data/simulation/${sim_design_id}/
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

step_sim() {
    echo "Running $path ..."
    julia --project=$JULIA_DIR $path ${prefix}_T${T}_p${p}
    echo "$path completed."
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

# Transform data
step_transform () {
    # Run Julia script for step 1
    echo "Running precalculations_and_write.jl.jl ..."
    julia --project=$JULIA_DIR src/compute/precalculations_and_write.jl.jl ${prefix}_T${T}_p${p}
    echo "precalculations_and_write.jl.jl completed."
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}


# # DESIGN A
# p=25
# T=500
# path=src/simulation/simulation_designA.jl
# prefix=designA

# DESIGN B
p=49 # m^2 
T=500
path=src/simulation/simulation_designB.jl
prefix=designB
sim_design_id=${prefix}_T${T}_p${p}

step_create_directories && step_sim && step_transform

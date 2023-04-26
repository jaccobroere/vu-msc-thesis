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

step_sim() {
    echo "Running $path ..."
    julia --project=$JULIA_DIR $path $p $T $h_A $h_B ${prefix}_T${T}_p${p}
    echo "$path completed."
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

# Transform data
step_transform () {
    bash scripts/transform_data.sh -prefix ${prefix}_T${T}_p${p}
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

step_detlam () {
    Rscript src/compute/determine_lambda.R ${prefix}_T${T}_p${p}
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

# DESIGN A
p=25
T=500
h_A=3
h_B=3
path=src/simulation/simulation_designA.jl
prefix=designA

# DESIGN B
# p=100 # m^2 
# T=500
# path=src/simulation/simulation_designB.jl
# prefix=designB

step_sim && step_transform && step_detlam



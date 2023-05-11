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
  echo "Creating directories ..."
  mkdir -p out/simulation/lambdas/${sim_design_id}/${uuidtag}
  mkdir -p data/simulation/${sim_design_id}/detlam/${uuidtag}
  current_step=$((current_step+1))
  print_progress_bar $current_step $total_steps 50
}

step_sim() {
    # Simulation of data step
    echo "Running $path ..."
    julia --project=$JULIA_DIR $path ${sim_design_id}/detlam $uuidtag
    echo "$path completed."
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

# Transform data
step_transform () {
    # Run the precalculation script only if it has not been done before for this design
    if [ ! -f "data/simulation/${sim_design_id}/detlam/Dtilde.mtx" ]; then
      echo "Running precalculations_and_write.jl ..."
      julia --project=$JULIA_DIR src/compute/precalculations_and_write.jl ${sim_design_id}/detlam $uuidtag
      echo "precalculations_and_write.jl completed."
    else 
      echo "precalculations_and_write.jl already completed at an earlier iteration for ${sim_design_id}."
    fi
    echo "precalculations_and_write.jl completed."
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

# Calculate performance for each lambda value once
step_detlam () {
    echo "Running determine_lambda.R ..."
    Rscript src/compute/determine_lambda.R ${sim_design_id} $uuidtag > /dev/null 2>&1
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

# Read in the arguments and parse it
inputarg=$1
design=$(echo $inputarg | sed -E 's/^([a-zA-Z]+)_T[0-9]+_p[0-9]+$/\1/')
T=$(echo $inputarg | sed -E 's/^[a-zA-Z]+_T([0-9]+)_p[0-9]+$/\1/')
p=$(echo $inputarg| sed -E 's/^[a-zA-Z]+_T[0-9]+_p([0-9]+)$/\
1/')
path=src/simulation/simulation_${design}.jl
sim_design_id=${inputarg}

uuidtag=$(uuidgen)
step_create_directories && step_sim && step_transform && step_detlam
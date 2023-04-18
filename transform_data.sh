#!/bin/bash

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

# Parse named arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
  -path_prefix)
    path_prefix="$2"
    shift # past argument
    shift # past value
    ;;
  *)
    echo "Unknown option: $1"
    exit 1
    ;;
  esac
done

# Define step 1 function
step1() {
  # Check if step1_script argument was provided
  if [[ -z "$path_prefix" ]]; then
    echo "Error: step1_script argument is missing"
    exit 1
  fi

  # Run Julia script for step 1
  echo "Running construct_transformation.jl ..."
  julia --project=/path/to/juliaenv/ src/compute/construct_transformations.jl $path_prefix
  echo "$step1_script completed."
}

total_steps=1
progress_bar_width=50

for ((step = 1; step <= total_steps; step++)); do
  # Execute a step
  "step${step}"

  # Update the progress bar
  print_progress_bar "$step" "$total_steps" "$progress_bar_width"
done

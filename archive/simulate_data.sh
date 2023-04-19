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
  -p)
    p="$2"
    shift # past argument
    shift # past value
    ;;
  -T)
    T="$2"
    shift # past argument
    shift # past value
    ;;
  -h_A)
    h_A="$2"
    shift # past argument
    shift # past value
    ;;
  -h_B)
    h_B="$2"
    shift # past argument
    shift # past value
    ;;
  -s1)
    s1="$2"
    shift # past argument
    shift # past value
    ;;
  -prefix)
    prefix="$2"
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
  # Run Julia script for step 1
  echo "Running $s1 ..."
  julia --project=juliaenv/ $s1 $p $T $h_A $h_B "${prefix}_T${T}_p${p}"
  echo "$s1 completed."
}

total_steps=1
progress_bar_width=50

for ((step = 1; step <= total_steps; step++)); do
  # Execute a step
  "step${step}"

  # Update the progress bar
  print_progress_bar "$step" "$total_steps" "$progress_bar_width"
done

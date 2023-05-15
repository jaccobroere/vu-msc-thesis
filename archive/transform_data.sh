#!/bin/bash

# Parse named arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
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
  # Check if step1_script argument was provided
  if [[ -z "$prefix" ]]; then
    echo "Error: step1_script argument is missing"
    exit 1
  fi

  # Run Julia script for step 1
  echo "Running precalculations_and_write.jl.jl ..."
  julia --project=$JULIA_DIR src/compute/jl/precalculations_and_write.jl.jl $prefix
  echo "precalculations_and_write.jl.jl completed."
}

step1
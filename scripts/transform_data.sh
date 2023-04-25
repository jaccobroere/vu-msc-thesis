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
  echo "Running transform_bootstrap_graph.jl ..."
  julia --project=juliaenv/ src/compute/transform_bootstrap_graph.jl $prefix
  echo "transform_bootstrap_graph.jl completed."
}

step1
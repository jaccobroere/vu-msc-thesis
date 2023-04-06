#!/bin/bash

# Define progress bar function
function progress {
    bar="                                                  "
    percent=$(($1 * 5))
    for i in $(seq 1 $percent); do
        bar=${bar:0:1}"="${bar:2}
    done
    echo -ne "\rProgress: [$bar]"
}

# Parse named arguments
while [[ $# -gt 0 ]]
do
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
        -path_prefix)
        path_prefix="$2"
        shift # past argument
        shift # past value
        ;;
        *)
        # unknown option
        echo "Unknown option: $1"
        exit 1
        ;;
    esac
done

# Run first Julia script
echo "Running simulation_design.jl..."
julia --project=/path/to/juliaenv/ src/simulation_design.jl $p $T $h_A $h_B $path_prefix
echo "simulation_design.jl completed."
progress 1

# Run second Julia script
echo "Running construct_V_simga.jl..."
julia --project=/path/to/juliaenv/ src/construct_V_simga.jl "out/${path_prefix}_y.csv" $path_prefix
echo "construct_V_simga.jl completed."
progress 2

# Run Python script
echo "Running construct_graph.py..."
python src/construct_graph.py $p
echo "construct_graph.py completed."
progress 3

# Run R script
echo "Running GSPLASH.R..."
Rscript src/GSPLASH.R $FOO $BAR
echo "GSPLASH.R completed."
progress 4

# Display progress bar
echo ""
echo "All scripts completed."
progress 5
echo ""

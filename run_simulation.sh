#!/bin/bash

# Define progress bar function
function progress {
    bar="                                                  "
    percent=$(($1 * 6))
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
        -gamma)
        gamma="$2"
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
julia --project=/path/to/juliaenv/ src/construct_V_sigma.jl $path_prefix
echo "construct_V_simga.jl completed."
progress 2

# Install PyJulia
echo "Installing PyJulia..."
python -m pip install julia
python src/install_pyjulia.py
echo "PyJulia installed."
progress 3

# Run Python script
echo "Running construct_graph.py..."
python src/construct_graph.py $p $path_prefix
echo "construct_graph.py completed."
progress 4

# Run R script
echo "Running GSPLASH.R..."
"C:\Program Files\R\R-4.2.1\bin\Rscript.exe" src/GSPLASH.R $path_prefix $gamma
echo "GSPLASH.R completed."
progress 5

# Display progress bar
echo ""
echo "All scripts completed."
progress 6
echo ""

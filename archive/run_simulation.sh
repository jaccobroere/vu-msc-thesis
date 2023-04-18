#!/bin/bash

# Define progress bar function
print_progress_bar() {
    local current_step=$1
    local total_steps=$2
    local width=$3

    # Calculate the percentage of completion
    local percentage=$((current_step * 100 / total_steps))

    # Calculate the number of completed and pending characters
    local completed_chars=$((width * current_step / total_steps))
    local pending_chars=$((width - completed_chars))

    # Create the progress bar string
    local progress_bar=""
    for ((i = 0; i < completed_chars; i++)); do
        progress_bar+="#"
    done
    for ((i = 0; i < pending_chars; i++)); do
        progress_bar+="-"
    done

    # Print the progress bar
    printf "\r[%s] %d%%" "$progress_bar" "$percentage"
    echo ""
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
        -fit)
        fit="$2"
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
step1() {
    echo "Running simulation_design.jl..."
    julia --project=/path/to/juliaenv/ src/simulation/simulation_design1.jl $p $T $h_A $h_B $path_prefix
    echo "simulation_design.jl completed."
}

step2() {
    # Run second Julia script
    echo "Running construct_V_simga.jl..."
    julia --project=/path/to/juliaenv/ src/compute/construct_V_sigma.jl $path_prefix
    echo "construct_V_simga.jl completed."
}

step3() {
    # Install PyJulia
    echo "Installing PyJulia..."
    python -m pip install julia igraph > /dev/null 2>&1
    python src/environment/install_pyjulia.py > /dev/null 2>&1
    echo "PyJulia installed."
}

step4() {
    # Run Python script
    echo "Running construct_graph.py..."
    python src/compute/construct_graph.py $p $path_prefix 
    echo "construct_graph.py completed."
}

step5() {
    # Run R script
    if [ "$fit" = true ] ; then
        echo "Running GSPLASH.R..."
        Rscript src/compute/GSPLASH.R $path_prefix 0.5 # > /dev/null 2>&1
        echo "GSPLASH.R completed."
    fi
}

total_steps=5
progress_bar_width=50

for ((step = 1; step <= total_steps; step++)); do
    # Execute a step
    "step${step}"

    # Update the progress bar
    print_progress_bar "$step" "$total_steps" "$progress_bar_width"
done



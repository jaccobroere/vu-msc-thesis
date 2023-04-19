# Progress bar
total_steps=2
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

# DESIGN A
p=49
T=500
h_A=3
h_B=7
s1=src/simulation/simulation_designA.jl
prefix=designA

# DESIGN B
p=25 # m^2 
T=500
path=src/simulation/simulation_designB.jl
prefix=designB


step_simA() {
    echo "Running $s1 ..."
    julia --project=juliaenv/ $s1 $p $T $h_A $h_B ${prefix}_T${T}_p${p}
    echo "$s1 completed."
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

step_simB() {
    echo "Running $path ..."
    julia --project=juliaenv/ $path $p $T ${prefix}_T${T}_p${p}
    echo "$path completed."
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

# Transform data
step_transform () {
    sh transform_data.sh -prefix ${prefix}_T${T}_p${p}
    current_step=$((current_step+1))
    print_progress_bar $current_step $total_steps 50
}

step_simB && step_transform
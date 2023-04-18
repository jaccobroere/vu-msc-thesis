p=10
T=500
h_A=2
h_B=2
s1=src/simulation/simulation_designA.jl
prefix=designA

sh simulate_data.sh -p $p -T $T -h_A $h_A -h_B $h_B -s1 $s1 -prefix $prefix
sh transform_data.sh -prefix ${prefix}_T${T}_p${p}
p=10
T=500
h_A=2
h_B=2
s1=designA

sh simulate_data.sh -p $p -T $T -h_A $h_A -h_B $h_B -s1 $s1
sh transform_data.sh -path_prefix "${s1}_T${T}_p${p}"
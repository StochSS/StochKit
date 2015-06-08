#!/bin/bash

cur_dir=$(pwd -P)

cd ..

output_options_array=("" "--keep-trajectories" "--keep-histograms")
interval_options_array=("" "-i 2")
models_array=("dimer_decay.xml" "heat_shock_mass_action.xml" "heat_shock_mixed.xml" "heat_shock_10x.xml" "schlogl.xml" "simple1p.xml" "stochkit_6_reaction_mixed.xml")
output_dir_array=("dimer_decay_output_case6-" "heat_shock_mass_action_output_case6-" "heat_shock_mixed_output_case6-" "heat_shock_10x_output_case6-" "schlogl_output_case6-" "simple1p_output_case6-" "stochkit_6_reaction_mixed_output_case6-")
seed_array=("524191051" "131598147" "1653256886" "1199494889" "1916221463" "1326376419" "100769862" "1154547651" "2000792645" "264702417" "189096838" "1677278027" "333957030" "1018928030" "1094973684" "1502831469" "193420206" "484280693" "1861941661" "1603636255" "426295216" "673464053" "1102473092" "1652131190" "972904451" "307168991" "1190100383" "1414477157" "1021671061" "475386278" "1041061928" "1693960895" "1454031231" "1305296388" "535798649" "58303751" "892943146" "1289944631" "1558231276" "2075758310" "173536695" "781867283")

count=0

for (( model_index = 0 ; model_index < ${#models_array[@]} ; model_index++ ))
do
    model_item=${models_array[$model_index]}
    output_dir_item=${output_dir_array[$model_index]}
    for (( output_index = 0 ; output_index < ${#output_options_array[@]} ; output_index++ ))
    do
       output_item=${output_options_array[$output_index]}
       for (( interval_index = 0 ; interval_index < ${#interval_options_array[@]} ; interval_index++ ))
       do
          interval_item=${interval_options_array[$interval_index]}
          printf -v count_pad "%02d" $count
          echo "case$count_pad" "$model_item" "$output_item" "$interval_item"
          if [ ! -d "$cur_dir/logs_tau/case$count_pad" ]
          then
             mkdir -m 700 "$cur_dir/logs_tau/case$count_pad"
          fi
          seed_item=${seed_array[$count]}
          ./tau_leaping -m "models/examples/$model_item" -t 10 -r 10 $interval_item --out-dir "models/examples/$output_dir_item$count_pad" "$output_item" -f --seed "$seed_item" -p 7 >"$cur_dir/logs_tau/case$count_pad/stdout3.txt" 2>&1
          let count=count+1
       done
    done
done

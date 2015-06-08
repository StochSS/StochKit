#!/bin/bash

cur_dir=$(pwd -P)

cd ..

output_options_array=("" "--keep-trajectories" "--keep-histograms")
interval_options_array=("" "-i 2")
models_array=("dimer_decay.xml" "heat_shock_mass_action.xml" "heat_shock_mixed.xml" "heat_shock_10x.xml" "schlogl.xml" "simple1p.xml" "stochkit_6_reaction_mixed.xml")
output_dir_array=("dimer_decay_output_case5-" "heat_shock_mass_action_output_case5-" "heat_shock_mixed_output_case5-" "heat_shock_10x_output_case5-" "schlogl_output_case5-" "simple1p_output_case5-" "stochkit_6_reaction_mixed_output_case5-")
seed_array=("1730909305" "1885376660" "1435457865" "1984852867" "1859755714" "1540707823" "1339599395" "248533812" "978973177" "943644076" "1998122333" "400779576" "416395489" "1110929690" "365988646" "108586539" "400913392" "1044242911" "880050353" "2117660381" "1810259562" "572851195" "186949734" "1935406109" "1167029832" "1788674320" "281055317" "751176089" "1484157281" "338606989" "1565750124" "146995529" "2043182840" "1026564428" "976634346" "2034650099" "64582304" "269688486" "1504660474" "750493361" "103041090" "370207993")

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
          ./tau_leaping -m "models/examples/$model_item" -t 10 -r 10 $interval_item --out-dir "models/examples/$output_dir_item$count_pad" "$output_item" -f --seed "$seed_item" >"$cur_dir/logs_tau/case$count_pad/stdout2.txt" 2>&1
          let count=count+1
       done
    done
done

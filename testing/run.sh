#!/bin/bash

cur_dir=$(pwd -P)

cd ..

output_options_array=("" "--keep-trajectories" "--keep-histograms")
interval_options_array=("" "-i 2")
models_array=("dimer_decay.xml" "heat_shock_mass_action.xml" "heat_shock_mixed.xml" "heat_shock_10x.xml" "schlogl.xml" "simple1p.xml" "stochkit_6_reaction_mixed.xml" "events1.xml" "stochkit_events.xml")
output_dir_array=("dimer_decay_output_case2-" "heat_shock_mass_action_output_case2-" "heat_shock_mixed_output_case2-" "heat_shock_10x_output_case2-" "schlogl_output_case2-" "simple1p_output_case2-" "stochkit_6_reaction_mixed_output_case2-" "events1_output_case2-" "stochkit_events_output_case2-")
seed_array=("75052381" "1398073729" "5680007" "1077743024" "2118936082" "664369650" "2055419909" "760559220" "1023622860" "2036910381" "1642909878" "1431171429" "1306538207" "2030777300" "888739240" "314621099" "61595042" "1815770255" "1865208005" "1975984741" "508256436" "1253109497" "622203109" "489785400" "1102076933" "1879248251" "1546058063" "503837780" "949504794" "997401045" "2050241623" "15124994" "2105646627" "1400994527" "1451888834" "1083425801" "1096413365" "1946545249" "1024568742" "826919755" "176600655" "660700426" "1756686300" "1028088496" "1559174225" "307935544" "817866424" "1525961881" "1411120843" "1614754526" "1954198377" "316433292" "1968413273" "1411934090")

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
          if [ ! -d "$cur_dir/logs/case$count_pad" ]
          then
             mkdir -m 700 "$cur_dir/logs/case$count_pad"
          fi
          seed_item=${seed_array[$count]}
          ./ssa -m "models/examples/$model_item" -t 10 -r 10 $interval_item --out-dir "models/examples/$output_dir_item$count_pad" "$output_item" -f --seed "$seed_item" >"$cur_dir/logs/case$count_pad/stdout2.txt" 2>&1
          let count=count+1
       done
    done
done

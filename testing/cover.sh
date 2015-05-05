#!/bin/bash

pushd ".." > /dev/null
par_dir=$(pwd -P)
popd > /dev/null

cp -p --update "$par_dir/src/"*.gcda "$par_dir/obj/"
cp -p --update "$par_dir/src/"*.gcno "$par_dir/obj/"
 
gcov --object-directory "$par_dir/obj" -p --source-prefix "$par_dir/src" --relative-only "$par_dir/src/solvers/BinHeap.cpp" "$par_dir/src/utility/Random.cpp" "$par_dir/src/model_parser/StringCalculator.cpp" "$par_dir/src/solvers/LDMTree.cpp" "$par_dir/src/drivers/CommandPass.cpp" "$par_dir/src/utility/CommandLineInterface.cpp" "$par_dir/src/drivers/MPISimulation.cpp" "$par_dir/src/drivers/CommandPassAux.cpp" "$par_dir/src/drivers/ParallelIntervalSimulation.cpp" "$par_dir/src/utility/AutomaticSelection.cpp" "$par_dir/src/utility/Parameter.cpp" "$par_dir/src/utility/StandardDriverUtilities.cpp" "$par_dir/src/drivers/tau_leaping_exp_adapt_serial.cpp" "$par_dir/src/drivers/tau_leaping_exp_adapt_mixed_serial.cpp" "$par_dir/src/drivers/ssa_serial.cpp" "$par_dir/src/utility/TimeBasedTrigger.cpp" "$par_dir/src/drivers/ssa.cpp" "$par_dir/src/drivers/tau_leaping.cpp"  "$par_dir/src/solvers/ConstantTimeGroup.cpp" "$par_dir/src/solvers/ConstantTimeGroupCollection.cpp" > stats.txt 2>&1
mv *.gcov coverage/

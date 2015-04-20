1. require libxml2 installed in system
2. for pure mass action input:
    make DM_MASS_ACTION
    ./test_mass_action_SSA_Direct $input_filename
3. for mixed input:
    make DM_MIXED
    ./test_mixed_SSA_Direct $input_filename
   once you have done this once, you can use
    ./test_mixed_SSA_Direct_recompile $input_filename
   to call the compiled version to save compilation time
4. stochkit_mass_action.xml is an example input file of pure mass action input,
   stochkit_mixed.xml is an exmaple input file of mixed input

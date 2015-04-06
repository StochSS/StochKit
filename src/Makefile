include ../.make/makefile.c
include ../.make/dependency.h
include ../.make/src4obj.h

stochkit_ssa: 	ssa \
		ssa_serial 
#		ssa_compiled

stochkit_tau :  tau_leaping \
		tau_leaping_exp_adapt_serial
#		tau_leaping_exp_adapt \
#		tau_leaping_exp_adapt_mixed \
#		tau_leaping_exp_adapt_mixed_compiled

stochkit_utility: 	precompiled_header \
			random \
			parameter \
			commandline \
			automatic_selection \
			input \
			std_driver \
			parallel \
			mpi \
                        commandexec \
                        commandpass \
                        commandpassaux \
			commandpassaux_serial

random: $(RANDOM_DEPS)
	$(CXX) $(RANDOM_SRC) -c -o $(STOCHKIT_OBJ)/Random.o $(COMPILE_ONLY_OPTIONS)

parameter: $(PARAMETER)
	$(CXX) $(PARAMETER_SRC) -c -o $(STOCHKIT_OBJ)/Parameter.o $(COMPILE_ONLY_OPTIONS)

commandline: $(COMMAND_LINE_INTERFACE)
	$(CXX) $(COMMAND_LINE_INTERFACE_SRC) -c -o $(STOCHKIT_OBJ)/CommandLineInterface.o $(COMPILE_ONLY_OPTIONS)

automatic_selection: $(AUTOMATIC_SELECTION)
	$(CXX) $(AUTOMATIC_SELECTION_SRC) -c -o $(STOCHKIT_OBJ)/AutomaticSelection.o $(COMPILE_ONLY_OPTIONS)

input: $(INPUT_DEPS)
	$(CXX) $(INPUT_SRC) -c -o $(STOCHKIT_OBJ)/StringCalculator.o $(COMPILE_ONLY_OPTIONS)

std_driver: $(STD_DRIVER_UT_DEPS)
	$(CXX) $(STD_DRIVER_UT_SRC) -c -o $(STOCHKIT_OBJ)/StandardDriverUtilities.o $(COMPILE_ONLY_OPTIONS)

parallel: $(PARALLEL_INTERVAL_DEPS)
	$(CXX) $(PARALLEL_INTERVAL_SRC) -c -o $(STOCHKIT_OBJ)/ParallelIntervalSimulation.o $(COMPILE_ONLY_OPTIONS)

mpi: $(MPI_DEPS)
ifeq ($(COMPILE_MPI),1)
	$(MPICXX) $(MPI_SRC) -c -o $(STOCHKIT_OBJ)/MPISimulation.o $(COMPILE_ONLY_OPTIONS) $(EXTRA_OPTIONS)
else
	$(CXX) $(MPI_SRC) -c -o $(STOCHKIT_OBJ)/MPISimulation.o $(COMPILE_ONLY_OPTIONS)
endif

commandexec: $(COMMAND_EXEC_DEPS)
	$(CXX) $(COMMAND_EXEC_SRC) -o $(STOCHKIT_BIN)/CommandExec $(ALLOPTIONS)

commandpass: $(COMMAND_PASS_DEPS)
	$(CXX) $(COMMAND_PASS_SRC) -c -o $(STOCHKIT_OBJ)/CommandPass.o $(COMPILE_ONLY_OPTIONS)

commandpassaux: $(COMMAND_PASS_AUX_DEPS)
ifeq ($(COMPILE_MPI),1)
	$(CXX) $(COMMAND_PASS_AUX_SRC) -c -o $(STOCHKIT_OBJ)/CommandPassAux.o $(COMPILE_ONLY_OPTIONS) $(EXTRA_OPTIONS)
else
	$(CXX) $(COMMAND_PASS_AUX_SRC) -c -o $(STOCHKIT_OBJ)/CommandPassAux.o $(COMPILE_ONLY_OPTIONS)
endif

commandpassaux_serial: $(COMMAND_PASS_AUX_DEPS)
	$(CXX) $(COMMAND_PASS_AUX_SRC) -c -o $(STOCHKIT_OBJ)/CommandPassAux_serial.o $(COMPILE_ONLY_OPTIONS)

precompiled_header: $(PRECOMPILED_HEADER_DEPS)
	$(CXX) $(PRECOMPILED_HEADER_SRC) $(COMPILE_ONLY_OPTIONS)

ldm_tree: $(LDMTREE_DEPS)
	$(CXX) $(LDM_TREE_SRC) -c -o $(STOCHKIT_OBJ)/LDMTree.o $(COMPILE_ONLY_OPTIONS)

constant_group: $(CONSTANT_GROUP_DEPS)
	$(CXX) $(CONSTANT_GROUP_SRC) -c -o $(STOCHKIT_OBJ)/ConstantTimeGroup.o $(COMPILE_ONLY_OPTIONS)

constant_group_collection: $(CONSTANT_GROUP_DEPS)
	$(CXX) $(CONSTANT_GROUP_COLLECTION_SRC) -c -o $(STOCHKIT_OBJ)/ConstantTimeGroupCollection.o $(COMPILE_ONLY_OPTIONS)

nrm_bin_heap: $(NRM_BIN_HEAP_DEPS)
	$(CXX) $(NRM_BIN_HEAP_SRC) -c -o $(STOCHKIT_OBJ)/BinHeap.o $(COMPILE_ONLY_OPTIONS)

ssa: constant_group constant_group_collection nrm_bin_heap ldm_tree $(SSA_DEPS)
ifeq ($(COMPILE_MPI),1)
	$(MPICXX) -DSTOCHKIT_SYSTEM_APP $(SSA_SRC) -o  $(STOCHKIT_BIN)/ssa $(ALLOPTIONS) $(EXTRA_OPTIONS)
else
	$(CXX) -DSTOCHKIT_SYSTEM_APP $(SSA_SRC) -o  $(STOCHKIT_BIN)/ssa $(ALLOPTIONS)
endif

ssa_serial: constant_group constant_group_collection nrm_bin_heap ldm_tree $(SSA_SERIAL_DEPS)
	$(CXX) -DSTOCHKIT_SYSTEM_APP $(SSA_SERIAL_SRC) -o  $(STOCHKIT_BIN)/ssa_serial $(ALLOPTIONS) 

ssa_compiled: constant_group constant_group_collection nrm_bin_heap ldm_tree $(SSA_COMPILED_DEPS)
	$(CXX) -DSTOCHKIT_SYSTEM_APP -DMIXED $(SSA_COMPILED_SRC) -o  "$(GENERATED_CODE_PATH)/bin/ssa_compiled" -include "$(GENERATED_CODE_PATH)/CustomPropensityFunctions.h" $(ALLOPTIONS)

#added for unified events-7/11/13 Arya
ssa_compiled_mixed: constant_group constant_group_collection nrm_bin_heap ldm_tree $(SSA_COMPILED_DEPS)
	$(CXX) -DEVENTS $(SSA_COMPILED_SRC) -o  "$(GENERATED_CODE_PATH)/bin/ssa_compiled_mixed" -include "$(GENERATED_CODE_PATH)/CustomPropensityFunctions.h" -include "$(GENERATED_CODE_PATH)/CustomStateBasedTriggerFunctions.h" -include "$(GENERATED_CODE_PATH)/CustomChangeSingleSpeciesFunctions.h" $(ALLOPTIONS)

tau_leaping: $(TAUL_DEPS)
ifeq ($(COMPILE_MPI),1)
	$(MPICXX) -DSTOCHKIT_SYSTEM_APP $(TAUL_SRC) -o  $(STOCHKIT_BIN)/tau_leaping $(ALLOPTIONS) $(EXTRA_OPTIONS)
else
	$(CXX) -DSTOCHKIT_SYSTEM_APP $(TAUL_SRC) -o  $(STOCHKIT_BIN)/tau_leaping $(ALLOPTIONS)
endif

tau_leaping_exp_adapt: $(TAUL_EXP_ADAPT_DEPS)
	$(CXX) -DSTOCHKIT_SYSTEM_APP $(TAUL_EXP_ADAPT_SRC) -o $(STOCHKIT_BIN)/tau_leaping_exp_adapt $(ALLOPTIONS)

tau_leaping_exp_adapt_serial: $(TLEAS_DEPS)
	$(CXX) -DSTOCHKIT_SYSTEM_APP -DMATRIX_STOICHIOMETRY $(TLEAS_SRC) -o $(STOCHKIT_BIN)/tau_leaping_exp_adapt_serial $(ALLOPTIONS)

tau_leaping_exp_adapt_mixed: $(TLEAM_DEPS)
	$(CXX) -DSTOCHKIT_SYSTEM_APP $(TLEAM_SRC) -o  $(STOCHKIT_BIN)/tau_leaping_exp_adapt_mixed $(ALLOPTIONS)

tau_leaping_exp_adapt_mixed_compiled: $(TEAMC_DEPS)
	$(CXX) -DSTOCHKIT_SYSTEM_APP -DMATRIX_STOICHIOMETRY -DMIXED $(TEAMC_SRC) -o "$(GENERATED_CODE_PATH)/bin/tau_leaping_exp_adapt_mixed_compiled" -include "$(GENERATED_CODE_PATH)/CustomPropensityFunctions.h" $(ALLOPTIONS)

clean:
	rm -f $(ODIR)/*.o core*   

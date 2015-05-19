#utility
_RANDOM_DEPS = Random.h Random.cpp
_PARAMETER = Parameter.h Parameter.cpp
_UTILITY_DEPS = StandardDriverUtilities.h StandardDriverUtilities.cpp CommandLineInterface.h CommandLineInterface.cpp StandardDriverTypes.h
_TRIGGER = TimeBasedTrigger.h TimeBasedTrigger.cpp
_COMMAND_LINE_INTERFACE = CommandLineInterface.h CommandLineInterface.cpp
_AUTOMATIC_SELECTION = AutomaticSelection.h AutomaticSelection.cpp
_STD_DRIVER_UT_DEPS = StandardDriverUtilities.h StandardDriverUtilities.cpp CommandLineInterface.h StandardDriverTypes.h
_PRECOMPILED_HEADER_DEPS = boost_headers.h

RANDOM_DEPS  =  $(patsubst %,$(STOCHKIT_SRC)/utility/%,$(_RANDOM_DEPS))
PARAMETER = $(patsubst %,$(STOCHKIT_SRC)/utility/%,$(_PARAMETER))
UTILITY_DEPS =  $(patsubst %,$(STOCHKIT_SRC)/utility/%,$(_UTILITY_DEPS))
TRIGGER      =  $(patsubst %,$(STOCHKIT_SRC)/utility/%,$(_TRIGGER))
COMMAND_LINE_INTERFACE  =  $(patsubst %,$(STOCHKIT_SRC)/utility/%,$(_COMMAND_LINE_INTERFACE))
AUTOMATIC_SELECTION  =  $(patsubst %,$(STOCHKIT_SRC)/utility/%,$(_AUTOMATIC_SELECTION))
STD_DRIVER_UT_DEPS  =  $(patsubst %,$(STOCHKIT_SRC)/utility/%,$(_STD_DRIVER_UT_DEPS))
PRECOMPILED_HEADER_DEPS  =  $(patsubst %,$(STOCHKIT_SRC)/utility/%,$(_PRECOMPILED_HEADER_DEPS))

#output
_OUTPUT_DEPS = StatsOutput.h  IntervalOutput.h StandardDriverOutput.h

OUTPUT_DEPS =  $(patsubst %,$(STOCHKIT_SRC)/output/%,$(_OUTPUT_DEPS))

#model_parser
_INPUT_DEPS = Input.h Input.ipp Input_mass_action.h StringCalculator.h StringCalculator.cpp 
_INPUT_MIXED_BEF_AFT_COM = Input_mixed_before_compile.h Input_mixed_after_compile.h
_INPUT_EVENTS_BEF_AFT_COM = Input_events_before_compile.h Input_events_after_compile.h
_INPUT_EVENTS_BEF_COM = Input_events_before_compile.h 
_INPUT_EVENTS_AFT_COM = Input_events_after_compile.h
_TAG_DEP = Input_tag.h ModelTag.h


INPUT_DEPS = $(patsubst %,$(STOCHKIT_SRC)/model_parser/%,$(_INPUT_DEPS))
INPUT_MIXED_BEF_AFT_COM = $(patsubst %,$(STOCHKIT_SRC)/model_parser/%,$(_INPUT_MIXED_BEF_AFT_COM))
INPUT_EVENTS_BEF_AFT_COM = $(patsubst %,$(STOCHKIT_SRC)/model_parser/%,$(_INPUT_EVENTS_BEF_AFT_COM))
INPUT_EVENTS_BEF_COM = $(patsubst %,$(STOCHKIT_SRC)/model_parser/%,$(_INPUT_EVENTS_BEF_COM))
INPUT_EVENTS_AFT_COM = $(patsubst %,$(STOCHKIT_SRC)/model_parser/%,$(_INPUT_EVENTS_AFT_COM))
TAG_DEP = $(patsubst %,$(STOCHKIT_SRC)/model_parser/%,$(_TAG_DEP))

#drivers
_SERIAL_DEPS = SerialIntervalSimulationDriver.h
_PARALLEL_INTERVAL_DEPS = ParallelIntervalSimulation.h ParallelIntervalSimulation.cpp
_MPI_DEPS = MPISimulation.h MPISimulation.cpp
_COMMAND_EXEC_DEPS = CommandExec.h CommandExec.cpp
_COMMAND_PASS_DEPS = CommandPass.h CommandPass.cpp
_COMMAND_PASS_AUX_DEPS = CommandPassAux.h CommandPassAux.cpp

SERIAL_DEPS =  $(patsubst %,$(STOCHKIT_SRC)/drivers/%,$(_SERIAL_DEPS))
PARALLEL_INTERVAL_DEPS =  $(patsubst %,$(STOCHKIT_SRC)/drivers/%,$(_PARALLEL_INTERVAL_DEPS))
MPI_DEPS = $(PARALLEL_INTERVAL_DEPS) $(patsubst %,$(STOCHKIT_SRC)/drivers/%,$(_MPI__DEPS))
COMMAND_EXEC_DEPS = $(patsubst %,$(STOCHKIT_SRC)/drivers/%,$(_COMMAND_EXEC_DEPS))
COMMAND_PASS_DEPS = $(patsubst %,$(STOCHKIT_SRC)/drivers/%,$(_COMMAND_PASS_DEPS))
COMMAND_PASS_AUX_DEPS = $(patsubst %,$(STOCHKIT_SRC)/drivers/%,$(_COMMAND_PASS_AUX_DEPS))

#solvers
#_SSA_DIRECT_DEPS = SSA_Base.h SSA_Base.ipp SSA_Direct.h SSA_Direct.ipp
_SSA_DIRECT_DEPS = SSA_Direct.h SSA_Direct.ipp

SSA_DIRECT_DEPS = $(patsubst %,$(STOCHKIT_SRC)/solvers/%,$(_SSA_DIRECT_DEPS))

_CONSTANT_TIME_HDEPS = SSA_ConstantTime.h SSA_ConstantTime.ipp
CONSTANT_TIME_HDEPS = $(patsubst %,$(STOCHKIT_SRC)/solvers/%,$(_CONSTANT_TIME_HDEPS))

_TAUL_EXP_ADP_DEPS = TauLeapingExplicitAdaptive.h TauLeapingExplicitAdaptive.ipp
TAUL_EXP_ADP_DEPS = $(patsubst %,$(STOCHKIT_SRC)/solvers/%,$(_TAUL_EXP_ADP_DEPS))

_CONSTANT_GROUP_DEPS = ConstantTimeGroup.cpp ConstantTimeGroupCollection.cpp ConstantTimeGroup.h ConstantTimeGroupCollection.h
CONSTANT_GROUP_DEPS = $(patsubst %,$(STOCHKIT_SRC)/solvers/%,$(_CONSTANT_GROUP_DEPS))

_NRM_BIN_HEAP_DEPS = BinHeap.cpp BinHeap.h
NRM_BIN_HEAP_DEPS = $(patsubst %,$(STOCHKIT_SRC)/solvers/%,$(_NRM_BIN_HEAP_DEPS))

_SSA_DEPS = $(_SSA_DIRECT_DEPS) $(_CONSTANT_TIME_HDEPS) \
			$(_CONSTANT_GROUP_DEPS) $(_NRM_BIN_HEAP_DEPS)
SSA_DEPS = $(SSA_DIRECT_DEPS) $(CONSTANT_TIME_HDEPS) \
			$(CONSTANT_GROUP_DEPS) $(NRM_BIN_HEAP_DEPS)

#for ssa
SSA_DEPS = $(PARALLEL_INTERVAL_DEPS) $(MPI_DEPS) $(PARAMETER) $(UTILITY_DEPS) $(COMMAND_LINE_INTERFACE) $(AUTOMATIC_SELECTION) $(TAG_DEP) $(STOCHKIT_SRC)/drivers/ssa.cpp

#for ssa_serial
SSA_SERIAL_DEPS= $(STOCHKIT_SRC)/drivers/ssa_serial.cpp $(OUTPUT_DEPS) $(SSA_DEPS) $(RANDOM_DEPS) $(INPUT_DEPS) $(SERIAL_DEPS) $(UTILITY_DEPS)

#for ssa_compiled
SSA_COMPILED_DEPS =  $(STOCHKIT_SRC)/drivers/ssa_serial.cpp $(OUTPUT_DEPS) $(SSA_DEPS) $(PARAMETER) $(RANDOM_DEPS) $(SERIAL_DEPS) $(INPUT_DEPS)  $(UTILITY_DEPS)

#for tau_leaping
TAUL_DEPS = $(STOCHKIT_SRC)/drivers/tau_leaping.cpp $(COMMAND_LINE_INTERFACE) $(AUTOMATIC_SELECTION) $(TAG_DEP) $(PARALLEL_INTERVAL_DEPS) $(MPI_DEPS) $(UTILITY_DEPS)

#for tau_leaping_exp_adapt
TAUL_EXP_ADAPT_DEPS = $(STOCHKIT_SRC)/drivers/tau_leaping_exp_adapt.cpp $(PARALLEL_INTERVAL_DEPS) $(UTILITY_DEPS)

#for tau_leaping_exp_adapt_serial
TLEAS_DEPS = $(STOCHKIT_SRC)/drivers/tau_leaping_exp_adapt_serial.cpp \
			$(OUTPUT_DEPS)  $(SSA_DIRECT_DEPS) $(RANDOM_DEPS) $(INPUT_DEPS) $(COMMAND_LINE_INTERFACE) \
			$(AUTOMATIC_SELECTION) $(TAUL_EXP_ADP_DEPS) $(PARALLEL_INTERVAL_DEPS)  $(SERIAL_DEPS) $(UTILITY_DEPS) 

#for tau_leaping_exp_adapt_mixed
TLEAM_DEPS = $(STOCHKIT_SRC)/drivers/tau_leaping_exp_adapt_mixed.cpp \
			$(OUTPUT_DEPS) $(TAUL_EXP_ADP_DEPS) $(SSA_DIRECT_DEPS) $(RANDOM_DEPS) $(INPUT_DEPS) $(COMMAND_LINE_INTERFACE) \
			$(AUTOMATIC_SELECTION) $(INPUT_MIXED_BEF_AFT_COM) $(UTILITY_DEPS) 

#for tau_leaping_exp_adapt_mixed_compiled
TEAM_DEPS = $(STOCHKIT_SRC)/drivers/tau_leaping_exp_adapt_mixed_compiled.cpp \
			$(OUTPUT_DEPS) $(TAUL_EXP_ADP_DEPS) $(SSA_DIRECT_DEPS) $(RANDOM_DEPS) $(INPUT_DEPS)  $(SERIAL_DEPS) $(COMMAND_LINE_INTERFACE) \
			$(AUTOMATIC_SELECTION) $(UTILITY_DEPS) 


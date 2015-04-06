include ../../.make/makefile.c
include ../../.make/dependency.h
include ../../.make/src4obj.h

all:	ssa \
		ssa_serial
#		ssa_compiled

ssa:
ifeq ($(COMPILE_MPI),1)
	$(MPICXX) -DUSER_DRIVER $(CUSTOM_SSA_SRC) -o ./bin/ssa -I. -include ./user_ensemble_simulate.h $(ALLOPTIONS) $(EXTRA_OPTIONS)
else
	$(CXX) -DUSER_DRIVER $(CUSTOM_SSA_SRC) -o  ./bin/ssa -I. -include ./user_ensemble_simulate.h $(ALLOPTIONS)
endif


ssa_serial:
	$(CXX) -DUSER_DRIVER $(CUSTOM_SSA_SERIAL_SRC) -o  ./bin/ssa_serial -I. -include ./user_ensemble_simulate.h $(ALLOPTIONS) 

ssa_compiled:
	$(CXX) -DUSER_DRIVER -DMIXED $(CUSTOM_SSA_COMPILED_SRC) -o  $(GENERATED_CODE_PATH)/bin/ssa_compiled -include $(GENERATED_CODE_PATH)/CustomPropensityFunctions.h -I. -include ./user_ensemble_simulate.h $(ALLOPTIONS)

clean:
	rm -f ./bin/*   

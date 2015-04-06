include ../../.make/makefile.c
include ../../.make/dependency.h
include ../../.make/src4obj.h

all:	ssa\
		ssa_serial
#		ssa_compiled

ssa:
ifeq ($(COMPILE_MPI),1)
	$(MPICXX) -DUSER_OUTPUT $(CUSTOM_SSA_SRC) -o ./bin/ssa -I. -include ./UserOutput.h $(ALLOPTIONS) $(EXTRA_OPTIONS)
else
	$(CXX) -DUSER_OUTPUT $(CUSTOM_SSA_SRC) -o  ./bin/ssa -I. -include ./UserOutput.h $(ALLOPTIONS)
endif


ssa_serial:
	$(CXX) -DUSER_OUTPUT $(CUSTOM_SSA_SERIAL_SRC) -o  ./bin/ssa_serial -I. -include ./UserOutput.h $(ALLOPTIONS) 

ssa_compiled:
	$(CXX) -DUSER_OUTPUT -DMIXED $(CUSTOM_SSA_COMPILED_SRC) -o  $(GENERATED_CODE_PATH)/bin/ssa_compiled -include $(GENERATED_CODE_PATH)/CustomPropensityFunctions.h -I. -include ./UserOutput.h $(ALLOPTIONS)

clean:
	rm -f ./bin/*   

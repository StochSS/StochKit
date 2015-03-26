#set environment variables needed by includes
STOCHKIT_HOME := ../../..
export STOCHKIT_HOME

#pull in variables from main StochKit makefiles
include ../../../.make/makefile.c
include ../../../.make/dependency.h
include ../../../.make/src4obj.h

constant_helper = ../../../obj/ConstantTimeGroup.o ../../../obj/ConstantTimeGroupCollection.o
nrm_helper =  obj/BinHeap.o

SOLVER_HELPERS = $(constant_helper) $(nrm_helper)

calibrateModel: staticCalibrator.o StaticCalibratorCommandLineInterface.o GetRunTime.o
	$(CXX) calibrateModel.cpp obj/staticCalibrator.o obj/StaticCalibratorCommandLineInterface.o obj/GetRunTime.o $(RANDOM_SRC) $(PARAMETER_SRC) $(INPUT_SRC) $(SOLVER_HELPERS) -o calibrateModel $(ALLOPTIONS)

#for non-mass action (AKA mixed) models
staticCalibratorCompiled: CustomPropensityFunctions.h StaticCalibratorCommandLineInterface.o ../GetRunTimeMixed.h ../GetRunTimeMixed.cpp
	$(CXX) -include CustomPropensityFunctions.h -DMIXED  $(RANDOM_SRC) $(PARAMETER_SRC) $(INPUT_SRC) $(SOLVER_HELPERS) obj/StaticCalibratorCommandLineInterface.o ../GetRunTimeMixed.cpp calibrateModelMixed.cpp -o calibrateModelMixed $(ALLOPTIONS)

staticCalibrator.o:  staticCalibrator.cpp StaticCalibratorCommandLineInterface.o BinHeap.o
	$(CXX) staticCalibrator.cpp -c -o obj/staticCalibrator.o $(COMPILE_ONLY_OPTIONS)

GetRunTime.o: ../GetRunTime.cpp ../GetRunTime.h
	$(CXX) ../GetRunTime.cpp -c -o obj/GetRunTime.o $(COMPILE_ONLY_OPTIONS)

StaticCalibratorCommandLineInterface.o : StaticCalibratorCommandLineInterface.cpp StaticCalibratorCommandLineInterface.h
	$(CXX) StaticCalibratorCommandLineInterface.cpp -c -o obj/StaticCalibratorCommandLineInterface.o $(COMPILE_ONLY_OPTIONS)

BinHeap.o : ../../solvers/BinHeap.cpp ../../solvers/BinHeap.h
	g++ ../../solvers/BinHeap.cpp -c -o obj/BinHeap.o

clean:
	rm -f obj/*.o calibrateModel
CC=g++
CFLAG=-c -g
BIN_=bin/
STLIB_PATH=../libs
BOOST_PATH=../libs/boost_1_41_0
BOOSTLIB_PATH=../libs/boost_1_41_0/stage/lib
BOOSTBIN_PATH=../bin/lib/
LIBXML_PATH=/usr/include/libxml2
LIBXML_LIB=xml2
INCLUDE_PATH=-I$(STLIB_PATH) -I$(BOOST_PATH) -I$(LIBXML_PATH) -I../src/
LINK_SEARCH_PATH=-L$(BOOSTBIN_PATH)
LINK_SEARCH_NAME=-l$(LIBXML_LIB)


#all: testrun.o
#	g++ bin/Data.o bin/Reaction.o bin/ReactionList.o bin/AdaptiveTauSSA.o bin/TauLeaping.o bin/SSA.o bin/BaseSystem.o bin/StoichiometryMatrix.o bin/StochVector.o testrun.o -o testrun 
#testrun.o: testrun.cpp
#	g++ -c testrun.cpp $(INCLUDE_PATH)
all: SSA_Direct \
	ODM

ODM: test_SSA_ODM.cpp TestModel.h ../src/SSA_ODM.h  ../src/Random.h ../src/Random.cpp 
	g++ -g -O3 -Wall -DVECTOR_OF_VECTORS_STOICHIOMETRY test_SSA_ODM.cpp ../src/Random.cpp ../src/StringCalculator.cpp -o ODM $(INCLUDE_PATH) $(INCLUDE_PATH) $(LINK_SEARCH_NAME)

NEWODM: test_ODM.cpp TestModel.h ../src/SSA_ODM.h  ../src/Random.h ../src/Random.cpp 
	g++ -g -O3 -Wall -DVECTOR_OF_VECTORS_STOICHIOMETRY test_ODM.cpp ../src/Random.cpp  -o new_ODM_test $(INCLUDE_PATH)  $(INCLUDE_PATH) $(LINK_SEARCH_NAME) 

SSA_Direct: ../src/SSA_Direct.h DimerDecay.h unit_SSA_Direct.cpp
	g++ -O3 -pg unit_SSA_Direct.cpp -o unit_SSA_Direct.o $(INCLUDE_PATH) $(LINK_SEARCH_PATH)

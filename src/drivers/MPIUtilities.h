#ifndef _MPI_UTILITIES_H_
#define _MPI_UTILITIES_H_

#include <string>
#include "CommandLineInterface.h"


namespace STOCHKIT
{

class MPIUtilities
{
public:
   static void MPIcompileMixed(std::string executableName, const CommandLineInterface& commandLine, const bool stochkit_system_app, const bool events=false);
   static void initialize(int * argc, char ***argv);
   static void finalize(void);
};

class PrintUtilities
{
private:
   static int rank;
   friend void MPIUtilities::initialize(int *argc, char ***argv);
public:
   static void printOnce(std::ostream &out, const char *const str);
};

}
#endif

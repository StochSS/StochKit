#include "CommandPassAux.h"

#ifdef COMPILE_MPI
#include "MPIUtilities.h"
#include "CommandPass.h"

static STOCHKIT::CommandPass cmdExecutor;

#endif

#include <cstdlib>

int cmdExec(const char *const cmd)
{
#ifdef COMPILE_MPI
   return cmdExecutor.execute(cmd);
#else
   return system(cmd);
#endif
}

void exitFunc(const int status)
{
#ifdef COMPILE_MPI
   STOCHKIT::MPIUtilities::finalize();
   exit(status);
#else
   exit(status);
#endif
}

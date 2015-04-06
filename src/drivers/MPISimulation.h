/******************************************************************************
 */

#ifndef _MPI_SIMULATION_H_
#define _MPI_SIMULATION_H_

#include "ParallelIntervalSimulation.h"
#ifndef WIN32
#include <unistd.h>
#endif

#ifdef COMPILE_MPI
#include <mpi.h>
#endif

#include "MPIUtilities.h"

namespace STOCHKIT
{

class MPISimulation : public ParallelIntervalSimulation
{
	
public:
  MPISimulation(int ac, char* av[]);

  void run(std::string executableName);
  
  void mergeOutput();
  
  ~MPISimulation();

protected:
  unsigned int processRank;
  bool active;
#ifdef COMPILE_MPI
  MPI_Comm communicator;
#else
  int communicator;
#endif
};

}

#endif

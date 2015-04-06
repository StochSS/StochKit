/**
  @file ssa.cpp

  @brief primary main for running Stochkit on multiple processors

  -Determines which solver to use
  -Determines the type of the model (mixed, mass action, events) and acts accordingly
  -Launches realizations of a simulation in parallel and creates output using ParallelIntervalSimulation.h
  -For mixed and events models StandardDriverUtilities::compileMixed() is called, creating an executable
  version of ssa_serial.cpp

  @see ParallelIntervalSimulation
  @see StandardDriverUtilities
*/

#include "boost_headers.h"

#include "ParallelIntervalSimulation.h"
#include "MPISimulation.h"
#include "CommandLineInterface.h"
#include "StandardDriverUtilities.h"
#include "CommandPassAux.h"
#include <cstdio>
#include <cstdlib>
#ifndef WIN32
#include <csignal>
#include <unistd.h>

#ifdef COMPILE_MPI
#include <mpi.h>
#endif

void signal_handler(int sig)
{
    if( sig == SIGTERM || sig == SIGHUP || sig == SIGABRT ){
        std::cout << "\nStochKit ERROR: signal received, exiting..." << std::endl;
        kill(0, sig);
        exit(1);
    }
}
#endif

using namespace STOCHKIT;

int main(int ac, char* av[])
{
#ifdef COMPILE_MPI
  MPIUtilities::initialize(&ac, &av);
#endif

#ifndef WIN32
  setpgid(0,0);

  if (signal(SIGTERM, signal_handler) == SIG_ERR){
    std::cout << "\nStochKit Error: can't catch SIGTERM." << std::endl;
    exitFunc(1);
  }

  if (signal(SIGHUP, signal_handler) == SIG_ERR){
    std::cout << "\nStochKit Error: can't catch SIGHUP." << std::endl;
    exitFunc(1);
  }

  if (signal(SIGABRT, signal_handler) == SIG_ERR){
    std::cout << "\nStochKit Error: can't catch SIGABRT." << std::endl;
    exitFunc(1);
  }
#endif

  CommandLineInterface commandLine(ac,av);

  PrintUtilities::printOnce(std::cout, "StochKit MESSAGE: determining appropriate solver...\n");

  int methodId;
  if( commandLine.shouldCalibrate() )
    methodId = 100;  
  else
    methodId = commandLine.getMethod();

  switch(methodId)
  {
     case 0:
       {
	PrintUtilities::printOnce(std::cout, "StochKit MESSAGE: Simulating using Direct solver...\n");
	break;
       }
     case 1:
       {
	PrintUtilities::printOnce(std::cout, "StochKit MESSAGE: Simulating using Direct solver...\n");
	break;
       }
     case 10:
       {
	PrintUtilities::printOnce(std::cout, "StochKit MESSAGE: Simulating using ODM solver...\n");
	break;
       }
     case 11:
       {
	PrintUtilities::printOnce(std::cout, "StochKit MESSAGE: Simulating using ODM solver...\n");
	break;
       }
     case 20:
       {
	PrintUtilities::printOnce(std::cout, "StochKit MESSAGE: Simulating using LDM solver...\n");
	break;
       }
     case 30:
       {
	PrintUtilities::printOnce(std::cout, "StochKit MESSAGE: Simulating using ConstantTime solver...\n");
	break;
       }
     case 40:
       {
	PrintUtilities::printOnce(std::cout, "StochKit MESSAGE: Simulating using NRM solver...\n");
	break;
       }
     case 100:
       {
	PrintUtilities::printOnce(std::cout, "StochKit MESSAGE: Running Calibrator...\n");
	break;
       }
     default:
       {
          PrintUtilities::printOnce(std::cout, "StochKit ERROR: Solver type not recognized. Possibly using a tau-leaping solver name.\n");
          return -1;
       }
  }

  //MASS ACTION
  if(commandLine.isMassAction()){
	boost::filesystem::path currentPath=boost::filesystem::path(av[0]).parent_path(); 
	std::string executablepath=boost::filesystem::system_complete(currentPath).string();
#ifdef WIN32
	executablepath+="\\ssa_serial.exe";

	//add quote for the executable
	executablepath="\""+executablepath+"\"";
#endif

#ifdef COMPILE_MPI
	MPISimulation parallelDriver(ac, av);
#else
	ParallelIntervalSimulation parallelDriver(ac, av);
#endif

  if(commandLine.shouldCalibrate())
  {
#ifdef WIN32
	parallelDriver.runCalibrator(executablepath);
#elif defined STOCHKIT_SYSTEM_APP
	parallelDriver.runCalibrator("\"$STOCHKIT_HOME/bin/ssa_serial\"");
#else
	executablepath+="/ssa_serial";
	parallelDriver.runCalibrator(executablepath);
#endif

  parallelDriver.mergeCalibratorOutput();
  }
  else
  {
#ifdef WIN32
	parallelDriver.run(executablepath);
#elif defined STOCHKIT_SYSTEM_APP
	parallelDriver.run("\"$STOCHKIT_HOME/bin/ssa_serial\"");
#else
	executablepath+="/ssa_serial";
	parallelDriver.run(executablepath);
#endif
	parallelDriver.mergeOutput();

   }

  }

  //MIXED
  else if(commandLine.isMixed()){
#ifdef WIN32
	boost::filesystem::path currentPath=boost::filesystem::path(av[0]).parent_path();
#endif

	std::string executableName="ssa_compiled";

#if defined(WIN32)
	StandardDriverUtilities::compileMixed(executableName,commandLine,currentPath);
#elif defined(COMPILE_MPI) && defined(STOCHKIT_SYSTEM_APP)
	MPIUtilities::MPIcompileMixed(executableName, commandLine, true);
#elif defined(COMPILE_MPI) && !defined(STOCHKIT_SYSTEM_APP)
	MPIUtilities::MPIcompileMixed(executableName, commandLine, false);
#elif !defined(COMPILE_MPI) && defined(STOCHKIT_SYSTEM_APP)
	StandardDriverUtilities::compileMixed(executableName,commandLine,true);
#elif !defined(COMPILE_MPI) && !defined(STOCHKIT_SYSTEM_APP)
	StandardDriverUtilities::compileMixed(executableName,commandLine,false);
#endif

#ifdef COMPILE_MPI
	MPISimulation parallelDriver(ac, av);
#else
	ParallelIntervalSimulation parallelDriver(ac, av);
#endif


  if(commandLine.shouldCalibrate())
  {
#ifdef WIN32
	parallelDriver.runCalibrator("\""+currentPath.string()+"\\"+executableName+"\"");
#else
	parallelDriver.runCalibrator("\""+commandLine.getGeneratedCodeDir()+"/bin/"+executableName+"\"");
#endif

	parallelDriver.mergeCalibratorOutput();
  }
  else
  {
#ifdef WIN32
	parallelDriver.run("\""+currentPath.string()+"\\"+executableName+"\"");
#else
	parallelDriver.run("\""+commandLine.getGeneratedCodeDir()+"/bin/"+executableName+"\"");
#endif

	parallelDriver.mergeOutput();
  }

   } 

   //EVENTS
   else if(commandLine.isEventsEnabled()){
#ifdef WIN32
  boost::filesystem::path currentPath=boost::filesystem::path(av[0]).parent_path();
#endif

  std::string executableName="ssa_compiled_mixed";

#if defined(WIN32)
	StandardDriverUtilities::compileMixed(executableName,commandLine,currentPath, true);
#elif defined(COMPILE_MPI) && defined(STOCHKIT_SYSTEM_APP)
	MPIUtilities::MPIcompileMixed(executableName, commandLine, true, true);
#elif defined(COMPILE_MPI) && !defined(STOCHKIT_SYSTEM_APP)
	MPIUtilities::MPIcompileMixed(executableName, commandLine, false, true);
#elif !defined(COMPILE_MPI) && defined(STOCHKIT_SYSTEM_APP)
	StandardDriverUtilities::compileMixed(executableName,commandLine,true, true);
#elif !defined(COMPILE_MPI) && !defined(STOCHKIT_SYSTEM_APP)
	StandardDriverUtilities::compileMixed(executableName,commandLine,false, true);
#endif


#ifdef COMPILE_MPI
	MPISimulation parallelDriver(ac, av);
#else
	ParallelIntervalSimulation parallelDriver(ac, av);
#endif

  if(commandLine.shouldCalibrate())
  {
#ifdef WIN32
  parallelDriver.runCalibrator("\""+currentPath.string()+"\\"+executableName+"\"");
#else
  parallelDriver.runCalibrator("\""+commandLine.getGeneratedCodeDir()+"/bin/"+executableName+"\"");
#endif
  parallelDriver.mergeCalibratorOutput();
  }
  else
  {
#ifdef WIN32
  parallelDriver.run("\""+currentPath.string()+"\\"+executableName+"\"");
#else
  parallelDriver.run("\""+commandLine.getGeneratedCodeDir()+"/bin/"+executableName+"\"");
#endif
  parallelDriver.mergeOutput();
  }

   }

   //ERROR
   else {
        PrintUtilities::printOnce(std::cerr, "StochKit ERROR: Model is other than mass-action, mixed, or events_enabled type. Please make sure there is no mistake in the model file.\n");
       return -1;
   }

#ifdef COMPILE_MPI
   MPIUtilities::finalize();
#endif
   return 0;
}

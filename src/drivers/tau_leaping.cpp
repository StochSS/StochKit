/*
*  tau_leaping.cpp
*  analyze model and choose which tau-leaping method to use
*/
#include "boost_headers.h"

#include "CommandLineInterface.h"
#include "Input_tag.h"
#include "ModelTag.h"
#include "boost/thread/thread.hpp"
#include "ParallelIntervalSimulation.h"//for unified purpose
#include "MPISimulation.h"
#include "CommandPassAux.h"
#ifndef WIN32
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <csignal>
//#include <unistd.h>

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

	{ //begin internal block

		//need to decide which solver (e.g. tau_leaping_exp_adapt or tau_leaping_exp_adapt_mixed)
		CommandLineInterface commandLine(ac,av);
#ifdef COMPILE_MPI
		MPISimulation parallelDriver(ac, av);
#else
		ParallelIntervalSimulation parallelDriver(ac, av);//for unified purpose
#endif

#ifdef WIN32 //it seems visual studio does not recognize the one line statement
		std::string filename=commandLine.getModelFileName();
        	char *const modelFileName = static_cast<char *>(_alloca(filename.length() + 2048) );
		memset(modelFileName, 0, filename.length() + 2048);
		strcpy(modelFileName, filename.c_str());
#else
        	char modelFileName[commandLine.getModelFileName().length() + 2048];
		memset(modelFileName, 0, commandLine.getModelFileName().length() + 2048);
		strcpy(modelFileName, commandLine.getModelFileName().c_str());
#endif

		Input_tag<ModelTag> input_model_tag(modelFileName);

		ModelTag model_tag = input_model_tag.writeModelTag();

		ModelTag::ModelType modelType = model_tag.Type;

#ifdef WIN32
		boost::filesystem::path currentPath=boost::filesystem::path(av[0]).parent_path();
#endif
		std::string executableName;//for unified purpose

		PrintUtilities::printOnce(std::cout, "StochKit MESSAGE: determining appropriate driver...");

		//if events, use event solver
		if (modelType==ModelTag::events_enabled) {
			PrintUtilities::printOnce(std::cout, "\nStochKit ERROR: events detected. tau_leaping is not event-enabled. Terminating\n");
			exitFunc(1);
		}
		else {

			if (modelType==ModelTag::mixed) {
#ifdef STATIC
				std::cout<<"\nStochKit ERROR: The Windows lite version of StochKit does not support custom propensity functions (try the full version). Terminating." << std::endl;
				return 0;
#else
				executableName="tau_leaping_exp_adapt_mixed_compiled";
				std::string printStr("running ");
				printStr += executableName;
				printStr += "...\n";
				PrintUtilities::printOnce(std::cout, printStr.c_str() );
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

#ifdef WIN32
				parallelDriver.run("\""+currentPath.string()+"\\"+executableName+"\"");
#else
				parallelDriver.run("\""+commandLine.getGeneratedCodeDir()+"/bin/"+executableName+"\"");
#endif
#endif
			}
			else
			{
				executableName="tau_leaping_exp_adapt_serial";
				std::string printStr("running ");
				printStr += executableName;
				printStr += "...\n";
				PrintUtilities::printOnce(std::cout, printStr.c_str() );
#ifdef WIN32
				parallelDriver.run("\""+currentPath.string()+"\\"+executableName+"\"");
#else
				parallelDriver.run("\"$STOCHKIT_HOME/bin/tau_leaping_exp_adapt_serial\"");
#endif
			}
			parallelDriver.mergeOutput();
		}
	} //end internal block
#ifdef COMPILE_MPI
   MPIUtilities::finalize();
#endif
	return 0;
}

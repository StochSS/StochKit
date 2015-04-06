#include "MPISimulation.h"
#include "CommandPassAux.h"

namespace STOCHKIT
{

#ifdef COMPILE_MPI

static inline void MPICheckError(const int retCode, const char *const msg)
{
   if(retCode != MPI_SUCCESS)
   {
      std::cerr << "MPI error: " << msg << std::endl;
      exit(1);
   }
}

static inline void syncHelper(void)
{
   sync();
   sync();
   sync();
}

class MPISerialize
{
private:
   char *ptr;
   const MPI_Comm comm;
   const int totalRank;
   static const int tag1 = 0, tag2 = 1;
public:
   inline MPISerialize(const MPI_Comm com, const int totalRank_);

   inline void enter(const int rank);

   inline void exit(const int rank);

   inline ~MPISerialize();
};

inline MPISerialize::MPISerialize(const MPI_Comm com, const int totalRank_) : comm(com), totalRank(totalRank_)
{
   ptr = new char[4];
   memset(ptr, 0, 4);
}

inline void MPISerialize::enter(const int rank)
{
   MPICheckError(MPI_Barrier(comm), "MPI_Barrier error");
   if(rank > 0)
   {
      MPICheckError(MPI_Ssend(ptr, 4, MPI_CHAR, 0, tag1, comm), "MPI_Ssend error");
   }
}

inline void MPISerialize::exit(const int rank)
{
   MPI_Status stat;
   int i;
   if(rank == 0)
   {
      for(i = 1; i < totalRank; i++)
      {
         MPICheckError(MPI_Recv(ptr, 4, MPI_CHAR, i, tag1, comm, &stat), "MPI_Recv error");
         MPICheckError(MPI_Recv(ptr, 4, MPI_CHAR, i, tag2, comm, &stat), "MPI_Recv error");
      }
   }
   else
   {
      MPICheckError(MPI_Ssend(ptr, 4, MPI_CHAR, 0, tag2, comm), "MPI_Ssend error");
   }
   MPICheckError(MPI_Barrier(comm), "MPI_Barrier error");
}

inline MPISerialize::~MPISerialize()
{
   delete [] ptr;
   ptr = NULL;
}


void MPIUtilities::MPIcompileMixed(std::string executableName, const CommandLineInterface& commandLine, const bool stochkit_system_app, const bool events)
{
   int processRank;
   syncHelper(); //flush
   MPICheckError(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier error");
   MPICheckError(MPI_Comm_rank(MPI_COMM_WORLD, &processRank), "MPI_Comm_rank error");
   if(0 == processRank) {
	StandardDriverUtilities::compileMixed(executableName,commandLine, stochkit_system_app, events);
   }
   syncHelper(); //flush
   MPICheckError(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier error");
   //repeated
   syncHelper(); //flush
   MPICheckError(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier error");
}

void MPIUtilities::initialize(int * argc, char ***argv)
{
   MPICheckError(MPI_Init(argc, argv), "MPI_Init error");
   MPICheckError(MPI_Comm_rank(MPI_COMM_WORLD, &PrintUtilities::rank), "MPI_Comm_rank error");
}

void MPIUtilities::finalize(void)
{
   MPICheckError(MPI_Finalize(), "MPI_Finalize error");
}

void PrintUtilities::printOnce(std::ostream &out, const char *const str)
{
   if(rank < 0)
   {
      std::cerr << "Not properly initialized " << std::endl;
      exit(1);
   }
   else if(rank == 0)
   {
      out << str << std::flush;
   }
}

int PrintUtilities::rank = -1;

#else
void MPIUtilities::MPIcompileMixed(std::string executableName, const CommandLineInterface& commandLine, const bool stochkit_system_app, const bool events)
{
   (void)executableName;
   (void)commandLine;
   (void)stochkit_system_app;
   (void)events;
   std::cerr << "Not compiled with MPI" << std::endl;
   exit(1);
}

void MPIUtilities::initialize(int * argc, char ***argv)
{
   (void)argc;
   (void)argv;
   std::cerr << "Not compiled with MPI" << std::endl;
   exit(1);
}

void MPIUtilities::finalize(void)
{
   std::cerr << "Not compiled with MPI" << std::endl;
   exit(1);
}

void PrintUtilities::printOnce(std::ostream &out, const char *const str)
{
   out << str << std::flush;
}

int PrintUtilities::rank = -1;

#endif

#ifdef COMPILE_MPI
	MPISimulation::MPISimulation(int ac, char* av[]) : ParallelIntervalSimulation(ac, av, 1)
	{
                numberOfProcesses=commandLine.getProcesses();
                int processSize, processId;
                unsigned int numprocs;
                MPICheckError(MPI_Comm_size(MPI_COMM_WORLD, &processSize), "MPI_Comm_size error");
                MPICheckError(MPI_Comm_rank(MPI_COMM_WORLD, &processId), "MPI_Comm_rank error");

                numprocs = processSize;
                processRank = processId;
		if (numberOfProcesses > numprocs) {
			if(0 == processRank) {
				std::cout << "StochKit MESSAGE (MPISimulation()): Requested number of processes exceeds number of processors." << std::endl;
			}
			numberOfProcesses=numprocs;
		}

		//if default number of processes, 0, is chosen, automatically determine
		if (numberOfProcesses==0) {
			numberOfProcesses=numprocs;
			if (numberOfProcesses==0) {
				if(0 == processRank) {
					std::cout << "StochKit MESSAGE (MPISimulation()): unable to detect number of processors.  Simulation will run on one processor." << std::endl;
				}
				numberOfProcesses=1;
			}
		}

                if(processRank < numberOfProcesses) {
                   active = true;
                }
                else {
                   active = false;
                }
                MPICheckError(MPI_Comm_split(MPI_COMM_WORLD, ( (active == true) ? 0 : 1), ( (active == true) ? processRank : (processRank - numberOfProcesses) ), &communicator), "MPI_Comm_split error");
                MPICheckError(MPI_Comm_rank(communicator, &processId), "MPI_Comm_rank error");
                processRank = processId;

		if (commandLine.getUseSeed()) {
			boost::mt19937 newEngine(commandLine.getSeed() );
			engine=newEngine;
		}
		//give a warning if they're about to generate a lot of data
                if(0 == processRank && active == true)
                {
			warnIfLargeOutput();
                }
		std::cout << std::flush;
	}

        MPISimulation::~MPISimulation() {
           MPICheckError(MPI_Comm_free(&communicator), "MPI_Comm_free error");
        }

	void MPISimulation::run(std::string executableName) {
                if(active == false) {
                   return;
                }
                int retCode;
		int newseed;  
		syncHelper();  //flush
		MPICheckError(MPI_Barrier(communicator), "MPI_Barrier error");
		if(0 == processRank) {
			boost::uniform_int<> distribution(1, RAND_MAX) ;
			boost::variate_generator<boost::mt19937, boost::uniform_int<> > seedOfNewThread(engine, distribution);
			int otherseed;

			StandardDriverUtilities::createOutputDirs(commandLine,true,numberOfProcesses);
			std::cout << "StochKit MESSAGE: created output directory \""+commandLine.getOutputDir()+"\"..." << std::endl;
			syncHelper();  //flush to disk

			newseed = seedOfNewThread();
			for (std::size_t i=1; i<numberOfProcesses; i++) {
				otherseed = seedOfNewThread();
				MPICheckError(MPI_Ssend(&otherseed, 1, MPI_INT, static_cast<int>(i), 1, communicator), "MPI_Ssend error");
			}
		}
		else {
			MPI_Status status;
 			MPICheckError(MPI_Recv(&newseed, 1, MPI_INT, 0, 1, communicator, &status), "MPI_Recv error");
		}
		syncHelper();  //flush
		MPICheckError(MPI_Barrier(communicator), "MPI_Barrier error");

		//record full command in a log file
		std::string command;

		//need to convert size_t to strings
		std::string threadNumString;
		std::size_t subRealizations;
		std::string subRealizationsString;
		std::string seedString;
		std::string assignmentCounterString;

		std::size_t numRuns=commandLine.getRealizations();
		if(0 == processRank) {
			std::cout << "running simulation...";
		}
		std::cout.flush();


		timeval timer;
		double startTime, endTime;
		gettimeofday(&timer,NULL);
		startTime=timer.tv_sec+(timer.tv_usec/1000000.0);
		
		std::vector<std::size_t> assignmentCounterVector;
		std::size_t assignmentCounter=0;

		for (std::size_t i=0; i<numberOfProcesses; i++) {
			subRealizations = assignment(numRuns, i);
			numRuns -= subRealizations;

			processorRealizationAssignments.push_back(subRealizations);
			if (i==masterProc && processorRealizationAssignments[i]==0) {
				++masterProc;
			}

			if (subRealizations>0) {
				//if keeping trajectories, need to supply a trajectory offset
				if (commandLine.getKeepTrajectories()) {
					assignmentCounterVector.push_back(assignmentCounter);
					//increment assignmentCounter
					assignmentCounter+=subRealizations;
				}
			}
		}


		seedString=StandardDriverUtilities::size_t2string(newseed);
		subRealizations = processorRealizationAssignments.at(processRank);
                command = "";
		std::string commandStr="";
		threadNumString=StandardDriverUtilities::size_t2string(processRank);

		if (subRealizations>0) {
			subRealizationsString=StandardDriverUtilities::size_t2string(subRealizations);

			command=executableName+commandLine.getCmdArgs()+" --serial --use-existing-output-dirs --histograms-dir histograms/.parallel/thread"+threadNumString+" --stats-dir stats/.parallel --means-file means"+threadNumString+".txt --variances-file variances"+threadNumString+".txt --stats-info-file .stats-info-"+threadNumString+".txt";

			//append seed to command
			if (!commandLine.getUseSeed()) {
				command+=" --seed "+seedString;
			}
			else {
				//seed already exists in the command line arguments, so we have to replace it with seedString
				command=ParallelIntervalSimulation::modifyCmdArgsSeed(command,seedString);
			}

			command=ParallelIntervalSimulation::modifyCmdArgsRealizations(command,subRealizationsString);

			//if keeping trajectories, need to supply a trajectory offset
			if (commandLine.getKeepTrajectories()) {
				assignmentCounter = assignmentCounterVector.at(processRank);
				assignmentCounterString=StandardDriverUtilities::size_t2string(assignmentCounter);
				command+=" --trajectories-offset "+assignmentCounterString;
			}

			//if keeping histograms, thread 0 needs to write a histogram info file
			if (commandLine.getKeepHistograms() && processRank==masterProc) {
				command+=" --histograms-info-file ../.histogram-info.txt ";//should create it in the histograms/.parallel directory
			}

			//redirect messages from executable to log file

			command+=" > "+commandLine.getOutputDir()+"/.StochKit/thread-log"+threadNumString+".txt 2>&1";


			commandStr="echo 'THREAD "+threadNumString+":"+command+"' >> "+commandLine.getOutputDir()+"/.StochKit/command-log.txt"; 

		}
		else {
			commandStr="echo 'THREAD "+threadNumString+":"+command+"' >> "+commandLine.getOutputDir()+"/.StochKit/command-log.txt"; 
		}


		MPISerialize serialObj(communicator, numberOfProcesses);

		serialObj.enter(processRank);
                retCode = cmdExec(commandStr.c_str()); //mutual exclusion
                if(retCode != 0)
                {
		   std::cerr << "System call error code: " << retCode << " on " << processRank << std::endl;
                   std::cerr << "Command log string: " << commandStr.c_str() << std::endl;
                }
		syncHelper(); //flush 
		serialObj.exit(processRank);

		//execute serial code
		if (subRealizations>0) {
			sthreads.create_thread(boost::bind(&ParallelIntervalSimulation::executable, this, command));
		}

		//wait for all threads to finish
		sthreads.join_all();
                syncHelper(); //flush
		MPICheckError(MPI_Barrier(communicator), "MPI_Barrier error");
                //repeated
                syncHelper(); //flush
                MPICheckError(MPI_Barrier(communicator), "MPI_Barrier error");
		double elapsedTime;

		gettimeofday(&timer,NULL);
		endTime=timer.tv_sec+(timer.tv_usec/1000000.0);
		elapsedTime=endTime-startTime;

		if(0 == processRank) {
			std::cout << "finished (simulation time approx. " << elapsedTime << " second";
			if (elapsedTime!=1) {
				std::cout << "s";
			}
			std::cout << ")"<<std::endl;
		}
	}

	void MPISimulation::mergeOutput() {
		if(active == false) {
			return;
		}
		syncHelper(); //flush
		MPICheckError(MPI_Barrier(communicator), "MPI_Barrier error");
                //repeated
                syncHelper(); //flush
                MPICheckError(MPI_Barrier(communicator), "MPI_Barrier error");
		if(0 == processRank) {
			ParallelIntervalSimulation::mergeOutput();
		}
		syncHelper(); //flush
		MPICheckError(MPI_Barrier(communicator), "MPI_Barrier error");
	}


#else

	MPISimulation::MPISimulation(int ac, char* av[]) : ParallelIntervalSimulation(ac, av, 1)
	{
		std::cerr << "Not compiled with MPI" << std::endl;
		exit(1);
	}

        MPISimulation::~MPISimulation() {
		std::cerr << "Not compiled with MPI" << std::endl;
		exit(1);
        }

        void MPISimulation::run(std::string executableName)
	{
		(void)executableName;
		std::cerr << "Not compiled with MPI" << std::endl;
		exit(1);
	}

	void MPISimulation::mergeOutput()
	{
		std::cerr << "Not compiled with MPI" << std::endl;
		exit(1);
	}
#endif
}

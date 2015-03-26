/******************************************************************************
*/
#include "StaticCalibratorCommandLineInterface.h"

namespace STOCHKIT
{
  StaticCalibratorCommandLineInterface::StaticCalibratorCommandLineInterface(int ac, char* av[]):
    visible("command line options")
  {
    std::string temp_dir;
    visible.add_options()
      ("model,m", boost::program_options::value<std::string>(&modelFileName),"**REQUIRED Model file name")
      ("time,t",boost::program_options::value<double>(&simulationTime),"**REQUIRED Simulation time (i.e. run each realization from t=0 to t=time)")
      ("replications,r",boost::program_options::value<std::size_t>(&replications),"Number of Replications")
      ("force,f","Overwrite existing output directory and output files without confirmation.")
      ("no-recompile","Use previously compiled model code. Applicable only to models with customized propensity functions.  See documentation for details.")
      ("seed",boost::program_options::value<int>(&seed),"Seed the random number generator")
      
      ("help,h","Use -h or --help to list all arguments");
    
    //HIDDEN OPTIONS
    hidden.add_options()
      ("use-existing-output-dirs","[hidden] used by parallel driver to tell a sub-process that the output directories already exist")
      ("stats-dir",boost::program_options::value<std::string>(&statsDir)->default_value("stats"),"[hidden]")
      ("means-file",boost::program_options::value<std::string>(&meansFileName)->default_value("means.txt"),"[hidden]")
      ("variances-file",boost::program_options::value<std::string>(&variancesFileName)->default_value("variances.txt"),"[hidden]")
      ("stats-info-file",boost::program_options::value<std::string>(&statsInfoFileName)->default_value(".stats-info.txt"),"[hidden]")
      ("trajectories-dir",boost::program_options::value<std::string>(&trajectoriesDir)->default_value("trajectories"),"[hidden]")
      ("trajectories-offset",boost::program_options::value<std::size_t>(&trajectoriesOffset)->default_value(0),"[hidden]")
      ("histograms-dir",boost::program_options::value<std::string>(&histogramsDir)->default_value("histograms"),"[hidden]")
      ("histograms-info-file",boost::program_options::value<std::string>(&histogramsInfoFileName)->default_value(""),"[hidden]")//only write an info file if this is non-empty (nonempty only in parallel simulation)
      ("ssa-steps",boost::program_options::value<std::size_t>(&SSASteps)->default_value(100),"[hidden] Number of SSA steps to take before trying another tau leap step (applicable to tau leaping only).");
    
    combined.add(visible).add(hidden);
    
    parse(ac, av);
    
    //save command-line arguments
    //skip the first (av[0] is the executable name)
    *av++;
    ac--;
    //pull off the arguments
    while (ac--) {
      cmdArgs+=" "+(std::string)*av++;
    }
    
#ifdef WIN32
	//find the model file path and add quote for window version
	//for "-m"
	std::string serchingString="-m "+modelFileName;
	std::string replaceString="-m "+(std::string)"\""+modelFileName+(std::string)"\"";
	int index=cmdArgs.find(serchingString, 0);
	if(index!=std::string::npos)
		cmdArgs.replace(index, serchingString.size(), replaceString); 
	//for "--model"
	serchingString="--model "+modelFileName;
	replaceString="-m "+(std::string)"\""+modelFileName+(std::string)"\"";
	index=cmdArgs.find(serchingString, 0);
	if(index!=std::string::npos)
		cmdArgs.replace(index, serchingString.size(), replaceString);
	//for "--out-dir"
	serchingString="--out-dir "+temp_dir;
	replaceString="--out-dir "+(std::string)"\""+temp_dir+(std::string)"\"";
	index=cmdArgs.find(serchingString, 0);
	if(index!=std::string::npos)
		cmdArgs.replace(index, serchingString.size(), replaceString); 
#endif
}


std::string StaticCalibratorCommandLineInterface::getModelFileName() const {
	return modelFileName;
}

std::string StaticCalibratorCommandLineInterface::getGeneratedCodeDir() const {
	return generatedCodeDir;
}

double StaticCalibratorCommandLineInterface::getSimulationTime() const {
	return simulationTime;
}

std::size_t StaticCalibratorCommandLineInterface::getReplications() const {
	return replications;
}

bool StaticCalibratorCommandLineInterface::getUseSeed() const {
	return useSeed;
}

int StaticCalibratorCommandLineInterface::getSeed() const {
	return seed;
}

bool StaticCalibratorCommandLineInterface::getForce() const {
	return force;
}

bool StaticCalibratorCommandLineInterface::getRecompile() const {
	return recompile;
}

std::string StaticCalibratorCommandLineInterface::getCmdArgs() const {
	return cmdArgs;
}

  

//HIDDEH OPTION GETTERS
bool StaticCalibratorCommandLineInterface::getUseExistingOutputDirs() const {
	return useExistingOutputDirs;
}

std::string StaticCalibratorCommandLineInterface::getStatsDir() const {
	return statsDir;
}
std::string StaticCalibratorCommandLineInterface::getMeansFileName() const {
	return meansFileName;
}
std::string StaticCalibratorCommandLineInterface::getVariancesFileName() const {
	return variancesFileName;
}

std::string StaticCalibratorCommandLineInterface::getStatsInfoFileName() const {
	return statsInfoFileName;
}

std::string StaticCalibratorCommandLineInterface::getHistogramsInfoFileName() const {
	return histogramsInfoFileName;
}

std::string StaticCalibratorCommandLineInterface::getTrajectoriesDir() const {
	return trajectoriesDir;
}
std::size_t StaticCalibratorCommandLineInterface::getTrajectoriesOffset() const {
	return trajectoriesOffset;
}
std::string StaticCalibratorCommandLineInterface::getHistogramsDir() const {
	return histogramsDir;
}

std::size_t StaticCalibratorCommandLineInterface::getSSASteps() const {
	return SSASteps;
}
//END HIDDEN OPTION GETTERS



void StaticCalibratorCommandLineInterface::parse(int ac, char* av[]) {
	try {
		boost::program_options::store(boost::program_options::parse_command_line(ac,av,combined), vm);
	}
	catch (...) {
		std::cout << "StochKit ERROR (CommandLineInterface::parse): unable to parse command-line arguments.  Run with --help for a list of required and optional parameters.\n";
		exit(1);
	}

	boost::program_options::notify(vm);

	if (vm.count("help")) {
		std::cout << visible;
		exit(0);
	}

	if (!(vm.count("model") && vm.count("time") )) {
		std::cout << "StochKit ERROR (CommandLineInterface::parse): missing required parameter(s).  Run with --help for a list of required and optional parameters.\n";
		exit(1);
	}

	//if replications aren't specified set to default 1
	if ( !(vm.count("replications")) ) {
	  replications = 1;
	}

	if (replications==0) {
		std::cout << "StochKit ERROR (CommandLineInterface::parse): number of replications = 0. Exiting.\n";
		exit(1);
	}

	if (vm.count("seed")) {
		useSeed=true;
	}
	else {
		useSeed=false;
	}

	if (vm.count("force")) {
		force=true;
	}
	else {
		force=false;
	}

	//create full path to output directory, default location, <full_path_to/model_filename (without extension)>_output
	std::string tmp_full_model_path=boost::filesystem::system_complete(modelFileName).string();
	outputDir=tmp_full_model_path.substr(0,tmp_full_model_path.find_last_of("."))+"_output";

	if (vm.count("no-recompile")) {
		recompile=false;
	}
	else {
		recompile=true;
	}

	std::string full_model_path=boost::filesystem::system_complete(modelFileName).string();

	
	generatedCodeDir=full_model_path.substr(0,full_model_path.find_last_of("."))+"_generated_code";
	
	if (vm.count("use-existing-output-dirs")) {
		useExistingOutputDirs=true;
	}
	else {
		useExistingOutputDirs=false;
	}

}//end parse

}

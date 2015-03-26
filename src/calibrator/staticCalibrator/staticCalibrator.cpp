#include "staticCalibrator.h"

bool resultCmpr(result a, result b)
{
  return (a.averageTimePerStep < b.averageTimePerStep);
}

staticCalibrator::staticCalibrator(int argc, char* argv[])
{
  commandLineInterface = new STOCHKIT::StaticCalibratorCommandLineInterface(argc,argv);

  modelFileName = commandLineInterface->getModelFileName();
  simulationTime = commandLineInterface->getSimulationTime();
  nReplications = commandLineInterface->getReplications();
  
  //create input tag and model tag to determine if model has events or mixed propensities
  Input_tag<ModelTag> input_model_tag( (char*)modelFileName.c_str() );
  ModelTag model_tag = input_model_tag.writeModelTag();
  ModelTag::ModelType modelType = model_tag.Type;

  //calibrator cannot currently handle events
  if(modelType==ModelTag::events_enabled)
    {
      std::cout << "StochKit ERROR (GetRunTime): calibrator does not support events. Terminating.\n";
      exit(1);
    }

  //if model type is mixed we must go through a separate compilation process before running trajectories
  if(modelType==ModelTag::mixed)
    performSimulationsMixed();
  else
    performSimulations();
}

staticCalibrator::~staticCalibrator()
{
  delete commandLineInterface;
}

void staticCalibrator::performSimulations()
{
  //get total running times for each method over all replications
  //remember getRunTime functions returns average time to fire a reaction
  double odmTime=0.0, nrmTime=0.0, constTime=0.0;
  for(int i=0; i<nReplications; i++)
    {
      odmTime += getODMRunTime((char*)modelFileName.c_str(), simulationTime);
      nrmTime += getNRMRunTime((char*)modelFileName.c_str(), simulationTime);
      constTime += getConstantRunTime((char*)modelFileName.c_str(), simulationTime);
    }
  std::cout << "Simulations Complete!" << std::endl;
  
  //average total running times for each method over the number
  //of replications to allow for random noise in timing data
  results.resize(3);
  results[0].ssaName="ODM";   results[0].averageTimePerStep=odmTime / nReplications;
  results[1].ssaName="NRM";   results[1].averageTimePerStep=nrmTime / nReplications;
  results[2].ssaName="CONST"; results[2].averageTimePerStep=constTime / nReplications;
  std::sort(results.begin(), results.end(), resultCmpr);
  
  //set class member optimalMethod
  setResultString();
  std::cout << "Optimal Method Found!" << std::endl;
      
  //save verdict to file
  createModelCalibrationDirectory();
  createTimeCalibrationFile();
}

void staticCalibrator::performSimulationsMixed()
{
  std::cout << "StochKit MESSAGE (calibrator): Model contains customized propensities, compiling (this will take a few moments)...\n";

  Input_mixed_before_compile<StandardDriverTypes::populationType,
			     StandardDriverTypes::stoichiometryType,
			     StandardDriverTypes::propensitiesType,
			     StandardDriverTypes::graphType>
    model( (char*)modelFileName.c_str() );

  //create necessary header of custom propensity functions
  std::string customPropensityFunctionsFileName="CustomPropensityFunctions.h";
  model.writeCustomPropensityFunctionsFile( (char*)customPropensityFunctionsFileName.c_str() );

  //compile new custom calibrator
  std::string command="make staticCalibratorCompiled --silent > mixed-compile-log.txt 2>&1";
  int compilationSuccessful=system(command.c_str());

  //run the compiled calibrator if compilation was successful
  if(compilationSuccessful != 0)
    {
      std::cout << "StochKit ERROR (calibrator): error compiling custom propensity functions.  Check mixed-compile-log.txt for details.\n";
      exit(1);
    }
  else
    {
      command="./calibrateModelMixed"+commandLineInterface->getCmdArgs();
      std::cout << "calling executable: "<<command<<"\n";
      system(command.c_str());
    }
}

void staticCalibrator::createModelCalibrationDirectory()
{
  modelCalibrationDirectory = modelFileName;

  //truncate ".xml" from model file then append to get output directory
  std::size_t dotXMLSubstringLocation = modelCalibrationDirectory.find_last_of(".");
  modelCalibrationDirectory = modelCalibrationDirectory.substr(0,dotXMLSubstringLocation);
  modelCalibrationDirectory += "_CalibrationInfo";

  //if CalibrationInfo directory has not been created before then create
  try{
    if ( !(boost::filesystem::exists(modelCalibrationDirectory)) ) 
      boost::filesystem::create_directories(modelCalibrationDirectory);
  }
  catch (...) {
    std::cerr << "StochKit ERROR (staticCalibrator::createModelCalibrationDirectory): error creating output directory.\n";
    exit(1);
  }

}

void staticCalibrator::createTimeCalibrationFile()
{
  //convert simulation time to string as file name
  std::string calibrationFileName = modelCalibrationDirectory;
  calibrationFileName += "/calib";
  std::ostringstream oss;
  oss << simulationTime;
  calibrationFileName += oss.str();
  calibrationFileName += ".txt";

  try{
    //check to see if calibration has already been done on this model
    if (boost::filesystem::exists(calibrationFileName))
      {
	//check to see if the user want to re-calibrate despite
	//already having calibrated the model in the past
	if (!commandLineInterface->getForce())
	  {
	    std::cout << "StochKit ERROR (staticCalibrator::createTimeCalibrationFile): calibration info already exists for this model and simulation time.\n";
	    std::cout << "Run with --force or -f to overwrite.\n";
	    std::cout << "Calibration terminated.\n";
	    exit(1);
	}
	else
	  {
	    //delete existing file
	    if( remove(calibrationFileName.c_str()) != 0 )
	      std::cerr << "StochKit ERROR (staticCalibrator::createTimeCalibrationFile): could not remove old calibration file\n";
	  }
      }

    //create file
    std::ofstream calibrationFile;
    calibrationFile.open( calibrationFileName.c_str() );
    
    if ( !calibrationFile )
      {
	std::cerr << "StochKit ERROR (staticCalibrator::createTimeCalibrationFile): Unable to open output file. Terminating.\n";
	exit(1);
      }

    calibrationFile << resultString;
    calibrationFile.close();
    std::cout << "Calibration Data Saved!" << std::endl;
  }
  catch (...) {
    std::cerr << "StochKit ERROR (staticCalibrator::createTimeCalibrationFile): error creating output directory.\n";
    exit(1);
  }
    
}

void staticCalibrator::setResultString()
{
  std::stringstream ret;
  ret << results[0].ssaName << "\t" << results[0].averageTimePerStep << std::endl;
  ret << results[1].ssaName << "\t" << results[1].averageTimePerStep << std::endl;
  ret << results[2].ssaName << "\t" << results[2].averageTimePerStep << std::endl;

  resultString = ret.str();
}

/*
std::string findFastestMethod(double odm, double nrm, double constant);

int main(int argc, char* argv[])
{
  STOCHKIT::StaticCalibratorCommandLineInterface commandLineInterface(argc,argv);
  
  //retrieve applicable data from commandLineInterface
  char* modelFileName = (char*)(commandLineInterface.getModelFileName().c_str());//cast as char* since c_str returns const char*
  double simulationTime = commandLineInterface.getSimulationTime();
  std::size_t nReplications = commandLineInterface.getReplications();
  

  //get total running times for each method over all replications
  double odmTime=0.0, nrmTime=0.0, constTime=0.0;
  for(int i=0; i<nReplications; i++)
    {
      odmTime += getODMRunTime(modelFileName, simulationTime);
      nrmTime += getNRMRunTime(modelFileName, simulationTime);
      constTime += getConstantRunTime(modelFileName, simulationTime);
    }

  //average total running times for each method over the number
  //of replications to allow for random noise in timing data
  double odmAverageTime = odmTime / nReplications;
  double nrmAverageTime = nrmTime / nReplications;
  double constAverageTime = constTime / nReplications;

  //determine optimal method
  std::string optimalMethod = findFastestMethod(odmAverageTime, nrmAverageTime, constAverageTime);


  //save optimal method to same directory as model
  //

  return 0;
}


std::string findFastestMethod(double odm, double nrm, double constant)
{
  if( nrm > odm )
    {
      if( nrm > constant ) return "NRM";
      else return "CONST";
    }
  else
    {
      if( odm > constant ) return "ODM";
      else return "CONST";
    }
}

*/

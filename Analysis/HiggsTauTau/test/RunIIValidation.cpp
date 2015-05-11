#include <iostream>
#include <vector>
#include <string>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/bind.hpp"
#include "boost/function.hpp"
#include "boost/format.hpp"
#include "TSystem.h"
// #include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnRootTools.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/AnalysisBase.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/CopyCollection.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/SimpleFilter.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/OverlapFilter.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/RunIIElectronValidation.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/RunIIMuonValidation.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/RunIITauValidation.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/RunIIJetValidation.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/RunIIGenJetValidation.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/RunIIReclusteredGenJetValidation.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/RunIIReclusteredJetValidation.h"

using boost::lexical_cast;
using boost::bind;
namespace po = boost::program_options;
using std::string;
using std::vector;
using namespace ic;

int main(int argc, char* argv[]){

  // Configurable parameters
  int max_events;                 // Maximum number of events to process
  string filelist;                // The file containing a list of files to use as input
  string input_prefix;            // A prefix that will be added to the path of each input file
  string output_name;             // Name of the ouput ROOT File
  string output_folder;           // Folder to write the output in
	bool from_miniaod;
	string skim_path="";
  
  po::options_description config("Configuration");
  po::variables_map vm;
  po::notify(vm);

  config.add_options()    
      ("max_events",          po::value<int>(&max_events)-> default_value(-1))
      ("filelist",            po::value<string>(&filelist)->required())
      ("input_prefix",        po::value<string>(&input_prefix)->default_value(""))
      ("output_name",         po::value<string>(&output_name)->required())
      ("output_folder",       po::value<string>(&output_folder)->default_value(""))
      ("from_miniaod",        po::value<bool>(&from_miniaod)->default_value(false))
			
    ;
    po::store(po::command_line_parser(argc, argv).
              options(config).allow_unregistered().run(), vm);
    po::notify(vm);
  
  
  std::cout << "-------------------------------------" << std::endl;
  std::cout << "RunIIValidationStudy" << std::endl;
  std::cout << "-------------------------------------" << std::endl;      string param_fmt = "%-25s %-40s\n";
  
 
  // Load necessary libraries for ROOT I/O of custom classes
  // gSystem->Load("libFWCoreFWLite.dylib");
  // gSystem->Load("libUserCodeICHiggsTauTau.dylib");
  // AutoLibraryLoader::enable();

  // Build a vector of input files
  vector<string> files = ParseFileLines(filelist);
  for (unsigned i = 0; i < files.size(); ++i) files[i] = input_prefix + files[i];
  
  // Create ROOT output fileservice
  fwlite::TFileService *fs = new fwlite::TFileService((output_folder+output_name).c_str());
  
  //Optional configurable parameters with which you can filter the collection you are interested in.
 //Defining an analysis using this class removes all the need for any for loops through the events. An "analysis" loops through all the events up to max_events and runs the modules which you define at the end of this script using analysis.AddModule
  
  ic::AnalysisBase analysis(
    "RunIIValidationStudy",// Analysis name
    files,                // Input files
    "icEventProducer/EventTree", // TTree name
    max_events);          // Max. events to process (-1 = all)
  analysis.SetTTreeCaching(true);
  analysis.StopOnFileFailure(true);
  analysis.RetryFileAfterFailure(7, 3);

  

 RunIIElectronValidation runIIElectronValidation = RunIIElectronValidation("runIIElectronValidation")
  .set_write_plots(false)
  .set_write_tree(true)
  .set_fs(fs);
	
RunIIMuonValidation runIIMuonValidation = RunIIMuonValidation("runIIMuonValidation")
  .set_write_plots(false)
  .set_write_tree(true)
  .set_fs(fs);

RunIITauValidation runIITauValidation = RunIITauValidation("runIITauValidation")
  .set_write_plots(false)
  .set_write_tree(true)
  .set_fs(fs);

RunIIJetValidation runIIJetValidation = RunIIJetValidation("runIIJetValidation")
  .set_write_plots(false)
  .set_write_tree(true)
  .set_fs(fs);

RunIIGenJetValidation runIIGenJetValidation = RunIIGenJetValidation("runIIGenJetValidation")
  .set_write_plots(false)
  .set_write_tree(true)
  .set_fs(fs);

RunIIReclusteredJetValidation runIIReclusteredJetValidation = RunIIReclusteredJetValidation("runIIReclusteredJetValidation")
  .set_write_plots(false)
  .set_write_tree(true)
  .set_fs(fs);

RunIIReclusteredGenJetValidation runIIReclusteredGenJetValidation = RunIIReclusteredGenJetValidation("runIIReclusteredGenJetValidation")
  .set_write_plots(false)
  .set_write_tree(true)
  .set_fs(fs);


 analysis.AddModule(&runIIElectronValidation);
 analysis.AddModule(&runIIMuonValidation);
 analysis.AddModule(&runIITauValidation);
 analysis.AddModule(&runIIJetValidation);
 analysis.AddModule(&runIIGenJetValidation);
if(from_miniaod){
 analysis.AddModule(&runIIReclusteredJetValidation);
 analysis.AddModule(&runIIReclusteredGenJetValidation);
}
 


  analysis.RunAnalysis();
  delete fs;
  return 0;
}


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/bind.hpp"
#include "boost/function.hpp"
#include "boost/format.hpp"
#include "boost/algorithm/string.hpp"
#include "TSystem.h"
// #include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnRootTools.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/AnalysisBase.h"

#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/CopyCollection.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/SimpleFilter.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/ComplexFilter.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/OverlapFilter.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/IDOverlapFilter.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/CompositeProducer.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/OneCollCompositeProducer.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/PileupWeight.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/EffectiveAreaIsolationFilter.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/MetSelection.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/MTSelection.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/MakeRunStats.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/JetEnergyCorrections.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/LumiMask.h"

#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/SingleMetMaker.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvDataTriggerFilter.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvWeights.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvWDecay.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/WDecaySeparator.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvZDecay.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/ModifyMet.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/JetMETModifier.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/JECStudy.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/MetEventFilters.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvPrint.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/CJVFilter.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/LightTree.h"
// #include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/LightTreeAM.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/LightTreeRDM.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/LightTreeTrig.h"

#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvConfig.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/SetDoW.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/SetDoDY.h"


using boost::lexical_cast;
using boost::bind;
namespace po = boost::program_options;
using std::string;
using std::vector;
using namespace ic;
int main(int argc, char* argv[]){

  // Configurable parameters
  string cfg;                     // The configuration file
  int max_events;                 // Maximum number of events to process
  string filelist;                // The file containing a list of files to use as input
  string input_prefix;            // A prefix that will be added to the path of each input file
  string output_name;             // Name of the ouput ROOT File
  string output_folder;           // Folder to write the output in
  bool doincludehighptz;          // For including high pt z sample

  string era_str;                 // Analysis data-taking era
  string mc_str;                  // Analysis MC production
  string prod;                    // Our production string

  string wstream;                 // W stream: enu, munu or taunu, or nunu for everything
  bool doICHEP2016;               // true = ICHEP dataset only
  bool doLastJSON;                // true = take last JSON 
  bool dojetEtaCut;               // true = cut eta region between 3 and 3.2
  bool is_data;                   // true = data, false = mc
  bool dojessyst;                 // Do Jet Energy Scale Systematic Run
  bool dodatajessyst;             // Do Alternate Data Jet Energy Scale Method Systematic Run
  bool jesupordown;               // If doing Jet Energy Scale Systematic Run, run with up or down correction (true for up, false for down)
  bool dosmear;                   // Do Smearing
  bool doType1Cor;                // Do type1 corrections (for single e from tautau group)
  bool doaltmatch;                // Do runmetuncertainties gen jet matching
  bool doetsmear;                 // Do runmetuncertainties smearing
  bool dogaus;                    // Do gaussian smearing for jets with no gen jet match
  bool dospring10gaus;            // Do gaussian smearing for jets with no gen jet match
  bool dojersyst;                 // Do Jet Energy Resolution Systematic Run
  bool jerbetterorworse;          // If doing Jet Energy Resolution Systematic Run, run with with better or worse (true for better, false for worse)
  bool douessyst;                 // Do Unclustered MET Systematic Run
  bool uesupordown;               // If doing Unclustered MET Systematic Run, run with up or down correction (true for up, false for down) 
  bool docrosschecktau;           // If doing cross check tau use alternate tau id discriminant
  bool do2015tauid;               // If doing cross check tau use alternate tau id discriminant
  bool do2016tauid;               // If doing cross check tau use alternate tau id discriminant
  bool taulepdiscrtight;          // Use tight electron and muon discriminants
  bool dojerdebug;                // Access runmetunc collections for debugging
  bool dotopreweighting;          // Do Top reweighting
  bool dopromptskim;              // Use prompt compatible light tree skimming
  bool donoskim;                  // Do no skimming

  string mettype;                 // MET input collection to be used
  string jettype;                 // JET input collection to be used
  string jesuncfile;              // File to get JES uncertainties from
  bool reapplyJEC;                // Reapply JEC from txt files
  string jecdata;                 // File list to get JEC to reapply offline

  bool doMetFilters;              // apply cleaning MET filters.
  string filters;
  //unsigned signal_region;       // DeltaPhi cut > 2.7
  bool dotrgeff;                  // Do trigger efficiency corrections
  bool do3dtrgeff;                // Do 3d trigger efficiency corrections
  bool do1dparkedtrgeff;          // Do 1d parked trigger efficiency corrections
  bool dofitted1dparkedtrgeff;    // Do 1d fitted parked trigger efficiency corrections
  bool dobinnedin2d1dtrgeff;      // Do 2d binned fitted 1d parked trigger efficiency corrections
  bool doidisoeff;                // Do lepton ID-iso efficiency corrections
  bool doidisoerr;                // Do lepton ID-iso efficiency correction error
  bool doidisoerrupordown;        // Do lepton ID-iso efficiency correction error up or down
  bool doidisoerrmuore;           // Do lepton ID-iso efficiency correction error for muons or electrons
  bool dolumixsweight;            // Do lumi*xs/evt weight online
  string inputparams;             // Params file to use for info on lumi xs and evt
  string trg_weight_file;         // Trigger weights to use
  string trg_to_use;              // Which trigger to use for data

  //CUTS
  double jet1ptcut;               //jet1ptcut
  double jet2ptcut;               //jet2ptcut
  double jetptprecut;             //jetptprecut, should be the CJV threshold
  double mjjcut;
  double detajjcut;

  unsigned debug;
  int randomseed;

  bool printEventList;  //print run,lumi,evt of events selected
  bool printEventContent; //print event content of events selected

  bool turnoffpuid;

  bool useOldLT;
  bool doTrigLT;
  bool doAllPairs;

 // Load the config
  po::options_description preconfig("Pre-Configuration");
  preconfig.add_options()("cfg", po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    ("max_events",            po::value<int>(&max_events)-> default_value(-1))
    ("filelist",              po::value<string>(&filelist)->required())
    ("input_prefix",          po::value<string>(&input_prefix)->default_value(""))
    ("output_name",           po::value<string>(&output_name)->required())
    ("output_folder",         po::value<string>(&output_folder)->default_value(""))
    ("doincludehighptz",      po::value<bool>(&doincludehighptz)->default_value(false))
    ("era",                   po::value<string>(&era_str)->required())
    ("mc",                    po::value<string>(&mc_str)->required())
    ("prod",                  po::value<string>(&prod)->required())
    ("wstream",               po::value<string>(&wstream)->default_value("nunu"))
    ("doICHEP2016",           po::value<bool>(&doICHEP2016)->default_value(false))
    ("doLastJSON",            po::value<bool>(&doLastJSON)->default_value(false))
    ("dojetEtaCut",           po::value<bool>(&dojetEtaCut)->default_value(false))
    ("is_data",               po::value<bool>(&is_data)->required())
    ("mettype",               po::value<string>(&mettype)->default_value("pfMetType1"))
    ("jet1ptcut",             po::value<double>(&jet1ptcut)->default_value(80.))
    ("jet2ptcut",             po::value<double>(&jet2ptcut)->default_value(70.))
    ("jetptprecut",           po::value<double>(&jetptprecut)->default_value(15.))
    ("mjjcut",                po::value<double>(&mjjcut)->default_value(800.))
    ("detajjcut",             po::value<double>(&detajjcut)->default_value(3.6))
    ("doMetFilters",          po::value<bool>(&doMetFilters)->default_value(false))
    ("filters",               po::value<string> (&filters)->default_value(""))
    ("dojessyst",             po::value<bool>(&dojessyst)->default_value(false))
    ("dodatajessyst",         po::value<bool>(&dodatajessyst)->default_value(false))
    ("jesupordown",           po::value<bool>(&jesupordown)->default_value(true))
    ("dojersyst",             po::value<bool>(&dojersyst)->default_value(false))
    ("jerbetterorworse",      po::value<bool>(&jerbetterorworse)->default_value(true))
    ("douessyst",             po::value<bool>(&douessyst)->default_value(false))
    ("uesupordown",           po::value<bool>(&uesupordown)->default_value(true))
    ("docrosschecktau",       po::value<bool>(&docrosschecktau)->default_value(false))
    ("do2015tauid",           po::value<bool>(&do2015tauid)->default_value(false))
    ("do2016tauid",           po::value<bool>(&do2016tauid)->default_value(true))
    ("taulepdiscrtight",      po::value<bool>(&taulepdiscrtight)->default_value(false))
    ("dojerdebug",            po::value<bool>(&dojerdebug)->default_value(false))
    ("dotrgeff",              po::value<bool>(&dotrgeff)->default_value(false))
    ("do3dtrgeff",            po::value<bool>(&do3dtrgeff)->default_value(false))
    ("do1dparkedtrgeff",      po::value<bool>(&do1dparkedtrgeff)->default_value(false))
    ("dofitted1dparkedtrgeff",po::value<bool>(&dofitted1dparkedtrgeff)->default_value(false))
    ("dobinnedin2d1dtrgeff",  po::value<bool>(&dobinnedin2d1dtrgeff)->default_value(false))
    ("doidisoeff",            po::value<bool>(&doidisoeff)->default_value(false))
    ("dotopreweighting",      po::value<bool>(&dotopreweighting)->default_value(false))
    ("dopromptskim",          po::value<bool>(&dopromptskim)->default_value(false))
    ("donoskim",              po::value<bool>(&donoskim)->default_value(false))
    ("doidisoerr",            po::value<bool>(&doidisoerr)->default_value(false))
    ("doidisoerrupordown",    po::value<bool>(&doidisoerrupordown)->default_value(true))
    ("doidisoerrmuore",       po::value<bool>(&doidisoerrmuore)->default_value(true))
    ("dolumixsweight",        po::value<bool>(&dolumixsweight)->default_value(false))
    ("inputparams",           po::value<string>(&inputparams)->default_value("filelists/161031/Params161031.dat"))
    ("jettype",               po::value<string>(&jettype)->default_value("pfJetsPFlow"))
    ("trg_weight_file",       po::value<string>(&trg_weight_file)->default_value("input/scale_factors/TrigEff2016_MET1DFitHFBinned_errors_12d9fb.root"))
    ("trg_to_use",            po::value<string>(&trg_to_use)->default_value("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140"))
    ("printEventList",        po::value<bool>(&printEventList)->default_value(false))
    ("printEventContent",     po::value<bool>(&printEventContent)->default_value(false))
    ("dosmear",               po::value<bool>(&dosmear)->default_value(false))
    ("doType1Cor",            po::value<bool>(&doType1Cor)->default_value(false))
    ("doaltmatch",            po::value<bool>(&doaltmatch)->default_value(false))
    ("doetsmear",             po::value<bool>(&doetsmear)->default_value(false))
    ("dogaus",                po::value<bool>(&dogaus)->default_value(false))
    ("dospring10gaus",        po::value<bool>(&dospring10gaus)->default_value(false))
    ("jesuncfile",            po::value<string>(&jesuncfile)->default_value("input/jec/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt"))
    ("reapplyJEC",            po::value<bool>(&reapplyJEC)->default_value(false))
    ("jecdata",               po::value<string>(&jecdata)->default_value("input/jec/Spring16_25nsV6_DATA_L1FastJet_AK4PFchs.txt,input/jec/Spring16_25nsV6_DATA_L2Relative_AK4PFchs.txt,input/jec/Spring16_25nsV6_DATA_L3Absolute_AK4PFchs.txt,input/jec/Spring16_25nsV6_DATA_L2L3Residual_AK4PFchs.txt"))
    ("turnoffpuid",           po::value<bool>(&turnoffpuid)->default_value(false))
    ("useOldLT",              po::value<bool>(&useOldLT)->default_value(false))
    ("doTrigLT",              po::value<bool>(&doTrigLT)->default_value(false))
    ("doAllPairs",            po::value<bool>(&doAllPairs)->default_value(false))
    ("debug",                 po::value<unsigned>(&debug)->default_value(0))
    ("randomseed",            po::value<int>(&randomseed)->default_value(4357));

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);

  vector<string> filtersVec;
  boost::split(filtersVec, filters, boost::is_any_of(","));

  std::vector<string> jecdatafiles;
  boost::split(jecdatafiles, jecdata, boost::is_any_of(","));

  // Some options must now be re-configured based on other options
  ic::era era           = String2Era(era_str);
  ic::mc mc             = String2MC(mc_str);

  std::cout << "**** HiggsNuNu Analysis *****" << std::endl;
  string param_fmt = "%-25s %-40s\n";
  std::cout << boost::format(param_fmt) % "max_events" % max_events;
  std::cout << boost::format(param_fmt) % "output" % (output_folder+output_name);
  std::cout << boost::format(param_fmt) % "dojessyst" % dojessyst;
  if (dojessyst) std::cout << boost::format(param_fmt) % "jesupordown" % jesupordown;
  std::cout << boost::format(param_fmt) % "dojersyst" % dojersyst;
  if (dojessyst) std::cout << boost::format(param_fmt) % "jerbettorworse" % jerbetterorworse;
  std::cout << boost::format(param_fmt) % "era" % era_str;
  std::cout << boost::format(param_fmt) % "mc" % mc_str;
  std::cout << boost::format(param_fmt) % "prod" % prod;
  std::cout << boost::format(param_fmt) % "wstream" % wstream;
  std::cout << boost::format(param_fmt) % "doICHEP2016" % doICHEP2016;
  std::cout << boost::format(param_fmt) % "doLastJSON" % doLastJSON;
  std::cout << boost::format(param_fmt) % "dojetEtaCut" % dojetEtaCut;
  std::cout << boost::format(param_fmt) % "is_data" % is_data;
  std::cout << boost::format(param_fmt) % "doMetFilters" % doMetFilters;
  std::cout << boost::format(param_fmt) % "filters" % filters;
  std::cout << boost::format(param_fmt) % "dotrgeff" % dotrgeff;
  std::cout << boost::format(param_fmt) % "doidisoeff" % doidisoeff;
  std::cout << boost::format(param_fmt) % "jet1ptcut" % jet1ptcut;
  std::cout << boost::format(param_fmt) % "jet2ptcut" % jet2ptcut;
  std::cout << boost::format(param_fmt) % "mjjcut" % mjjcut;
  std::cout << boost::format(param_fmt) % "detajjcut" % detajjcut;


  // Load necessary libraries for ROOT I/O of custom classes
  // gSystem->Load("libFWCoreFWLite.dylib");
  // gSystem->Load("libUserCodeICHiggsTauTau.dylib");
  // AutoLibraryLoader::enable();
  bool is2012 = false;
  if (era_str.find("2012") != era_str.npos) is2012 = true;

  // Build a vector of input files
  vector<string> files = ParseFileLines(filelist);
  if (files.size()==0) {
    std::cout << " -- ERROR ! No input file provided, filelist " << filelist << " doesn't exist or is empty, please check." << std::endl;
    return 1;
  }
  for (unsigned i = 0; i < files.size(); ++i) files[i] = input_prefix + files[i];

  // Create ROOT output fileservice
  std::ofstream checkfile;
  checkfile.open((output_folder+"tmp.txt").c_str());
  if(!(checkfile.good())){
    std::cout<<"OUTPUT FOLDER DOESN'T EXIST: EXITING"<<std::endl;
    return 1;
  }
  checkfile.close();

  bool isW = output_name.find("JetsToLNu") != output_name.npos ||
    output_name.find("EWKW") != output_name.npos;


  //rename with _wstream.root properly
  std::string fullpath = output_folder+output_name;
  if (isW) {
    std::string suffix=".root";
    std::size_t found = fullpath.find(suffix);
    if(found!=std::string::npos){
      fullpath.erase(found,5);
    }
    fullpath = fullpath+"_"+wstream+".root";
  }
  fwlite::TFileService *fs = new fwlite::TFileService(fullpath.c_str());

  bool ignoreLeptons=false;
  if (output_name.find("iglep") != output_name.npos) {
    ignoreLeptons = true;
  }

  double elec_dz, elec_dxy;
  double muon_dz, muon_dxy;
  double veto_elec_dz, veto_elec_dxy;
  double veto_muon_dz, veto_muon_dxy;
  double elec_pt, elec_eta, muon_pt, muon_eta;
  double veto_elec_pt, veto_elec_eta, veto_muon_pt, veto_muon_eta;
  double loose_photon_pt, loose_photon_eta, medium_photon_pt, medium_photon_eta, tight_photon_pt, tight_photon_eta;

  double muon_iso = is2012 ? 0.12 : 0.15;
  double veto_muon_iso = is2012 ? 0.2 : 0.25;

  elec_dz = 0.1;
  elec_dxy = 0.02;
  veto_elec_dz = 0.2;
  veto_elec_dxy = 0.04;
  muon_dz = is2012 ? 0.2 : 0.5;
  muon_dxy = is2012 ? 0.045 : 0.2;
  veto_muon_dz = 0.5;
  veto_muon_dxy = 0.2;

  elec_pt = 20.0;
  elec_eta = 2.5;//2.4;
  muon_pt = 20.0;
  muon_eta = 2.4;//2.1;

  veto_elec_pt = 10.0;
  veto_elec_eta = 2.5;//2.4;
  veto_muon_pt = 10.0;
  veto_muon_eta = 2.4;//2.1;

  loose_photon_pt = 15;
  loose_photon_eta = 2.5;
  medium_photon_pt = 15;
  medium_photon_eta = 2.5;
  tight_photon_pt = 15;
  tight_photon_eta = 2.5;

  std::cout << "----------PARAMETERS----------" << std::endl;
  std::cout << boost::format("%-15s %-10s\n") % "elec_pt:" % elec_pt;
  std::cout << boost::format("%-15s %-10s\n") % "elec_eta:" % elec_eta;
  std::cout << boost::format("%-15s %-10s\n") % "veto_elec_pt:" % veto_elec_pt;
  std::cout << boost::format("%-15s %-10s\n") % "veto_elec_eta:" % veto_elec_eta;
  std::cout << boost::format("%-15s %-10s\n") % "elec_dxy:" % elec_dxy;
  std::cout << boost::format("%-15s %-10s\n") % "elec_dz:" % elec_dz;
  std::cout << boost::format("%-15s %-10s\n") % "veto_elec_dxy:" % veto_elec_dxy;
  std::cout << boost::format("%-15s %-10s\n") % "veto_elec_dz:" % veto_elec_dz;
  std::cout << boost::format("%-15s %-10s\n") % "muon_pt:" % muon_pt;
  std::cout << boost::format("%-15s %-10s\n") % "muon_eta:" % muon_eta;
  std::cout << boost::format("%-15s %-10s\n") % "veto_muon_pt:" % veto_muon_pt;
  std::cout << boost::format("%-15s %-10s\n") % "veto muon_eta:" % veto_muon_eta;
  std::cout << boost::format("%-15s %-10s\n") % "muon_dxy:" % muon_dxy;
  std::cout << boost::format("%-15s %-10s\n") % "muon_dz:" % muon_dz;
  std::cout << boost::format("%-15s %-10s\n") % "veto_muon_dxy:" % veto_muon_dxy;
  std::cout << boost::format("%-15s %-10s\n") % "veto_muon_dz:" % veto_muon_dz;
  std::cout << boost::format("%-15s %-10s\n") % "muon_iso:" % muon_iso;
  std::cout << boost::format("%-15s %-10s\n") % "veto_muon_iso:" % veto_muon_iso;



  // Create analysis object
  ic::AnalysisBase analysis(
    "HiggsNuNu",        // Analysis name
    files,                // Input files
    "icEventProducer/EventTree", // TTree name
    max_events);          // Max. events to process (-1 = all)

    analysis.StopOnFileFailure(false);
  // ------------------------------------------------------------------------------------
  // Misc Modules
  // ------------------------------------------------------------------------------------



  //print the event content
  HinvPrint hinvPrint("HinvPrint",is_data);

  //print run,lumi,evt of events selected
  HinvPrint hinvPrintList("HinvPrintList",is_data,false,true);

  //if (output_name.find("DYJJ01") != output_name.npos) {



  string data_json;
//   if (era == era::data_2011) data_json           =  "input/json/json_data_2011_et_mt.txt";
//   if (era == era::data_2012_ichep) data_json     =  "input/json/data_2012_ichep.txt";
//   if (era == era::data_2012_hcp) data_json       =  "input/json/data_2012_hcp.txt";
//   if (era == era::data_2012_moriond) data_json   =  "input/json/data_2012_moriond.txt";
//   if (era == era::data_2012_donly) data_json     =  "input/json/data_2012_donly.txt";

  std::string mydebugoutput("/home/hep/rd1715/CMSSW_8_0_20/src/UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/mydebugJSON");
  std::string suffix = output_name.substr( 0 , output_name.find(".root") );
  mydebugoutput.append(suffix);

  if (era == era::data_2016){
    if (doICHEP2016){ // 12d9 /fb
      data_json     =  "input/json/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt";
    }
    else if (doLastJSON){ // 36d8 /fb ReReco
      data_json     =  "input/json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
    }
  }

  LumiMask lumiMask = LumiMask("LumiMask")
    .set_produce_output_jsons(mydebugoutput.c_str())
    .set_input_file(data_json);

  MakeRunStats runStats = MakeRunStats("RunStats")
    .set_output_name(fullpath+".runstats");


  string mc_pu_file;
  if (mc == mc::fall11_42X) mc_pu_file   = "input/pileup/MC_Fall11_PU_S6-500bins.root";
  if (mc == mc::summer12_53X) mc_pu_file = "input/pileup/MC_Summer12_PU_S10-600bins.root";
  if (mc == mc::summer12_52X) mc_pu_file = "input/pileup/MC_Summer12_PU_S7-600bins.root";
  if (mc == mc::phys14_72X) mc_pu_file   = "input/pileup/MC_Summer12_PU_S10-600bins.root";//!!FIX WITH NEW PU
  if (mc == mc::spring15_74X) mc_pu_file = "input/pileup/MC_Spring15_PU25_Startup.root";
  if (mc == mc::fall15_76X) mc_pu_file   = "input/pileup/MC_Fall15_PU25_V1.root";
  if (mc == mc::spring16_80X) mc_pu_file = "input/pileup/MC_Spring16_PU25ns_V1.root";

  string data_pu_file;
  if (era == era::data_2012_rereco) data_pu_file  =  "input/pileup/Data_Pileup_2012_ReRecoPixel-600bins.root";
  if (era == era::data_2011) data_pu_file         =  "input/pileup/Data_Pileup_2011_HCP-500bins.root";
  if (era == era::data_2012_ichep) data_pu_file   =  "input/pileup/Data_Pileup_2012.root";
  if (era == era::data_2012_hcp) data_pu_file     =  "input/pileup/Data_Pileup_2012_HCP-600bins.root";
  if (era == era::data_2012_moriond) data_pu_file =  "input/pileup/Data_Pileup_2012_Moriond-600bins.root";
  if (era == era::data_2012_donly) data_pu_file   =  "input/pileup/Data_Pileup_2012_DOnly-600bins.root";
  if (era == era::data_2015_50ns) data_pu_file    =  "input/pileup/Data_Pileup_2012_ReRecoPixel-600bins.root";//!!FIX WITH NEW PU
  if (era == era::data_2015_25ns) data_pu_file    =  "input/pileup/Data_Pileup_mb69_2015D_246908-260627-600bins.root";
  if (era == era::data_2016){
    if (doICHEP2016){ // 12d9 /fb
      data_pu_file         =  "input/pileup/12d9/Data_Pileup_mb69d2_2016-600bins.root";
    }
    else if (doLastJSON){ // 36d8 /fb ReReco
      data_pu_file         =  "input/pileup/ReReco_Moriond17/Data_Pileup_mb69d2_2016-600bins.root";
    }
  }

  TH1D data_pu  = GetFromTFile<TH1D>(data_pu_file, "/", "pileup");
  TH1D mc_pu    = GetFromTFile<TH1D>(mc_pu_file, "/", "pileup");

  TH1D data_pu_up;
  TH1D data_pu_down;

  if(era==era::data_2012_moriond){
    data_pu_up  = GetFromTFile<TH1D>("input/pileup/Data_Pileup_2012_Moriond-600bins-Up.root", "/", "pileup");
    data_pu_down  = GetFromTFile<TH1D>("input/pileup/Data_Pileup_2012_Moriond-600bins-Down.root", "/", "pileup");
  }
  else if(era==era::data_2012_rereco){
    data_pu_up  = GetFromTFile<TH1D>("input/pileup/Data_Pileup_2012_ReRecoPixel-600bins-Up.root", "/", "pileup");
    data_pu_down  = GetFromTFile<TH1D>("input/pileup/Data_Pileup_2012_ReRecoPixel-600bins-Down.root", "/", "pileup");
  }
  else if(era == era::data_2015_50ns){
    data_pu_up  = GetFromTFile<TH1D>("input/pileup/Data_Pileup_2012_ReRecoPixel-600bins-Up.root", "/", "pileup");
    data_pu_down  = GetFromTFile<TH1D>("input/pileup/Data_Pileup_2012_ReRecoPixel-600bins-Down.root", "/", "pileup");
  }
  else if(era == era::data_2015_25ns){
    data_pu_up  = GetFromTFile<TH1D>("input/pileup/Data_Pileup_mb72d5_2015D_246908-260627-600bins.root", "/", "pileup");
    data_pu_down  = GetFromTFile<TH1D>("input/pileup/Data_Pileup_mb65d6_2015D_246908-260627-600bins.root", "/", "pileup");
  }
  else if(era == era::data_2016){
    if (doICHEP2016){ // 12d9 /fb
      data_pu_up    = GetFromTFile<TH1D>("input/pileup/12d9/Data_Pileup_mb72d4_2016-600bins.root", "/", "pileup");
      data_pu_down  = GetFromTFile<TH1D>("input/pileup/12d9/Data_Pileup_mb66_2016-600bins.root", "/", "pileup");
    }
    else if (doLastJSON){ // 36d77 /fb
      data_pu_up    = GetFromTFile<TH1D>("input/pileup/ReReco_Moriond17/Data_Pileup_mb72d4_2016-600bins.root", "/", "pileup");
      data_pu_down  = GetFromTFile<TH1D>("input/pileup/ReReco_Moriond17/Data_Pileup_mb66_2016-600bins.root", "/", "pileup");
    }
  }

  if (!is_data) {
    std::cout << "** Pileup Files **" << std::endl;
    std::cout << boost::format(param_fmt) % "mc_pu_file" % mc_pu_file;
    std::cout << boost::format(param_fmt) % "data_pu_file" % data_pu_file;
  }
  PileupWeight pileupWeight = PileupWeight("PileupWeight")
    .set_data(&data_pu)
    .set_mc(&mc_pu)
    .set_print_weights(false);
  PileupWeight pileupWeight_up = PileupWeight("PileupWeight_up","pileup_up")
    .set_weight_is_active(false)
    .set_data(&data_pu_up)
    .set_mc(&mc_pu)
    .set_print_weights(false);
  PileupWeight pileupWeight_down = PileupWeight("PileupWeight_down","pileup_down")
    .set_weight_is_active(false)
    .set_data(&data_pu_down)
    .set_mc(&mc_pu)
    .set_print_weights(false);

  //MAKE ele and mu eff weights like this
  // This is needed when you want to study the MET trigger efficiencies (Filter Flag_* as well) at the LT stage. Once the study is done, comment it again!
   // HinvDataTriggerFilter dataMCTriggerPathFilter("TriggerPathFilter");
   // dataMCTriggerPathFilter.set_is_data(is_data);
   // //dataMCTriggerPathFilter.set_trigger_path("HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v"); 
   // dataMCTriggerPathFilter.set_trigger_path("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140");
   // //dataMCTriggerPathFilter.set_trig_obj_label("triggerObjectsDiPFJet40PFMETnoMu65MJJ800VBFAllJets");
   // dataMCTriggerPathFilter.set_trig_obj_label("triggerObjectsDiPFJet40DEta3p5MJJ600PFMETNoMu140");
   
  // JetEnergyCorrections<PFJet> jetEnergyCorrections = JetEnergyCorrections<PFJet>
//   ("JetEnergyCorrections")
//   .set_input_label(jettype)
//   .set_is_data(is_data)
//   .set_l1_file("input/jec/START53_V10_L1FastJet_AK5PF.txt")
//   .set_l2_file("input/jec/START53_V10_L2Relative_AK5PF.txt")
//   .set_l3_file("input/jec/START53_V10_L3Absolute_AK5PF.txt");
  
  std::vector<string> inputVec;
  inputVec.push_back("input/halofilters/allevents.txt");
  MetEventFilters cscTightHaloFilter = MetEventFilters("CscTightHaloFilter",inputVec,doMetFilters);


  SimpleFilter<Vertex> goodVertexFilter = SimpleFilter<Vertex>("goodVertexFilter")
    .set_input_label("vertices")
    .set_predicate(bind(DummyFunction<Vertex>, _1))
    .set_min(1)
    .set_max(999);

  // ------------------------------------------------------------------------------------
  // Photon Modules
  // ------------------------------------------------------------------------------------
  //loose photons
  CopyCollection<Photon>  loosePhotonCopyCollection("CopyToLoosePhotons","photons","loosePhotons");
  
  ComplexFilter<Photon,EventInfo,double> loosePhotonFilter = ComplexFilter<Photon,EventInfo,double>
    ("LoosePhotonFilter")
    .set_primary_input_label("loosePhotons").set_predicate(bind(MinPtMaxEta, _1, loose_photon_pt, loose_photon_eta)&&
							   bind(LoosePhotonIDSpring15, _1,_2))
    .set_secondary_input_label("eventInfo").set_secondary_predicate(bind(&EventInfo::jet_rho,_1))
    .set_min(0)
    .set_max(999);
  
  //medium photons
  CopyCollection<Photon>  mediumPhotonCopyCollection("CopyToMediumPhotons","photons","mediumPhotons");
  
  ComplexFilter<Photon,EventInfo,double> mediumPhotonFilter = ComplexFilter<Photon,EventInfo,double>
    ("MediumPhotonPtEtaFilter")
    .set_primary_input_label("mediumPhotons").set_predicate(bind(MinPtMaxEta, _1, medium_photon_pt, medium_photon_eta)&&
							    bind(MediumPhotonIDSpring15, _1,_2))
    .set_secondary_input_label("eventInfo").set_secondary_predicate(bind(&EventInfo::jet_rho,_1))
    .set_min(0)
    .set_max(999);

  //tight photons
  CopyCollection<Photon>  tightPhotonCopyCollection("CopyToTightPhotons","photons","tightPhotons");
  
  ComplexFilter<Photon,EventInfo,double> tightPhotonFilter = ComplexFilter<Photon,EventInfo,double>
    ("TightPhotonPtEtaFilter")
    .set_primary_input_label("tightPhotons").set_predicate(bind(MinPtMaxEta, _1, tight_photon_pt, tight_photon_eta)&&
							   bind(TightPhotonIDSpring15, _1,_2))
    .set_secondary_input_label("eventInfo").set_secondary_predicate(bind(&EventInfo::jet_rho,_1))						   
    .set_min(0)
    .set_max(999);

  // ------------------------------------------------------------------------------------
  // Electron Modules
  // ------------------------------------------------------------------------------------

  // Electron Veto
  CopyCollection<Electron>  vetoElectronCopyCollection("CopyToVetoElectrons","electrons","vetoElectrons");

  ComplexFilter<Electron,EventInfo,double> vetoElectronFilter = ComplexFilter<Electron,EventInfo,double>
    ("VetoElectronPtEtaFilter")
    .set_primary_input_label("vetoElectrons").set_predicate(bind(MinPtMaxSCEta, _1, veto_elec_pt, veto_elec_eta) &&
							    //bind(VetoElectronFullID16, _1, _2)// && 
							    bind(VetoElectronFullIDSpring15, _1, _2)// && 
							    //bind(fabs, bind(&Electron::dxy_vertex, _1)) < veto_elec_dxy && 
							    //bind(fabs, bind(&Electron::dz_vertex, _1)) < veto_elec_dz
							    )
    .set_secondary_input_label("eventInfo").set_secondary_predicate(bind(&EventInfo::jet_rho,_1))						   
    .set_min(0)
    .set_max(999);

  EffectiveAreaIsolationFilter vetoElectronIso = EffectiveAreaIsolationFilter("VetoElectronIso","vetoElectrons",0.15);
 
  //electron selection 
  CopyCollection<Electron>  selElectronCopyCollection("CopyToSelElectrons","electrons","selElectrons");
  ComplexFilter<Electron,EventInfo,double> selElectronFilter = ComplexFilter<Electron,EventInfo,double>
    ("SelElectronPtEtaFilter")
    .set_primary_input_label("selElectrons").set_predicate(bind(MinPtMaxSCEta, _1, elec_pt, elec_eta) &&
							   //bind(TightElectronFullID16, _1, _2) //&&
							   bind(TightElectronFullIDSpring15, _1, _2) //&&
							   //bind(fabs, bind(&Electron::dxy_vertex, _1)) < elec_dxy && 
							   //bind(fabs, bind(&Electron::dz_vertex, _1)) < elec_dz
							   )
    .set_secondary_input_label("eventInfo").set_secondary_predicate(bind(&EventInfo::jet_rho,_1))						   
    .set_min(0)
    .set_max(999);

  EffectiveAreaIsolationFilter selElectronIso = EffectiveAreaIsolationFilter("SelElectronIso","selElectrons",0.10);

  OverlapFilter<Electron, Muon> vetoelecMuonOverlapFilter = OverlapFilter<Electron, Muon>("VetoElecMuonOverlapFilter")
    .set_input_label("vetoElectrons")
    .set_reference_label("vetoMuons")
    .set_min_dr(0.3);
  OverlapFilter<Electron, Muon> selelecMuonOverlapFilter = OverlapFilter<Electron, Muon>("SelElecMuonOverlapFilter")
    .set_input_label("selElectrons")
    .set_reference_label("vetoMuons")
    .set_min_dr(0.3);

  // ------------------------------------------------------------------------------------
  // Muon Modules
  // ------------------------------------------------------------------------------------

  // Muon Veto
  CopyCollection<Muon> vetoMuonCopyCollection("CopyToVetoMuons","muons","vetoMuons");

  SimpleFilter<Muon> vetoMuonFilter = SimpleFilter<Muon>
    ("VetoMuonPtEtaFilter")
    .set_input_label("vetoMuons")
    .set_predicate(bind(MinPtMaxEta, _1, veto_muon_pt, veto_muon_eta) &&
		   bind(MuonLoose, _1) &&
		   bind(PF04IsolationVal<Muon>, _1, 0.5, false) < veto_muon_iso
		   && bind(&Muon::HasIdIso,_1,"muon_idiso_badGlobalMuonTagger")==false
		   //&& bind(fabs, bind(&Muon::dxy_vertex, _1)) < veto_muon_dxy 
		   //&& bind(fabs, bind(&Muon::dz_vertex, _1)) < veto_muon_dz
		   )
    .set_min(0)
    .set_max(999);

  CopyCollection<Muon> vetoMuonNoIsoCopyCollection("CopyToVetoMuonsNoIso","muons","vetoMuonsNoIso");
  SimpleFilter<Muon> vetoMuonNoIsoFilter = SimpleFilter<Muon>
    ("VetoMuonNoIsoPtEtaFilter")
    .set_input_label("vetoMuonsNoIso").set_predicate(bind(MinPtMaxEta, _1, veto_muon_pt, veto_muon_eta) &&
						     bind(MuonLoose, _1)
						//&& bind(PF04Isolation<Muon>, _1, 0.5, 0.2)
						//&& bind(fabs, bind(&Muon::dxy_vertex, _1)) < veto_muon_dxy 
						//&& bind(fabs, bind(&Muon::dz_vertex, _1)) < veto_muon_dz
						)
    .set_min(0)
    .set_max(999);

  //sel muons
  CopyCollection<Muon> selMuonCopyCollection("CopyToSelMuons","muons","selMuons");
  SimpleFilter<Muon> selMuonFilter = SimpleFilter<Muon>
    ("SelMuonPtEtaFilter")
    .set_input_label("selMuons").set_predicate(bind(MinPtMaxEta, _1, muon_pt, muon_eta) &&
					       bind(MuonTight, _1) && 
					       bind(PF04IsolationVal<Muon>, _1, 0.5, false) < muon_iso
					       && bind(fabs, bind(&Muon::dxy_vertex, _1)) < muon_dxy && 
					       bind(fabs, bind(&Muon::dz_vertex, _1)) < muon_dz
					       && bind(&Muon::HasIdIso,_1,"muon_idiso_badGlobalMuonTagger")==false
					       )
    .set_min(0)
    .set_max(999);


  OneCollCompositeProducer<Muon> mumuLeadingPairProducer = OneCollCompositeProducer<Muon>
    ("MuMuLeadingPairProducer")
    .set_input_label("selMuons")
    .set_candidate_name_first("muon1")
    .set_candidate_name_second("muon2")
    .set_select_leading_pair(true)
    .set_output_label("mumuLeadingCandidates");                                                        

  SimpleFilter<CompositeCandidate> MassMuonPairFilter = SimpleFilter<CompositeCandidate>("MassMuonPairFilter")
    .set_input_label("mumuLeadingCandidates")
    .set_predicate( bind(PairMassInRange, _1,60,120) )
    .set_min(1)
    .set_max(999);    

  // ------------------------------------------------------------------------------------
  // Tau modules
  // ------------------------------------------------------------------------------------

  SimpleFilter<Tau> tauPtEtaFilter = SimpleFilter<Tau>("TauPtEtaFilter")
    .set_input_label("taus")
    .set_predicate(bind(MinPtMaxEta, _1, 18, 2.3))
    .set_min(0);

  SimpleFilter<Tau> tauDzFilter = SimpleFilter<Tau>("TauDzFilter")
    .set_input_label("taus")
    .set_predicate(bind(fabs, bind(&Tau::lead_dz_vertex, _1)) < 0.2)
    .set_min(0);

  CopyCollection<Tau>  vetoTauCopyCollection("CopyToVetoTaus","taus","vetoTaus");

  std::string tau_id_discr, tau_iso_discr, vetotau_iso_discr, tau_anti_elec_discr_1, tau_anti_elec_discr_2, tau_anti_muon_discr;
  if(do2016tauid){
    tau_id_discr          = "decayModeFinding";
    tau_iso_discr         = "byMediumIsolationMVArun2v1DBoldDMwLT";
    //tau_iso_discr         = "byLooseCombinedIsolationDeltaBetaCorr3Hits";
    vetotau_iso_discr         = "byCombinedIsolationDeltaBetaCorrRaw3Hits";
    tau_anti_muon_discr   = "againstMuonLoose3";
    tau_anti_elec_discr_1 = "againstElectronVLooseMVA6";
    tau_anti_elec_discr_2 = "againstElectronVLooseMVA6";
  }
  else if(do2015tauid){
    tau_id_discr          = "decayModeFindingNewDMs";
    tau_iso_discr         = "byTightCombinedIsolationDeltaBetaCorr3Hits";
    tau_anti_muon_discr   = "againstMuonTight3";
    tau_anti_elec_discr_1 = "againstElectronTightMVA6";
    tau_anti_elec_discr_2 = "againstElectronTightMVA6";
  }
  else{
    tau_id_discr          = "decayModeFinding";
    if(!docrosschecktau){
      tau_iso_discr       = "byTightCombinedIsolationDeltaBetaCorr3Hits";
    }
    else{
      tau_iso_discr       = "byMediumIsolationMVA2";
    }
    if(taulepdiscrtight){
      tau_anti_muon_discr   = "againstMuonTight";
      tau_anti_elec_discr_1 = "againstElectronTight";
      tau_anti_elec_discr_2 = "againstElectronTight";
    }
    else{
      tau_anti_muon_discr   = "againstMuonLoose";
      tau_anti_elec_discr_1 = "againstElectronLoose";
      tau_anti_elec_discr_2 = "againstElectronLoose";
    }
  }
  
  //tau_iso_discr         = "byLooseCombinedIsolationDeltaBetaCorr3Hits";
  
  std::cout << "** Tau Discriminators **" << std::endl;
  std::cout << boost::format(param_fmt) % "isolation" %  tau_iso_discr;
  std::cout << boost::format(param_fmt) % "isolation for veto taus" %  tau_iso_discr;
  std::cout << boost::format(param_fmt) % "anti-electron1" % tau_anti_elec_discr_1;
  std::cout << boost::format(param_fmt) % "anti-electron2" % tau_anti_elec_discr_2;
  std::cout << boost::format(param_fmt) % "anti-muon" % tau_anti_muon_discr;


  SimpleFilter<Tau> tauIsoFilter = SimpleFilter<Tau>("TauIsoFilter")
    .set_input_label("taus")
    .set_predicate((bind(&Tau::GetTauID, _1, tau_iso_discr) > 0.5) && (bind(&Tau::GetTauID, _1, tau_id_discr) > 0.5))
    .set_min(0);
  
  SimpleFilter<Tau> vetotauIsoFilter = SimpleFilter<Tau>("VetoTauIsoFilter")
    .set_input_label("vetoTaus")
    .set_predicate((bind(&Tau::GetTauID, _1, vetotau_iso_discr) < 5) && (bind(&Tau::GetTauID, _1, tau_id_discr) > 0.5))
    .set_min(0);

  SimpleFilter<Tau> tauElRejectFilter = SimpleFilter<Tau>("TauElRejectFilter")
    .set_predicate( (bind(&Tau::GetTauID, _1, tau_anti_elec_discr_1) > 0.5) && (bind(&Tau::GetTauID, _1, tau_anti_elec_discr_2) > 0.5) )        
    .set_input_label("taus")
    .set_min(0); 

  SimpleFilter<Tau> tauMuRejectFilter = SimpleFilter<Tau>("TauMuRejectFilter")
    .set_predicate(bind(&Tau::GetTauID, _1, tau_anti_muon_discr) > 0.5)
    .set_input_label("taus")
    .set_min(0);

  OverlapFilter<Tau, Muon> tauMuonOverlapFilter = OverlapFilter<Tau, Muon>("TauMuonOverlapFilter")
    .set_input_label("taus")
    .set_reference_label("vetoMuons")
    .set_min_dr(0.4);

  OverlapFilter<Tau, Electron> tauElecOverlapFilter = OverlapFilter<Tau, Electron>("TauElecOverlapFilter")
    .set_input_label("taus")
    .set_reference_label("vetoElectrons")
    .set_min_dr(0.4);

  OverlapFilter<Tau, Muon> vetotauMuonOverlapFilter = OverlapFilter<Tau, Muon>("VetoTauMuonOverlapFilter")
    .set_input_label("vetoTaus")
    .set_reference_label("vetoMuons")
    .set_min_dr(0.4);

  OverlapFilter<Tau, Electron> vetotauElecOverlapFilter = OverlapFilter<Tau, Electron>("VetoTauElecOverlapFilter")
    .set_input_label("vetoTaus")
    .set_reference_label("vetoElectrons")
    .set_min_dr(0.4);

  // ------------------------------------------------------------------------------------
  // Jet Modules
  // ------------------------------------------------------------------------------------  
  /*JECStudy jecStudy = JECStudy
    ("jecStudy")
    .set_input_label(jettype)
    .set_jec_old_files("input/jec/Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt,input/jec/Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt,input/jec/Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt,input/jec/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt")
    .set_jec_new_files("input/jec/Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt,input/jec/Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt,input/jec/Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt,input/jec/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt")
    .set_is_data(is_data)
    .set_reapplyJEC(true);
  */

  JetMETModifier ModifyJetMET = JetMETModifier
    ("ModifyJetMET")
    .set_input_label(jettype)
    .set_jec_data_files(jecdatafiles)
    .set_jec_mc_files(jecdatafiles)
    .set_met_label(mettype)
    .set_dogaus(dogaus)
    .set_dospring10gaus(dospring10gaus)
    .set_is_data(is_data)
    .set_jesuncfile(jesuncfile)
    .set_randomseed(randomseed)
    .set_run2(!is2012)
    .set_fs(fs);
  std::vector<JetMETModifier::jetmetCor> corVec;
  if (reapplyJEC && is_data) corVec.push_back(JetMETModifier::jetmetCor::jecData);
  if (reapplyJEC && !is_data) corVec.push_back(JetMETModifier::jetmetCor::jecMC);
  if (dosmear)    corVec.push_back(JetMETModifier::jetmetCor::smearMC);
  if (doType1Cor) corVec.push_back(JetMETModifier::jetmetCor::type1cor);
  ModifyJetMET.set_corVec(corVec);

  if (dojessyst) ModifyJetMET.set_syst(jesupordown?JetMETModifier::jetmetSyst::jesUp:JetMETModifier::jetmetSyst::jesDown);
  if (dojersyst) ModifyJetMET.set_syst(jerbetterorworse?JetMETModifier::jetmetSyst::jerBetter:JetMETModifier::jetmetSyst::jerWorse);
  if (douessyst) ModifyJetMET.set_syst(uesupordown?JetMETModifier::jetmetSyst::uesUp:JetMETModifier::jetmetSyst::uesDown);

  CopyCollection<PFJet> alljetsCopyCollection("copytoalljets",jettype,"AllpfJetsPFlow");

  SimpleFilter<PFJet> jetIDFilter = SimpleFilter<PFJet>
    ("JetIDFilter")
    .set_input_label(jettype);
    if(!turnoffpuid){
      jetIDFilter.set_predicate((bind(PFJetID2016, _1)) && bind(PileupJetID, _1,is2012?2:4,true));
      //jetIDFilter.set_predicate((bind(PFJetID2015, _1)) && bind(PileupJetID, _1,is2012?2:3,true));
    }
    else{
      jetIDFilter.set_predicate(bind(PFJetID2016, _1));
      //jetIDFilter.set_predicate(bind(PFJetID2015, _1));
    }

  // Jet pT eta filter
  SimpleFilter<PFJet> jetPtEtaFilter = SimpleFilter<PFJet>
    ("JetPtEtaFilter")
    .set_input_label(jettype).set_predicate(bind(MinPtMaxEta, _1, jetptprecut, 4.7));

  // Jet eta filter - cut on 3.0<eta<3.2
  SimpleFilter<PFJet> jetEtaFilter = SimpleFilter<PFJet>
    ("JetEtaFilter")
    .set_input_label(jettype).set_predicate(bind(EtaOutsideRange, _1, 3.0, 3.2));


  CJVFilter FilterCJV = CJVFilter("FilterCJV")
    .set_jetsinput_label(jettype)
    .set_pairinput_label("jjLeadingCandidates")
    .set_ptcut(30);
  
  //in principle no need to remove overlap because events with leptons are rejected...
  //except for specific e/mu selection for W background estimation.
  // and to calculate efficiencies of cuts correctly
  //we want the leading jets to really be jets even before vetoing leptons...
  //so: removal of jets matched with veto electrons and muons
  OverlapFilter<PFJet, Electron> jetElecOverlapFilter = OverlapFilter<PFJet, Electron>("jetElecOverlapFilter")
    .set_input_label(jettype)
    .set_reference_label("vetoElectrons")
    .set_min_dr(0.4);

  OverlapFilter<PFJet, Muon> jetMuonOverlapFilter = OverlapFilter<PFJet, Muon>("jetMuonOverlapFilter")
    .set_input_label(jettype)
    .set_reference_label("vetoMuons") //NoIso")
    .set_min_dr(0.4);

  OverlapFilter<PFJet, Tau> jetTauOverlapFilter = OverlapFilter<PFJet, Tau>("jetTauOverlapFilter")
    .set_input_label(jettype)
    .set_reference_label("taus") //NoIso")
    .set_min_dr(0.4);

  OneCollCompositeProducer<PFJet> jjLeadingPairProducer = OneCollCompositeProducer<PFJet>
    ("JetJetLeadingPairProducer")
    .set_input_label(jettype)
    .set_candidate_name_first("jet1")
    .set_candidate_name_second("jet2")
    .set_select_leading_pair(true)
    .set_output_label("jjLeadingCandidates");

  OneCollCompositeProducer<PFJet> jjAllPairProducer = OneCollCompositeProducer<PFJet>
    ("JetJetAllPairProducer")
    .set_input_label(jettype)
    .set_candidate_name_first("jet1")
    .set_candidate_name_second("jet2")
    .set_select_leading_pair(false)
    .set_output_label("jjAllCandidates");

  int npairs=1;
  if(donoskim)npairs=0;
  bool cutaboveorbelow=true;
  SimpleFilter<CompositeCandidate> jetPairFilter = SimpleFilter<CompositeCandidate>("JetPairFilter")
    .set_input_label("jjLeadingCandidates")
    .set_predicate( bind(OrderedPairPtSelection, _1,jet1ptcut, jet2ptcut, cutaboveorbelow) )
    .set_min(npairs)
    .set_max(999);

  if(!donoskim && !doAllPairs){
    jetPairFilter.set_predicate(bind(OrderedPairPtSelection, _1,jet1ptcut, jet2ptcut, cutaboveorbelow) && !bind(PairDEtaLessThan, _1, detajjcut) && bind(PairMassInRange, _1,mjjcut,100000) );
  }

  if (doAllPairs) {
    jetPairFilter.set_input_label("jjAllCandidates");
    jetPairFilter.set_predicate(!bind(PairDEtaLessThan, _1, 3.2) && bind(PairMassInRange, _1,600,100000) && bind(PairPtSelection, _1,30,30)  );
  }

  // ------------------------------------------------------------------------------------
  // Met Modules
  // ------------------------------------------------------------------------------------  

  SingleMetMaker singleMet = SingleMetMaker(mettype,mettype+"Collection");

  MetSelection metFilters = MetSelection("MetFilters",mettype,doMetFilters,filtersVec,0);
  

  unsigned nLepToAdd = 100;
  //no need to be explicit in the number of leptons: will take the number 
  // that is asked in the selection (i.e. exactly one for the W->e,mu selections...)
  
  ModifyMet metNoMuons     = ModifyMet("metNoMuons",mettype,"vetoMuons",2,nLepToAdd);
  ModifyMet metNoElectrons = ModifyMet("metNoElectrons",mettype,"vetoElectrons",1,nLepToAdd);
  ModifyMet metNoENoMu     = ModifyMet("metNoENoMu","metNoMuons","vetoElectrons",1,nLepToAdd);


  // ------------------------------------------------------------------------------------
  // Weight Modules
  // ------------------------------------------------------------------------------------  



  
  HinvWeights hinvWeights = HinvWeights("HinvWeights")
    .set_era(era)
    .set_mc(mc)
    .set_save_weights(true)
    .set_do_top_reweighting(dotopreweighting)
    .set_do_trg_weights(false)
    .set_trg_applied_in_mc(false)
    .set_do_idiso_tight_weights(false)
    .set_do_idiso_veto_weights(false)
    .set_do_idiso_err(doidisoerr)
    .set_do_idiso_errupordown(doidisoerrupordown)
    .set_do_idiso_errmuore(doidisoerrmuore)
    .set_fs(fs)
    .set_input_met("metNoMuons");
  if (!is_data) {
    std::vector<double> jptbinning;
    //jptbinning.push_back(70);
    //jptbinning.push_back(80);
    //jptbinning.push_back(1000);
    //for metmht
    jptbinning.push_back(0);
    jptbinning.push_back(7000);
    std::vector<double> mjjbinning;
    //mjjbinning.push_back(800);
    //mjjbinning.push_back(1000);
    //mjjbinning.push_back(10000);
    mjjbinning.push_back(0);
    mjjbinning.push_back(14000);
    hinvWeights.set_do_trg_weights(dotrgeff)
      .set_do_3dtrg_weights(do3dtrgeff)
      .set_do_1dparkedtrg_weights(do1dparkedtrgeff)
      .set_do_fitted1dparkedtrg_weights(dofitted1dparkedtrgeff)
      .set_do_binnedin2d1dfittedtrg_weights(dobinnedin2d1dtrgeff)
      .set_binnedin2d1dfitweightvar1binning(jptbinning)
      .set_binnedin2d1dfitweightvar2binning(mjjbinning)
      .set_do_run2(true)
      .set_do_metmht(true)
      .set_trg_weight_file(trg_weight_file)
      .set_trg_applied_in_mc(false);
    if(do3dtrgeff){
      hinvWeights.set_Alumi(0.889)
	.set_BClumi(11.023)
	.set_Dlumi(7.315);
    }
    hinvWeights.set_do_idiso_veto_weights(false);
    hinvWeights.set_do_idiso_tight_weights(false);
  }

  HinvWeights xsWeights = HinvWeights("XSWeights")
    .set_era(era)
    .set_mc(mc)
    .set_save_weights(false)
    .set_do_trg_weights(false)
    .set_trg_applied_in_mc(false)
    .set_do_idiso_tight_weights(false)
    .set_do_idiso_veto_weights(false)
    .set_do_lumixs_weights(dolumixsweight)
    .set_input_params(inputparams)
    .set_sample_name(output_name)
    .set_fs(fs);

//   Do not weight anymore, since we are dropping the inclusive sample due to MIT sync
//   If you switch it back, remember to change filelists/*/Params*.dat 
//   if (output_name.find("JetsToLNu") != output_name.npos &&
//       output_name.find("EWKW") == output_name.npos &&
//       output_name.find("nlo") == output_name.npos) {
// 
//     SetDoW( mc, &xsWeights );
//   }

  if (output_name.find("JetsToLL-mg-m50") != output_name.npos &&
      output_name.find("PtZ-100-madgraph") == output_name.npos &&
      output_name.find("Zpt150") == output_name.npos &&
      output_name.find("DYJJ01") == output_name.npos &&
      output_name.find("m50-ht") == output_name.npos) {

    SetDoDY( mc, &xsWeights );
  }

  if (output_name.find("JetsToLL-mg-m50-ht") != output_name.npos ||
      output_name.find("JetsToNuNu") != output_name.npos) {
    xsWeights.set_do_dy_reweighting(true);
  }

  if (output_name.find("JetsToLNu-mg-ht") != output_name.npos) {
    xsWeights.set_do_w_reweighting(true);
  }

  /*if (output_name.find("JetsToLL") != output_name.npos &&
      output_name.find("m50-ht") != output_name.npos){
    xsWeights.set_do_dy_soup_htbinned(false);
    xsWeights.set_do_dy_reweighting(false);
    xsWeights.SetDYTargetFractions(1,1,1,1,1);
    xsWeights.SetDYInputYields(9004328,0,962195,1069003,1031103);
    }*/

  // ------------------------------------------------------------------------------------
  // Gen particle selection modules
  // ------------------------------------------------------------------------------------

  int lFlavour = 1;
  if (wstream == "enu") lFlavour = 11;
  else if (wstream == "munu") lFlavour = 13;
  else if (wstream == "taunu") lFlavour = 15;

  WDecaySeparator WtoLeptonFilter = WDecaySeparator("WtoLeptonSelector",lFlavour);

  HinvWDecay WtoLeptonFilter2012 = HinvWDecay("WtoLeptonSelector2012",lFlavour);
  WtoLeptonFilter2012.set_do_newstatuscodes(false);
  WtoLeptonFilter2012.set_do_statusflags(false);

  HinvZDecay ZhighptFilter = HinvZDecay("ZhighptFilter",is2012?13:0,0,150,HinvZDecay::cutVar::pt,is2012);
  HinvZDecay ZhighhtFilter = HinvZDecay("ZhighhtFilter",is2012?13:0,0,100,HinvZDecay::cutVar::ht,is2012);

  // ------------------------------------------------------------------------------------
  // Plot Modules
  // ------------------------------------------------------------------------------------  

  LightTree lightTree = LightTree("LightTree")
    .set_fs(fs)
    .set_met_label(mettype)
    .set_jet_label(jettype)
    .set_dijet_label("jjLeadingCandidates")
    .set_sel_label("JetPair")
    .set_is_data(is_data)
    .set_dotrigskim(false)
    .set_do_promptskim(dopromptskim)
    .set_do_noskim(donoskim)
    .set_ignoreLeptons(ignoreLeptons)
    .set_trigger_path("HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v");
  //.set_trig_obj_label("triggerObjectsPFMET170NoiseCleaned");

  LightTreeTrig lightTreeTrig = LightTreeTrig("LightTreeTrig")
    .set_fs(fs)
    .set_met_label(mettype)
    .set_jet_label(jettype)
    .set_dijet_label("jjLeadingCandidates")
    .set_is_data(is_data)
    .set_dotrigskim(false)
    .set_do_noskim(donoskim)
    .set_trigger_path("HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v");
  //.set_trig_obj_label("triggerObjectsPFMET170NoiseCleaned");

//  LightTreeAM lightTreeNew = LightTreeAM("LightTreeNew")
  LightTreeRDM lightTreeNew = LightTreeRDM("LightTreeNew")
    .set_fs(fs)
    .set_debug(debug)
    .set_met_label(mettype)
    .set_jet_label(jettype)
    .set_dijet_label("jjLeadingCandidates")
    .set_is_data(is_data)
    .set_do_trigskim(false)
    .set_do_promptskim(dopromptskim)
    .set_do_noskim(donoskim)
    .set_ignoreLeptons(ignoreLeptons);

  if (doAllPairs) lightTreeNew.set_dijet_label("jjAllCandidates");

  // ------------------------------------------------------------------------------------
  // Build Analysis Sequence
  // ------------------------------------------------------------------------------------  

  if (is_data) analysis.AddModule(&lumiMask);

  //fix for single elec files
  if (!doType1Cor) analysis.AddModule(&singleMet);
  if (!is_data) {
    //do W streaming to e,mu,tau
    if (isW) {
      if (wstream != "nunu") {
	if (is2012) analysis.AddModule(&WtoLeptonFilter2012);
	else analysis.AddModule(&WtoLeptonFilter);
      }
    }
    //Do PU Weight
    analysis.AddModule(&pileupWeight);
    analysis.AddModule(&pileupWeight_up);
    analysis.AddModule(&pileupWeight_down);
    //just apply W and Z weights
    analysis.AddModule(&xsWeights);
    //Z pt <100 GeV cut for inclusive DY samples
    if(doincludehighptz && output_name.find("JetsToLL-mg-m50") != output_name.npos && output_name.find("Zpt150") == output_name.npos){
      analysis.AddModule(&ZhighptFilter);
    }
    if(output_name.find("JetsToLL") != output_name.npos && output_name.find("nlo") == output_name.npos && output_name.find("-ht") == output_name.npos && output_name.find("Zpt150") == output_name.npos){
      analysis.AddModule(&ZhighhtFilter);
    }
  }

  //if (printEventList) analysis.AddModule(&hinvPrintList);
  if (true) { //it was "is_data"; now changed to apply on MC badChargedHadronFilter & badMuonFilter; but be careful to change correctly the cfg file
    // see line 459
    //analysis.AddModule(&dataMCTriggerPathFilter);
    analysis.AddModule(&metFilters);
    //74X only
    //analysis.AddModule(&cscTightHaloFilter);
  }

  if(!donoskim) analysis.AddModule(&goodVertexFilter);
  //jet modules
  //Module to do jet smearing and systematics
  //analysis.AddModule(&jecStudy);
  analysis.AddModule(&ModifyJetMET);

  analysis.AddModule(&jetIDFilter);
  //don't want pile-up jets to calculate HT,MHT...
  analysis.AddModule(&alljetsCopyCollection);

  //prepare collections of photons
  //no photons in tautau single ele data files...
  if (!doType1Cor){
    analysis.AddModule(&loosePhotonCopyCollection);
    analysis.AddModule(&loosePhotonFilter);
    analysis.AddModule(&mediumPhotonCopyCollection);
    analysis.AddModule(&mediumPhotonFilter);
    analysis.AddModule(&tightPhotonCopyCollection);
    analysis.AddModule(&tightPhotonFilter);
  }

  //prepare collections of veto leptons
  analysis.AddModule(&vetoElectronCopyCollection);
  analysis.AddModule(&vetoElectronFilter);
  //analysis.AddModule(&vetoElectronIso);
  analysis.AddModule(&vetoMuonCopyCollection);
  analysis.AddModule(&vetoMuonFilter);

  //filter leptons before making jet pairs and changing MET...
  analysis.AddModule(&selElectronCopyCollection);
  analysis.AddModule(&selElectronFilter);
  //analysis.AddModule(&selElectronIso);
  analysis.AddModule(&selMuonCopyCollection);
  analysis.AddModule(&selMuonFilter);
  analysis.AddModule(&vetoelecMuonOverlapFilter);
  analysis.AddModule(&selelecMuonOverlapFilter);

  //filter taus for plots
  analysis.AddModule(&tauPtEtaFilter);

  //if (printEventList) analysis.AddModule(&hinvPrintList);


  //deal with removing overlap with selected leptons
  analysis.AddModule(&jetMuonOverlapFilter);
  analysis.AddModule(&jetElecOverlapFilter);
  //no need to clean taus, we don't do it in the signal selection.
  //analysis.AddModule(&jetTauOverlapFilter);


  //add met without leptons for plots
  analysis.AddModule(&metNoMuons);
  analysis.AddModule(&metNoElectrons);
  analysis.AddModule(&metNoENoMu);

  //filter taus
  analysis.AddModule(&vetoTauCopyCollection);
  analysis.AddModule(&vetotauIsoFilter);
  analysis.AddModule(&vetotauMuonOverlapFilter);
  analysis.AddModule(&vetotauElecOverlapFilter);

  analysis.AddModule(&tauDzFilter);
  analysis.AddModule(&tauIsoFilter);
  analysis.AddModule(&tauElRejectFilter);
  analysis.AddModule(&tauMuRejectFilter);
  analysis.AddModule(&tauMuonOverlapFilter);
  analysis.AddModule(&tauElecOverlapFilter);

  //filter jets
  analysis.AddModule(&jetPtEtaFilter);
  if(dojetEtaCut){
    analysis.AddModule(&jetEtaFilter);
  }

  //if (printEventContent) analysis.AddModule(&hinvPrint);
  //two-leading jet pair production before plotting
  analysis.AddModule(&jjLeadingPairProducer);
  if (doAllPairs) analysis.AddModule(&jjAllPairProducer);

  //if (printEventContent) analysis.AddModule(&hinvPrint);

  //analysis.AddModule(&hinvPrint);

  //Need two jets and metnomuons to apply trigger weights.
  //need sel leptons to apply idiso weights
  if (!is_data) analysis.AddModule(&hinvWeights);

    //if (printEventList) analysis.AddModule(&hinvPrintList);

  //record the number of jets in the gap
  analysis.AddModule(&jetPairFilter);
  analysis.AddModule(&FilterCJV);

  //jet pair selection
  //if (printEventList) analysis.AddModule(&hinvPrintList);

  //write tree
  if (useOldLT) analysis.AddModule(&lightTree);
  else analysis.AddModule(&lightTreeNew);

  if (!useOldLT && doTrigLT) analysis.AddModule(&lightTreeTrig);

  // Run analysis
  analysis.RetryFileAfterFailure(60,5);// int <pause between attempts in seconds>, int <number of retry attempts to make> );
  analysis.StopOnFileFailure(true);
  analysis.SetTTreeCaching(true); 

  analysis.RunAnalysis();
  delete fs;
  return 0;

}

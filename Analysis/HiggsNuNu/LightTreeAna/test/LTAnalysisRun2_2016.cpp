#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HiggsNuNuAnalysisTools.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeAnalyser.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeModule.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/LightTreeFiles.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataNormShape.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataZNormShape.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataShape.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/TrigEff.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataZEst.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/NormPlots.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/SimplePlots.h"
//#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/MVATrain.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/AddFriends.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/Plotter.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/HistPlotter.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/SummaryTable.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/EventList.h"
//#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/deltaphi.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/BkgSubDataNormShape.h"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "TColor.h"
#include "TMath.h"

namespace po=boost::program_options;
using namespace ic;

double deltaPhi(double phi1,double phi2)
{
  double result = phi1 - phi2;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi()) result += 2*TMath::Pi();
  return result;
};
double deltaR(double eta1,double phi1,double eta2,double phi2)
{
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
};
double deltaRmin(double eta1,double phi1,double eta2,double phi2,double eta3,double phi3)
{
  return std::min(deltaR(eta1,phi1,eta3,phi3),deltaR(eta2,phi2,eta3,phi3));
};

double getPostFitSF(std::string channel, std::string process){

  if (process=="error") return -1; // to be removed after test
  else return 1; // to be removed after test

  if (channel=="enu"){
    if (process=="wewk") return 7.824/7.31107;
    else if (process=="wqcd") return 16.32675/15.2317;
    else if (process=="top") return 1.42998/1.43447;
    else if (process=="VV") return 0.37785/0.375499;
    else if (process=="QCD") return 1.;
    else if (process=="error") return 2.82507/26.0124;
  } else if (channel=="munu"){
    if (process=="wewk") return 27.464/26.2531;
    else if (process=="wqcd") return 52.87799/50.5378;
    else if (process=="top") return 6.8237/6.82765;
    else if (process=="VV") return 0.79000/0.775874;
    else if (process=="QCD") return 4.5360/5.58472;
    else if (process=="error") return 8.6861979/92.4918;
  }else if (channel=="taunu"){
    if (process=="wewk") return 5.23699/5.09101;
    else if (process=="wqcd") return 11.4885/11.1967;
    else if (process=="top") return 7.13003/7.10643;
    else if (process=="VV") return 1.;
    else if (process=="QCD") return 0.4229/0.519936;
    else if (process=="error") return 3.210465/25.063774;
  } else if (channel=="mumu"){
    if (process=="zewk") return 2.10208/1.94073;
    else if (process=="zqcd") return 3.7217/3.37843;
    else if (process=="top") return 0.218995/0.219786;
    else if (process=="VV") return 0.018013/0.0173571;
    else if (process=="QCD") return 1.;
    else if (process=="error") return 1.3153597/6.060841;
  } else if (channel=="nunu"){
    if (process=="wewk") return (5.53635+4.37434+5.617556)/(5.20467+4.17524+5.41536);
    else if (process=="wqcd") return (9.504317+13.400+12.8982)/(8.81025+12.836+12.3952);
    else if (process=="zewk") return 21.9927864/20.3972;
    else if (process=="zqcd") return 40.42113/36.7616;
    else if (process=="top") return 2.2525/2.25912;
    else if (process=="VV") return 0.6799356/0.680839;
    else if (process=="QCD") return 3.08484/3.80019;
    else if (process=="error") return 27.1030/119.76226;
  }

  return 1;
};


int main(int argc, char* argv[]){
  /*##########################################
  #                                          #
  #            SET UP OPTIONS                #
  #                                          #
  ##########################################*/

  //Input output and config options
  std::string cfg;
  std::string outputname;
  std::string inputfolder;
  std::string eos_path_mc;
  std::string eos_path_data;
  std::string inputparams;
  std::string filelist;
  std::string basesel;
  std::string baseselele;

  std::string channel;
  std::string syst;
  std::string dataset;

  bool do_mettrig;
  bool apply_trig_in_mc;
  bool runblind;
  bool runblindreg;
  bool blindcutreg;

  bool dataonly;
  bool datalist;
  bool do_list;

  bool do_latex;
  bool do_logy;
  bool do_mcbkg;
  bool use_nlo;
  bool do_tauveto;
  bool do_bveto;
  double met_cutval;
  bool do_lep_mt_cut;
  bool do_MIT_UCSD_sync;
  bool do_new_muVeto;

  std::string jetmetdphicut;
  std::string metsigcut;


  std::string histTitlePar;
  std::string shapePar;

  unsigned debug;

  double lumiSF;

  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options
    ("output_name,o",            po::value<std::string>(&outputname)->default_value("tmp.root"))
    ("input_folder,i",           po::value<std::string>(&inputfolder)->default_value(""))
    ("eos_path_mc",              po::value<std::string>(&eos_path_mc)->default_value(""))
    ("eos_path_data",            po::value<std::string>(&eos_path_data)->default_value(""))
    ("syst,s",                   po::value<std::string>(&syst)->default_value(""))
    ("input_params,p",           po::value<std::string>(&inputparams)->default_value(""))
    ("filelist,f",               po::value<std::string>(&filelist)->default_value("filelists/run2filelist80X.dat"))
    ("dataset,d",                po::value<std::string>(&dataset)->default_value("MET"))
    ("dataonly",                 po::value<bool>(&dataonly)->default_value(false))
    ("datalist",                 po::value<bool>(&datalist)->default_value(false))
    ("do_list",                  po::value<bool>(&do_list)->default_value(false))
    ("basesel",                  po::value<std::string>(&basesel)->default_value(""))
    ("baseselele",               po::value<std::string>(&baseselele)->default_value(""))
    ("channel",                  po::value<std::string>(&channel)->default_value("nunu"))

    ("runblind",                 po::value<bool>(&runblind)->default_value(true))

    ("do_latex",                 po::value<bool>(&do_latex)->default_value(false))

    ("do_mettrig",               po::value<bool>(&do_mettrig)->default_value(false))
    ("apply_trig_in_mc",         po::value<bool>(&apply_trig_in_mc)->default_value(false))

    ("do_logy",                  po::value<bool>(&do_logy)->default_value(false))
    ("blindcutreg",              po::value<bool>(&blindcutreg)->default_value(true))
    ("runblindreg",              po::value<bool>(&runblindreg)->default_value(true))
    ("debug",                    po::value<unsigned>(&debug)->default_value(0))
    ("do_mcbkg",                 po::value<bool>(&do_mcbkg)->default_value(true))
    ("use_nlo",                  po::value<bool>(&use_nlo)->default_value(false))
    ("jetmetdphicut",            po::value<std::string>(&jetmetdphicut)->default_value(""))
    ("metsigcut",                po::value<std::string>(&metsigcut)->default_value(""))
    ("histTitlePar",             po::value<std::string>(&histTitlePar)->default_value(";#Delta#phi(E_{T}^{miss},j);Events"))
    ("shapePar",                 po::value<std::string>(&shapePar)->default_value("fourjetsmetnomu_mindphi(32,0.,3.1416)"))
    ("lumiSF",                   po::value<double>(&lumiSF)->default_value(1.0))
    ("do_tauveto",               po::value<bool>(&do_tauveto)->default_value(false))
    ("do_bveto",                 po::value<bool>(&do_bveto)->default_value(false))
    ("met_cutval",               po::value<double>(&met_cutval)->default_value(0))
    ("do_lep_mt_cut",            po::value<bool>(&do_lep_mt_cut)->default_value(false))
    ("do_MIT_UCSD_sync",         po::value<bool>(&do_MIT_UCSD_sync)->default_value(false))
    ("do_new_muVeto",            po::value<bool>(&do_new_muVeto)->default_value(false))
    ;

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);

  //this is making plots in scripts
  std::vector<std::string> histTitle;
  boost::split(histTitle, histTitlePar, boost::is_any_of("!"));
  std::vector<std::string> shape;
  boost::split(shape, shapePar, boost::is_any_of("!"));

  std::cout << " -- input shape string: " << shapePar << std::endl;
  std::cout << " -- Processing " << shape.size() << " shapes." << std::endl;
  for (unsigned iS(0);iS<shape.size();++iS){
    //shape[iS] = "\""+shape[iS]+"\"";
    std::cout << shape[iS] << std::endl;
    //histTitle[iS] = "\""+histTitle[iS]+"\"";
  }


  /*  //trigger systematics
  TFile *mettrigsyst = TFile::Open("../input/scale_factors/mettrigger.root");
  if (!mettrigsyst) {
    std::cout << " -- Error, file for met trigger systematics not found." << std::endl;
    return 1;
  }
  mettrigsyst->cd();
  TH1D *hmm = (TH1D*)gDirectory->Get("zmm_sys");
  TH1D *hvv = (TH1D*)gDirectory->Get("zvv_sys");
  TH1D *htrigsyst = (TH1D*)hvv->Clone("htrigsyst");
  for(int iBin = 0; iBin < htrigsyst->GetNbinsX()+1; iBin++){
    if(htrigsyst->GetBinCenter(iBin+1) < 500)
      htrigsyst->SetBinContent(iBin+1,1-hmm->GetBinContent(hmm->FindBin(htrigsyst->GetBinCenter(iBin+1)))/hvv->GetBinContent(iBin+1));
    else
      htrigsyst->SetBinContent(iBin+1,0.);
  }
  */

  /*##########################################
  #                                          #
  #          INSTANTIATE ANALYSER            #
  #                                          #
  ##########################################*/

  LTAnalyser* analysis = new LTAnalyser(outputname);


  if (debug) std::cout << "syst=" << syst << " size " << syst.size() << std::endl;

  analysis->SetEosFolders(eos_path_data,eos_path_mc);

  analysis->AddFiles(filelist);
  if(syst.find("TAU")==syst.npos&&syst.find("BTAG")==syst.npos&&syst.find("LEPEFF")==syst.npos&&syst!="PUUP"&&syst!="PUDOWN"&&syst.find("TRIG")==syst.npos&&syst.size()!=0){
    std::cout<<"Syst, taking input from: "<<inputfolder<<"/"<<syst<<std::endl;
    analysis->SetInFolder(inputfolder+"/"+syst);
  }
  else{
    std::cout<<"Taking input from: "<<inputfolder<<std::endl;
    analysis->SetInFolder(inputfolder);
  }
  analysis->SetInputParams(inputparams);

  std::cout << " -- Base selection:               " << basesel    << std::endl;
  std::cout << " -- Base selection for electrons: " << baseselele << std::endl;

  //add metsig cut
  std::string metsigsel;
  if(channel=="mumu") metsigsel="&& metnomuons/sqrt(sumet-mu1_pt-mu2_pt)>"+metsigcut;
  else if (channel=="munu") metsigsel="&& metnomuons/sqrt(sumet-mu1_pt)>"+metsigcut;
  else metsigsel="&& metnomuons/sqrt(sumet)>"+metsigcut;

  if (channel=="ee" || channel=="enu") basesel = baseselele;//+metsigsel;
  else basesel = basesel;//+metsigsel;

  analysis->set_baseselection(basesel);

  /*##########################################
  #                                          #
  #            DEFINE MODULES                #
  #                                          #
  ##########################################*/

  std::string dataextrasel;
  std::string mcextrasel;

  if (channel=="ee" || channel=="enu"){
    dataextrasel="&&(pass_singleEltrigger==1)";
  }
  else if (channel=="qcdA" || channel=="qcdC"){
    dataextrasel="&&(pass_alljethttrigger>=1)";
    dataset="JetHT";
  }
  else if (channel!="gamma"){
    //if(!do_mettrig) dataextrasel="&&(pass_sigtrigger==1)";
    //for metmht trigger
    //HACK
    dataextrasel="&&(pass_metmht90trigger==1 || pass_metmht100trigger==1 || pass_metmht110trigger==1 || pass_metmht120trigger==1 || pass_mettrigger==1)";
//     if(!do_mettrig){
//       dataextrasel="&&(pass_metmht90trigger==1 || pass_metmht100trigger==1 || pass_metmht110trigger==1 || pass_metmht120trigger==1)";
//     }
//     else {
//       dataextrasel="&&(pass_mettrigger==1)";
//     }
  }
  else {
    dataextrasel="&&(pass_photontrigger==1)";
  }

  if (apply_trig_in_mc) mcextrasel=dataextrasel;

  std::string sigcat;
  std::string sigcat_novetotaus;
  std::string zextrasigcat;

  std::string tauveto,tauvetoweight;
  std::string bveto,bvetoweight;
  std::string met_cut;
  std::string lep_mt_cut;
  std::string jet1_ID;


  if (do_tauveto){
    if (channel!="taunu") tauveto="&&nvetotaus==0";
    else tauveto="";
//     tauveto="";
    if (channel!="taunu") dataextrasel += "&&nvetotaus==0";
    if (syst=="TAUUP") tauvetoweight="*TMath::Power(1-0.99,nvetotaus)*(1+5*nvetotaus)";//"*1.";"*TMath::Power(1-0.97,nvetotaus)";
    else if (syst=="TAUDOWN") tauvetoweight="*TMath::Power(1-0.99,nvetotaus)*(1-5*nvetotaus)";//"*TMath::Power(1-0.9,nvetotaus)";
    else tauvetoweight="*TMath::Power(1-0.99,nvetotaus)";
  } else {
    tauveto="";
    tauvetoweight="";
  }
  if (do_bveto){
    //bveto="&&n_jets_csv2medium==0";
    //bveto="&&( n_jets_30 < 3 || ( n_jets_30>2 && jet3_csv<0.8484 && ( n_jets_30<4 || ( n_jets_30>3 && jet4_csv<0.8484 ) ) ) )";
    bveto="";
    dataextrasel+="&&n_jets_csv2medium==0";
    if (syst=="BTAGUP") bvetoweight="*weight_0bup_alljets";
    else if (syst=="BTAGDOWN") bvetoweight="*weight_0bdown_alljets";
    else bvetoweight="*weight_0b_alljets";
  } else {
    bveto="";
    bvetoweight="";
  }

  if (met_cutval>0){
    std::ostringstream ltmp;
    ltmp << "&&met>" << met_cutval;
    met_cut=ltmp.str();
  } else {
    met_cut="";
  }
  if (do_lep_mt_cut){
    lep_mt_cut="&&lep_mt<160";
  } else {
    lep_mt_cut="";
  }
  if (do_MIT_UCSD_sync){
    jet1_ID="&&( !( abs(jet1_eta) < 2.4 && !(jet1_chargedhadfrac>0.1 && jet1_neutralhadfrac<0.8) ) )";
  } else {
    jet1_ID="";
  }

  //HACK
  std::string nunucat;
  if (do_new_muVeto) nunucat  = "nvetoelectrons==0&&"+jetmetdphicut+bveto+jet1_ID;
  else               nunucat  = "nvetomuons==0&&nvetoelectrons==0&&"+jetmetdphicut+bveto+jet1_ID;

  std::string enucat   = "nselelectrons==1&&nvetomuons==0&&nvetoelectrons==1&&ele1_pt>40&&"+jetmetdphicut+bveto+lep_mt_cut+met_cut+jet1_ID;
  std::string munucat  = "nselmuons==1&&nvetomuons==1&&nvetoelectrons==0&&lep_mt>=0&&"+jetmetdphicut+bveto+lep_mt_cut+jet1_ID;

  std::string eecat    = "nselelectrons>=1&&nvetoelectrons==2&&nvetomuons==0&&m_ee>60&&m_ee<120&&oppsign_ee&&ele1_pt>40&&"+jetmetdphicut+bveto+jet1_ID;
  std::string mumucat  = "nselmuons>=1&&nvetomuons==2&&nvetoelectrons==0&&m_mumu>60&&m_mumu<120&&oppsign_mumu&&"+jetmetdphicut+bveto+jet1_ID;

  std::string taunucat = "ntaus==1&&nvetomuons==0&&nvetoelectrons==0&&"+jetmetdphicut+bveto+lep_mt_cut+jet1_ID;
  std::string gammacat = "ntightphotons==1&&nvetomuons==0&&nvetoelectrons==0&&"+jetmetdphicut+bveto+jet1_ID;
  std::string toplcat  = "nvetomuons==1&&nvetoelectrons==1&&nselmuons==1&&nselelectrons==1"+bveto+jet1_ID;
  std::string topbcat  = "(nselmuons>=1 || nselelectrons>=1)&&(jet1_csv>0.679||jet2_csv>0.679)&&(forward_tag_eta>2.8||forward_tag_eta<-2.8)"+bveto+jet1_ID;

  if(channel=="nunu" || channel.find("qcd")!=channel.npos){//nunu
    sigcat=nunucat+tauveto;
    sigcat_novetotaus=nunucat;
  }
  else {
    if(channel=="mumu"){//zmumu
      sigcat=mumucat+tauveto;
      sigcat_novetotaus=mumucat;
    }
    else if(channel=="ee"){//zee
      sigcat=eecat+tauveto;
      sigcat_novetotaus=eecat;
    }
    else if(channel=="munu"){//wmu
      sigcat=munucat+tauveto;
      sigcat_novetotaus=munucat;
    }
    else if(channel=="enu"){//wel
      sigcat=enucat+tauveto;
      sigcat_novetotaus=enucat;
    }
    else if(channel=="taunu"){//wtau
      sigcat=taunucat+tauveto;
      sigcat_novetotaus=taunucat;
    }
    else if(channel=="gamma"){
      sigcat=gammacat+tauveto;
      sigcat_novetotaus=gammacat;
    }
    else if(channel=="topl"){
      sigcat=toplcat+tauveto;
      sigcat_novetotaus=toplcat;
    }
    else if(channel=="topb"){
      sigcat=topbcat+tauveto;
      sigcat_novetotaus=topbcat;
    }
     else{
      std::cout<<"Error: Channel "<<channel<<" not recognised, exiting"<<std::endl;
      return 1;
    }
  }

  std::string sigmcweight;
  std::string dataweight;

  std::ostringstream mcweightsystfactor;
  mcweightsystfactor << "*" << lumiSF;
  mcweightsystfactor << bvetoweight;
//   if (channel!="taunu") mcweightsystfactor << tauvetoweight;
  if(syst=="PUUP") mcweightsystfactor << "*puweight_up_scale";
  if(syst=="PUDOWN") mcweightsystfactor << "*puweight_down_scale";

  //if (syst=="TRIGUP" && channel!="ee" && channel!="enu") mcweightsystfactor<<"*weight_trig_1/weight_trig_0";
  //if (syst=="TRIGDOWN" && channel!="ee" && channel!="enu") mcweightsystfactor<<"*weight_trig_2/weight_trig_0";

  if (syst=="TRIGUP" && (channel=="ee" || channel=="enu")) mcweightsystfactor<<"*weight_eletrigEff_up";
  else if (syst=="TRIGDOWN" && (channel=="ee" || channel=="enu")) mcweightsystfactor<<"*weight_eletrigEff_down";
  else if (channel=="ee" || channel=="enu") mcweightsystfactor<<"*weight_eletrigEff";

  if(channel=="taunu"||channel=="gamma"||channel=="nunu"||channel.find("qcd")!=channel.npos){
    if (syst=="LEPEFF_ELEUP") mcweightsystfactor<<"*weight_eleVeto_up/weight_eleVeto";
    if (syst=="LEPEFF_ELEDOWN") mcweightsystfactor<<"*weight_eleVeto_down/weight_eleVeto";
    if (syst=="LEPEFF_GSFUP") mcweightsystfactor<<"*weight_eleVeto_gsfup/weight_eleVeto";
    if (syst=="LEPEFF_GSFDOWN") mcweightsystfactor<<"*weight_eleVeto_gsfdown/weight_eleVeto";
    //HACK
    if (do_new_muVeto) {
      dataextrasel+="&&nvetomuons==0";
      if      (syst=="LEPEFF_MUIDUP")    mcweightsystfactor<<"*weight_0muloose_idup";
      else if (syst=="LEPEFF_MUIDDOWN")  mcweightsystfactor<<"*weight_0muloose_iddown";
      else if (syst=="LEPEFF_MUISOUP")   mcweightsystfactor<<"*weight_0muloose_isoup";
      else if (syst=="LEPEFF_MUISODOWN") mcweightsystfactor<<"*weight_0muloose_isodown";
      else if (syst=="LEPEFF_MUTKUP")    mcweightsystfactor<<"*weight_0muloose_tkup";
      else if (syst=="LEPEFF_MUTKDOWN")  mcweightsystfactor<<"*weight_0muloose_tkdown";
      else                               mcweightsystfactor<<"*weight_0muloose";
    } else {
      if (syst=="LEPEFF_MUIDUP") mcweightsystfactor<<"*weight_muVeto_idup/weight_muVeto";
      if (syst=="LEPEFF_MUIDDOWN") mcweightsystfactor<<"*weight_muVeto_iddown/weight_muVeto";
      if (syst=="LEPEFF_MUISOUP") mcweightsystfactor<<"*weight_muVeto_isoup/weight_muVeto";
      if (syst=="LEPEFF_MUISODOWN") mcweightsystfactor<<"*weight_muVeto_isodown/weight_muVeto";
      if (syst=="LEPEFF_MUTKUP") mcweightsystfactor<<"*weight_muVeto_tkup/weight_muVeto";
      if (syst=="LEPEFF_MUTKDOWN") mcweightsystfactor<<"*weight_muVeto_tkdown/weight_muVeto";

    }
  }
  else if (channel=="ee" || channel=="enu") {
    if (syst=="LEPEFF_ELEUP") mcweightsystfactor<<"*weight_eleTight_up/weight_eleTight";
    if (syst=="LEPEFF_ELEDOWN") mcweightsystfactor<<"*weight_eleTight_down/weight_eleTight";
    if (syst=="LEPEFF_GSFUP") mcweightsystfactor<<"*weight_gsfTight_up/weight_gsfTight";
    if (syst=="LEPEFF_GSFDOWN") mcweightsystfactor<<"*weight_gsfTight_down/weight_gsfTight";
  }
  else if (channel=="mumu" || channel=="munu") {
    if (syst=="LEPEFF_MUIDUP") mcweightsystfactor<<"*weight_muidTight_up/weight_muidTight";
    if (syst=="LEPEFF_MUIDDOWN") mcweightsystfactor<<"*weight_muidTight_down/weight_muidTight";
    if (syst=="LEPEFF_MUISOUP") mcweightsystfactor<<"*weight_muisoTight_up/weight_muisoTight";
    if (syst=="LEPEFF_MUISODOWN") mcweightsystfactor<<"*weight_muisoTight_down/weight_muisoTight";
    if (syst=="LEPEFF_MUTKUP") mcweightsystfactor<<"*weight_mutkTight_up/weight_mutkTight";
    if (syst=="LEPEFF_MUTKDOWN") mcweightsystfactor<<"*weight_mutkTight_down/weight_mutkTight";
  }

  //if (syst=="TRIG0UP") mcweightsystfactor<<"*weight_trig_1/weight_trig_0";
  //if (syst=="TRIG0DOWN") mcweightsystfactor<<"*weight_trig_2/weight_trig_0";
  //if (syst=="TRIG1UP") mcweightsystfactor<<"*weight_trig_3/weight_trig_0";
  //if (syst=="TRIG1DOWN") mcweightsystfactor<<"*weight_trig_4/weight_trig_0";
  //if (syst=="TRIG2UP") mcweightsystfactor<<"*weight_trig_5/weight_trig_0";
  //if (syst=="TRIG2DOWN") mcweightsystfactor<<"*weight_trig_6/weight_trig_0";

  if(channel=="taunu"||channel=="gamma"||channel=="nunu"||channel=="qcdB") {
    //remove weight_lepveto, only mu part
    if (do_new_muVeto) sigmcweight="weight_nolepnotrig*weight_mettrig*weight_eleVeto"+mcweightsystfactor.str();//"total_weight_lepveto"+mcweightsystfactor.str();
    else               sigmcweight="total_weight_lepveto"+mcweightsystfactor.str();
  }
  //remove trigger weight for e channels which do not use signal trigger
  else if (channel=="qcdA" || channel=="qcdC"){
    dataweight="pass_alljethttrigger";
    sigmcweight="weight_nolepnotrig"+mcweightsystfactor.str();
  }
  else if (channel=="ee" || channel=="enu") sigmcweight="weight_leptight*weight_nolepnotrig"+mcweightsystfactor.str();
  else if (channel=="mumu") sigmcweight="weight_leptight*weight_nolepnotrig*weight_mettrig_zmm"+mcweightsystfactor.str();
  else sigmcweight="total_weight_leptight"+mcweightsystfactor.str();

  //lepton veto weight
//   std::string lepveto_weight;
//   if(channel=="nunu") lepveto_weight="*weight_lepveto";
//   else lepveto_weight="*1";

  //add NLO reweighting
  sigmcweight=sigmcweight+"*v_nlo_Reweight*ewk_v_nlo_Reweight";

  //veto taus reweighting
  std::string sigmcweight_wtau;
  sigmcweight_wtau=sigmcweight+tauvetoweight;

  if (channel=="ee" || channel == "enu") dataset="SingleElectron";



  std::string bothcentral="TMath::Abs(jet1_eta)<3&&TMath::Abs(jet2_eta)<3";
  std::string bothforward="TMath::Abs(jet1_eta)>=3&&TMath::Abs(jet2_eta)>=3";
  std::string j2forwardj1central="TMath::Abs(jet1_eta)<3&&TMath::Abs(jet2_eta)>=3";
  std::string j1forwardj2central="TMath::Abs(jet1_eta)>=3&&TMath::Abs(jet2_eta)<3";

  std::string additionalcut=(syst=="PUUP")?("&&abs(puweight_up_scale)<200"): (syst=="PUDOWN")?("&&abs(puweight_down_scale)<10") : ("");
  if (syst.find("TRIG")!=syst.npos && channel!= "ee" && channel != "enu"){
    if (channel=="mumu") additionalcut="&&weight_mettrig_zmm>0";
    else additionalcut="&&weight_mettrig>0";
  }
//   if (syst.find("LEPEFF_MUIDUP")!=syst.npos && channel== "nunu") additionalcut="&&weight_muVeto_idup>0";
  analysis->set_baseselection(basesel+additionalcut);

  //DATA SHAPE GENERATION
  DataShape data("data");
  data.set_dataset(dataset)
    .set_dirname("data_obs")
    .set_shape(shape)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+dataextrasel);
  if (channel=="qcdA" || channel=="qcdC") {
    data.set_dataweight(dataweight);
  }

  DataShape totsignal125("totsignal125");
  totsignal125.set_dataset("H125")
    .set_dirname("sig125")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape qqH110("qqH110");
  qqH110.set_dataset("VBFH110")
    .set_dirname("qqH110")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape qqH125("qqH125");
  qqH125.set_dataset("VBFH125")
    .set_dirname("qqH125")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape qqH150("qqH150");
  qqH150.set_dataset("VBFH150")
    .set_dirname("qqH150")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape qqH200("qqH200");
  qqH200.set_dataset("VBFH200")
    .set_dirname("qqH200")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape qqH300("qqH300");
  qqH300.set_dataset("VBFH300")
    .set_dirname("qqH300")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape qqH400("qqH400");
  qqH400.set_dataset("VBFH400")
    .set_dirname("qqH400")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape qqH500("qqH500");
  qqH500.set_dataset("VBFH500")
    .set_dirname("qqH500")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape qqH600("qqH600");
  qqH600.set_dataset("VBFH600")
    .set_dirname("qqH600")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape qqH800("qqH800");
  qqH800.set_dataset("VBFH800")
    .set_dirname("qqH800")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape qqH1000("qqH1000");
  qqH1000.set_dataset("VBFH1000")
    .set_dirname("qqH1000")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape ggH110("ggH110");
  ggH110.set_dataset("GluGluH110")
    .set_dirname("ggH110")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape ggH125("ggH125");
  ggH125.set_dataset("GluGluH125")
    .set_dirname("ggH125")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape ggH150("ggH150");
  ggH150.set_dataset("GluGluH150")
    .set_dirname("ggH150")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape ggH200("ggH200");
  ggH200.set_dataset("GluGluH200")
    .set_dirname("ggH200")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape ggH300("ggH300");
  ggH300.set_dataset("GluGluH300")
    .set_dirname("ggH300")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape ggH400("ggH400");
  ggH400.set_dataset("GluGluH400")
    .set_dirname("ggH400")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape ggH500("ggH500");
  ggH500.set_dataset("GluGluH500")
    .set_dirname("ggH500")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape ggH600("ggH600");
  ggH600.set_dataset("GluGluH600")
    .set_dirname("ggH600")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape ggH800("ggH800");
  ggH800.set_dataset("GluGluH800")
    .set_dirname("ggH800")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape ggH1000("ggH1000");
  ggH1000.set_dataset("GluGluH1000")
    .set_dirname("ggH1000")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape vv("vv");
  vv.set_dataset("VV")
    .set_dirname("vv")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  DataShape topraw("topraw");
  topraw.set_dataset("Top")
    .set_dirname("top")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

//   DataShape znunuraw("znunuraw");
//   znunuraw.set_dataset("ZJets_nunu")
//     .set_dirname("zvv")
//     .set_shape(shape)
//     .set_dataweight(sigmcweight)
//     .set_basesel(analysis->baseselection())
//     .set_cat(sigcat+mcextrasel);


  DataShape qcdznunuraw("qcdznunuraw");
  qcdznunuraw.set_dataset("ZJets_nunu")
    .set_dirname("zvvqcd")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);


  DataShape ewkznunuraw("ewkznunuraw");
  ewkznunuraw.set_dataset("EWK_ZJets_nunu")
    .set_dirname("zvvewk")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

//   DataShape zmumuraw("zmumuraw");
//   zmumuraw.set_dataset("ZJets_ll")
//     .set_dirname("zmumu")
//     .set_shape(shape)
//     .set_dataweight(sigmcweight)
//     .set_basesel(analysis->baseselection())
//     .set_cat(sigcat+mcextrasel);
//   if (use_nlo)  zmumuraw.set_dataset("ZJets_ll_nlo");

  DataShape qcdzmumuraw("qcdzmumuraw");
  qcdzmumuraw.set_dataset("ZJets_ll")
    .set_dirname("zmumuqcd")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);
  if (use_nlo)  qcdzmumuraw.set_dataset("ZJets_ll_nlo");

  DataShape ewkzmumuraw("ewkzmumuraw");
  ewkzmumuraw.set_dataset("EWK_ZJets_ll")
    .set_dirname("zmumuewk")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

//   DataShape zeeraw("zeeraw");
//   zeeraw.set_dataset("ZJets_ll")
//     .set_dirname("zee")
//     .set_shape(shape)
//     .set_dataweight(sigmcweight)
//     .set_basesel(analysis->baseselection())
//     .set_cat(sigcat+mcextrasel);
//   if (use_nlo)  zeeraw.set_dataset("ZJets_ll_nlo");

  DataShape qcdzeeraw("qcdzeeraw");
  qcdzeeraw.set_dataset("ZJets_ll")
    .set_dirname("zeeqcd")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);
  if (use_nlo)  qcdzeeraw.set_dataset("ZJets_ll_nlo");

  DataShape ewkzeeraw("ewkzeeraw");
  ewkzeeraw.set_dataset("EWK_ZJets_ll")
    .set_dirname("zeeewk")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

//   DataShape zllraw("zllraw");
//   zllraw.set_dataset("ZJets_ll")
//     .set_dirname("zll")
//     .set_shape(shape)
//     .set_dataweight(sigmcweight)
//     .set_basesel(analysis->baseselection())
//     .set_cat(sigcat+mcextrasel);
//     if (use_nlo)  zllraw.set_dataset("ZJets_ll_nlo");

  DataShape qcdzllraw("qcdzllraw");
  qcdzllraw.set_dataset("ZJets_ll")
    .set_dirname("zllqcd")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);
    if (use_nlo)  qcdzllraw.set_dataset("ZJets_ll_nlo");

  DataShape ewkzllraw("ewkzllraw");
  ewkzllraw.set_dataset("EWK_ZJets_ll")
    .set_dirname("zllewk")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

//   DataShape qcdwraw("qcdwraw");
//   qcdwraw.set_dataset("QCDWJets")
//     .set_dirname("wqcd")
//     .set_shape(shape)
//     .set_dataweight(sigmcweight)
//     .set_basesel(analysis->baseselection())
//     .set_cat(sigcat+mcextrasel);
//   if (use_nlo)  qcdwraw.set_dataset("QCDWJets_nlo");
// 
//   DataShape ewkwraw("ewkwraw");
//   ewkwraw.set_dataset("EWKWJets")
//     .set_dirname("wewk")
//     .set_shape(shape)
//     .set_dataweight(sigmcweight)
//     .set_basesel(analysis->baseselection())
//     .set_cat(sigcat+mcextrasel);

//   DataShape wmunuraw("wmunuraw");
//   wmunuraw.set_dataset("WJets_munu")
//     .set_dirname("wmu")
//     .set_shape(shape)
//     .set_dataweight(sigmcweight) //+lepveto_weight)
//     .set_basesel(analysis->baseselection())
//     .set_cat(sigcat+mcextrasel);
//   if (use_nlo)  wmunuraw.set_dataset("WJets_nlo_munu");

  DataShape qcdwmunuraw("qcdwmunuraw");
  qcdwmunuraw.set_dataset("WJets_munu")
    .set_dirname("wmuqcd")
    .set_shape(shape)
    .set_dataweight(sigmcweight) //+lepveto_weight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);
  if (use_nlo)  qcdwmunuraw.set_dataset("WJets_nlo_munu");

  DataShape ewkwmunuraw("ewkwmunuraw");
  ewkwmunuraw.set_dataset("EWK_WJets_munu")
    .set_dirname("wmuewk")
    .set_shape(shape)
    .set_dataweight(sigmcweight) //+lepveto_weight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

//   DataShape wenuraw("wenuraw");
//   wenuraw.set_dataset("WJets_enu")
//     .set_dirname("wel")
//     .set_shape(shape)
//     .set_dataweight(sigmcweight) //+lepveto_weight)
//     .set_basesel(analysis->baseselection())
//     .set_cat(sigcat+mcextrasel);
//   if (use_nlo)  wenuraw.set_dataset("WJets_nlo_enu");

  DataShape qcdwenuraw("qcdwenuraw");
  qcdwenuraw.set_dataset("WJets_enu")
    .set_dirname("welqcd")
    .set_shape(shape)
    .set_dataweight(sigmcweight) //+lepveto_weight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);
  if (use_nlo)  qcdwenuraw.set_dataset("WJets_nlo_enu");

  DataShape ewkwenuraw("ewkwenuraw");
  ewkwenuraw.set_dataset("EWK_WJets_enu")
    .set_dirname("welewk")
    .set_shape(shape)
    .set_dataweight(sigmcweight) //+lepveto_weight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

//   DataShape wtaunuraw("wtaunuraw");
//   wtaunuraw.set_dataset("WJets_taunu")
//     .set_dirname("wtau")
//     .set_shape(shape)
//     .set_dataweight(sigmcweight_wtau)
//     .set_basesel(analysis->baseselection())
//     .set_cat(sigcat_novetotaus+mcextrasel);
//   if (use_nlo) wtaunuraw.set_dataset("WJets_nlo_taunu");

  DataShape qcdwtaunuraw("qcdwtaunuraw");
  qcdwtaunuraw.set_dataset("WJets_taunu")
    .set_dirname("wtauqcd")
    .set_shape(shape)
    .set_dataweight(sigmcweight_wtau)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat_novetotaus+mcextrasel);
  if (use_nlo) qcdwtaunuraw.set_dataset("WJets_nlo_taunu");

  DataShape ewkwtaunuraw("ewkwtaunuraw");
  ewkwtaunuraw.set_dataset("EWK_WJets_taunu")
    .set_dirname("wtauewk")
    .set_shape(shape)
    .set_dataweight(sigmcweight_wtau)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat_novetotaus+mcextrasel);

  DataShape qcdraw("qcdraw");
  qcdraw.set_dataset("QCD")
    .set_dirname("qcd")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  //if (channel=="qcdA" || channel=="qcdC") qcdraw.set_cat(sigcat+mcextrasel+dataextrasel);


  DataShape gjetsraw("gjetsraw");
  gjetsraw.set_dataset("GJets")
    .set_dirname("gjets")
    .set_shape(shape)
    .set_dataweight(sigmcweight)
    .set_basesel(analysis->baseselection())
    .set_cat(sigcat+mcextrasel);

  //HISTPLOTTER
  std::vector<LTShapeElement> shapevec;
  for(unsigned ishape=0;ishape<shape.size();ishape++){
    //do not plot 2D hists....
    if (shape[ishape].find(":")!=shape[ishape].npos && shape[ishape].find("::")==shape[ishape].npos) continue;

    std::string strs = extractShapeName(shape[ishape]);

    LTShapeElement thisshape;
    thisshape.set_name(strs);
    thisshape.set_legleft(0.65);
    thisshape.set_legright(0.97);


//     if (strs=="forward_tag_eta"||strs=="central_tag_eta"){
//       thisshape.set_legleft(0.39);
//       thisshape.set_legright(0.71);
//     }
    if (strs.find("fourjetsmetnomu")!=strs.npos) thisshape.set_axisrangemultiplier((channel=="mumu"||channel=="nunu")?2.2:1.6);
    if (strs.find("dijet_deta")!=strs.npos && channel=="nunu") thisshape.set_axisrangemultiplier(1.7);

    thisshape.set_histtitle(histTitle[ishape]);
    if(do_logy) thisshape.set_dology(true);
    shapevec.push_back(thisshape);
  }

  std::vector<LTPlotElement> elementvec;

  LTPlotElement dataele;
  dataele.set_is_data(true)
    .set_scale(1)
    .set_legname("Data")
    .set_is_inrationum(true)
    .set_sample("data_obs");
  if(runblindreg&&(channel=="nunu"||channel=="qcdD")){
    std::vector<std::string> blindvars;
    std::vector<std::pair<double,double> > blindrange;
    blindvars.push_back("fourjetsmetnomu_mindphi");
    blindvars.push_back("metnomu_significance");
    blindvars.push_back("metnomuons");
    blindvars.push_back("dijet_deta");
    blindvars.push_back("metnomuons*pow");
    std::pair<double,double>mindphirange(1.,10.);
    std::pair<double,double>metsigrange(150.,1000.);
    std::pair<double,double>metnomurange(350.,1000.);
    std::pair<double,double>detarange(5.,10.);
    blindrange.push_back(mindphirange);
    blindrange.push_back(metsigrange);
    blindrange.push_back(metnomurange);
    blindrange.push_back(detarange);
    blindrange.push_back(std::pair<double,double>(3,100));
    dataele.set_blindvar(blindvars)
      .set_blindrange(blindrange);
  }

  if(blindcutreg&&!runblindreg){
    std::vector<std::string> blindvars;
    std::vector<std::pair<double,double> > blindrange;
    blindvars.push_back("forward_tag_eta");
    std::pair<double,double>forwardetarange(-1.8,1.8);
    blindrange.push_back(forwardetarange);
    dataele.set_blindvar(blindvars)
      .set_blindrange(blindrange);
  }

  //All W
//   LTPlotElement qcdwele;
//   qcdwele.set_is_data(false)
//     .set_scale(getPostFitSF(channel,"wqcd"))
//     .set_color(kOrange-4)
//     .set_in_stack(true)
//     .set_is_inratioden(true)
//     .set_legname("QCD W#rightarrow l#nu")
//     .set_sample("wqcd");
//   if(!do_mcbkg)qcdwele.set_has_dderrors(1);
// 
//   LTPlotElement ewkwele;
//   ewkwele.set_is_data(false)
//     .set_scale(getPostFitSF(channel,"wewk"))
//     .set_color(kOrange+2)
//     .set_in_stack(true)
//     .set_is_inratioden(true)
//     .set_legname("EWK W#rightarrow l#nu")
//     .set_sample("wewk");
//   if(!do_mcbkg)ewkwele.set_has_dderrors(1);

  //separate W
//   LTPlotElement wmunuele;
//   wmunuele.set_is_data(false)
//     .set_scale(1)
//     .set_color(kOrange-4)
//     .set_in_stack(true)
//     .set_is_inratioden(true)
//     .set_legname("W#rightarrow#mu#nu")
//     .set_sample("wmu");
//   if(!do_mcbkg)wmunuele.set_has_dderrors(1);

  LTPlotElement qcdwmunuele;
  qcdwmunuele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"wqcd"))
    .set_color(TColor::GetColor("#E19D07"))
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("W(#mu#nu)+jets QCD")
    .set_sample("wmuqcd");
  if(!do_mcbkg)qcdwmunuele.set_has_dderrors(1);

  LTPlotElement ewkwmunuele;
  ewkwmunuele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"wewk"))
    .set_color(kAzure+2)
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("W(#mu#nu)+jets EWK")
    .set_sample("wmuewk");
  if(!do_mcbkg)ewkwmunuele.set_has_dderrors(1);

//   LTPlotElement wenuele;
//   wenuele.set_is_data(false)
//     .set_scale(1)
//     .set_color(kOrange  + 2)
//     .set_in_stack(true)
//     .set_is_inratioden(true)
//     .set_legname("W#rightarrow e#nu")
//     .set_sample("wel");
//   if(!do_mcbkg)wenuele.set_has_dderrors(1);

  LTPlotElement qcdwenuele;
  qcdwenuele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"wqcd"))
    .set_color(TColor::GetColor("#FAAF08"))
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("W(e#nu)+jets QCD")
    .set_sample("welqcd");
  if(!do_mcbkg)qcdwenuele.set_has_dderrors(1);

  LTPlotElement ewkwenuele;
  ewkwenuele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"wewk"))
    .set_color(kAzure+1)
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("W(e#nu)+jets EWK")
    .set_sample("welewk");
  if(!do_mcbkg)ewkwenuele.set_has_dderrors(1);

//   LTPlotElement wtaunuele;
//   wtaunuele.set_is_data(false)
//     .set_scale(1)
//     .set_color(kOrange + 4)
//     .set_in_stack(true)
//     .set_is_inratioden(true)
//     .set_legname("W#rightarrow#tau#nu")
//     .set_sample("wtau");
//   if(!do_mcbkg)wtaunuele.set_has_dderrors(1);

  LTPlotElement qcdwtaunuele;
  qcdwtaunuele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"wqcd"))
    .set_color(TColor::GetColor("#C88C06"))
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("W(#tau#nu)+jets QCD")
    .set_sample("wtauqcd");
  if(!do_mcbkg)qcdwtaunuele.set_has_dderrors(1);

  LTPlotElement ewkwtaunuele;
  ewkwtaunuele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"wewk"))
    .set_color(kAzure+3)
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("W(#tau#nu)+jets EWK")
    .set_sample("wtauewk");
  if(!do_mcbkg)ewkwtaunuele.set_has_dderrors(1);

  //separate Z
//   LTPlotElement zmumuele;
//   zmumuele.set_is_data(false)
//     .set_scale(1)
//     .set_color(TColor::GetColor("#9A9EAB"))
//     .set_in_stack(true)
//     .set_is_inratioden(true)
//     .set_legname("Z#rightarrow#mu#mu")
//     .set_sample("zmumu");
//   if(!do_mcbkg)zmumuele.set_has_dderrors(1);

  LTPlotElement qcdzmumuele;
  qcdzmumuele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"zqcd"))
    .set_color(TColor::GetColor("#9A9EAB"))
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("Z(#mu#mu) QCD")
    .set_sample("zmumuqcd");
  if(!do_mcbkg)qcdzmumuele.set_has_dderrors(1);

  LTPlotElement ewkzmumuele;
  ewkzmumuele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"zewk"))
    .set_color(TColor::GetColor("#91ABC4"))
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("Z(#mu#mu) EWK")
    .set_sample("zmumuewk");
  if(!do_mcbkg)ewkzmumuele.set_has_dderrors(1);

//   LTPlotElement zeeele;
//   zeeele.set_is_data(false)
//     .set_scale(1)
//     .set_color(TColor::GetColor("#9A9EAB"))
//     .set_in_stack(true)
//     .set_is_inratioden(true)
//     .set_legname("Z#rightarrow ee")
//     .set_sample("zee");
//   if(!do_mcbkg)zeeele.set_has_dderrors(1);

  LTPlotElement qcdzeeele;
  qcdzeeele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"zqcd"))
    .set_color(TColor::GetColor("#9A9EAB"))
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("Z(ee) QCD")
    .set_sample("zeeqcd");
  if(!do_mcbkg)qcdzeeele.set_has_dderrors(1);

  LTPlotElement ewkzeeele;
  ewkzeeele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"zewk"))
    .set_color(TColor::GetColor("#91ABC4"))
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("Z(ee) EWK")
    .set_sample("zeeewk");
  if(!do_mcbkg)ewkzeeele.set_has_dderrors(1);

  //Zll for SR
//   LTPlotElement zllele;
//   zllele.set_is_data(false)
//   .set_scale(1)
//   .set_color(TColor::GetColor("#9A9EAB"))
//   .set_in_stack(true)
//   .set_is_inratioden(true)
//   .set_legname("Z#rightarrow ll")
//   .set_sample("zll");
//   if(!do_mcbkg)zllele.set_has_dderrors(1);

  LTPlotElement qcdzllele;
  qcdzllele.set_is_data(false)
  .set_scale(1)
  .set_color(TColor::GetColor("#9A9EAB"))
  .set_in_stack(true)
  .set_is_inratioden(true)
  .set_legname("Z(ll) QCD")
  .set_sample("zllqcd");
  if(!do_mcbkg)qcdzllele.set_has_dderrors(1);

  LTPlotElement ewkzllele;
  ewkzllele.set_is_data(false)
  .set_scale(1)
  .set_color(TColor::GetColor("#91ABC4"))
  .set_in_stack(true)
  .set_is_inratioden(true)
  .set_legname("Z(ll) EWK")
  .set_sample("zllewk");
  if(!do_mcbkg)ewkzllele.set_has_dderrors(1);

  //Znunu
//   LTPlotElement znunuele;
//   znunuele.set_is_data(false)
//     .set_scale(1)
//     .set_color(TColor::GetColor("#4D975D"))
//     .set_in_stack(true)
//     .set_is_inratioden(true)
//     .set_legname("Z(#nu#nu)+jets")
//     .set_sample("zvv");
//   if(!do_mcbkg)znunuele.set_has_dderrors(1);

  LTPlotElement qcdznunuele;
  qcdznunuele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"zqcd"))
    .set_color(TColor::GetColor("#4D975D"))
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("Z(#nu#nu)+jets QCD")
    .set_sample("zvvqcd");
  if(!do_mcbkg)qcdznunuele.set_has_dderrors(1);

  LTPlotElement ewkznunuele;
  ewkznunuele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"zewk"))
    .set_color(kCyan+1)
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("Z(#nu#nu)+jets EWK")
    .set_sample("zvvewk");
  if(!do_mcbkg)ewkznunuele.set_has_dderrors(1);

  LTPlotElement vvele;
  vvele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"VV"))
    .set_color(TColor::GetColor("#4897D8"))
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("VV")
    .set_sample("vv");

  LTPlotElement topele;
  topele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"top"))
    .set_color(TColor::GetColor("#CF3721"))
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("Top quark")
    .set_sample("top");

  LTPlotElement qcdele;
  qcdele.set_is_data(false)
    .set_scale(getPostFitSF(channel,"QCD"))
    //AMM uncomment for QCD in signal region plots mindphi.
    //.set_scale(3.8/7.9)
    //AMM uncomment for QCD in signal region plots all but mindphi.
    //.set_scale(3.8/6.5)
    .set_color(TColor::GetColor("#F1F1F2"))
    .set_in_stack(true)
    .set_is_inratioden(true)
    .set_legname("QCD")
    .set_sample("qcd");

//   LTPlotElement gjetsele;
//   gjetsele.set_is_data(false)
//     .set_scale(1)
//     .set_color(kYellow+2)
//     .set_in_stack(true)
//     .set_is_inratioden(true)
//     .set_legname("#gamma+jets")
//     .set_sample("gjets");


  LTPlotElement sigele;
  sigele.set_is_data(false)
    .set_scale(1)
    .set_color(kRed)
    .set_line_style(1)
    .set_in_stack(false)
    .set_legname("Signal")
    .set_sample("sig125");

  LTPlotElement qqHele;
  qqHele.set_is_data(false)
    .set_scale(1)
    .set_color(kBlack)
    .set_line_style(1)
    .set_in_stack(false)
    .set_legname("qqH_{inv} m_{h} = 125 GeV")
    .set_sample("qqH125");

  LTPlotElement ggHele;
  ggHele.set_is_data(false)
    .set_scale(1)
    .set_color(kBlack)
    .set_line_style(2)
    .set_in_stack(false)
    .set_legname("ggH_{inv} m_{h} = 125 GeV")
    .set_sample("ggH125");

  if(!((channel=="nunu"||channel=="qcdD")&&runblind))elementvec.push_back(dataele);
  if(!dataonly){
    elementvec.push_back(qcdele);

    if(channel=="mumu") {
      elementvec.push_back(ewkzmumuele);
      elementvec.push_back(qcdzmumuele);
    }
    else if(channel=="ee") {
      elementvec.push_back(ewkzeeele);
      elementvec.push_back(qcdzeeele);
    }
    if (channel=="nunu" || channel=="enu" || channel=="munu"|| channel=="taunu"){
      elementvec.push_back(ewkzllele);
      elementvec.push_back(qcdzllele);
    }

    elementvec.push_back(topele);
    elementvec.push_back(vvele);

//     if (channel=="gamma"){
//       elementvec.push_back(gjetsele);
//     }

    elementvec.push_back(ewkwtaunuele);
    elementvec.push_back(ewkwmunuele);
    elementvec.push_back(ewkwenuele);

    if(channel=="nunu" || channel=="qcd" || channel=="taunu") {
      elementvec.push_back(ewkznunuele);
    }

    elementvec.push_back(qcdwtaunuele);
    elementvec.push_back(qcdwmunuele);
    elementvec.push_back(qcdwenuele);

    if(channel=="nunu" || channel=="qcd" || channel=="taunu") {
      elementvec.push_back(qcdznunuele);
    }
    if(channel=="nunu" || channel=="qcd" || channel=="taunu") {
//       elementvec.push_back(sigele);
      if (channel=="nunu"){
        elementvec.push_back(qqHele);
        elementvec.push_back(ggHele);
      }
    }
  }

  if (debug){
    std::cout << " number of elements: " << elementvec.size() << std::endl;
    for(unsigned iElement=0;iElement<elementvec.size();iElement++){
      std::cout<<" --- Ele " << iElement << " sample " << elementvec[iElement].sample() << std::endl;
    }
  }

  HistPlotter plotter("plotter");
  plotter.set_dirname("ControlPlots")
    .set_add_underflows(true)
    .set_add_overflows(true)
    .set_elements(elementvec)
    //.set_histTitles(histTitle)
    .set_shapes(shapevec)
    .set_toterror(getPostFitSF(channel,"error"));
  if(!dataonly) plotter.set_do_ratio(true);
  else plotter.set_do_ratio(false);
  if(channel=="nunu"&&runblind)plotter.set_do_ratio(false);

  if (debug) plotter.set_do_debug(true);

  std::vector<std::string> dirvec;
  if(!dataonly){
    dirvec.push_back("welqcd");
    dirvec.push_back("wmuqcd");
    dirvec.push_back("wtauqcd");
    dirvec.push_back("welewk");
    dirvec.push_back("wmuewk");
    dirvec.push_back("wtauewk");
    dirvec.push_back("zmumuqcd");
    dirvec.push_back("zmumuewk");
    dirvec.push_back("zeeqcd");
    dirvec.push_back("zeeewk");
    dirvec.push_back("zllqcd");
    dirvec.push_back("zllewk");
//     dirvec.push_back("wqcd");
//     dirvec.push_back("wewk");
    dirvec.push_back("zvvqcd");
    dirvec.push_back("zvvewk");
    dirvec.push_back("gjets");
    dirvec.push_back("qcd");
    dirvec.push_back("vv");
    //dirvec.push_back("wg");  
    dirvec.push_back("top");
    dirvec.push_back("sig125");
    dirvec.push_back("qqH125");
    dirvec.push_back("ggH125");
  }
  if(!(channel=="nunu"&&runblind))dirvec.push_back("data_obs");

  SummaryTable summary("summary");
  summary.set_shape(shapevec)
    .set_dirs(dirvec);
  
  /*##########################################
  #                                          #
  #   SET UP ANALYSIS SEQUENCE AND RUN       #
  #                                          #
  ##########################################*/

  if(!dataonly){

    if(do_mcbkg){
      analysis->AddModule(&qcdraw);
      analysis->AddModule(&qcdwmunuraw);
      analysis->AddModule(&qcdwenuraw);
      analysis->AddModule(&qcdwtaunuraw);
      analysis->AddModule(&ewkwmunuraw);
      analysis->AddModule(&ewkwenuraw);
      analysis->AddModule(&ewkwtaunuraw);
//       analysis->AddModule(&qcdwraw);
//       analysis->AddModule(&ewkwraw);
      analysis->AddModule(&topraw);
      analysis->AddModule(&gjetsraw);
      analysis->AddModule(&qcdzmumuraw);
      analysis->AddModule(&ewkzmumuraw);
      analysis->AddModule(&qcdzeeraw);
      analysis->AddModule(&ewkzeeraw);
      analysis->AddModule(&qcdzllraw);
      analysis->AddModule(&ewkzllraw);
      analysis->AddModule(&qcdznunuraw);
      analysis->AddModule(&ewkznunuraw);
    }
    analysis->AddModule(&vv);
  }
  if(!(channel=="nunu"&&runblind))analysis->AddModule(&data);

  if(!dataonly){
    analysis->AddModule(&ggH110);
    analysis->AddModule(&ggH150);
    analysis->AddModule(&ggH200);
    analysis->AddModule(&ggH300);
    analysis->AddModule(&ggH400);
    analysis->AddModule(&ggH500);
    analysis->AddModule(&ggH600);
    analysis->AddModule(&ggH800);
    analysis->AddModule(&ggH1000);
    analysis->AddModule(&qqH110);
    analysis->AddModule(&qqH150);
    analysis->AddModule(&qqH200);
    analysis->AddModule(&qqH300);
    analysis->AddModule(&qqH400);
    analysis->AddModule(&qqH500);
    analysis->AddModule(&qqH600);
    analysis->AddModule(&qqH800);
    analysis->AddModule(&qqH1000);
    analysis->AddModule(&qqH125);
    analysis->AddModule(&ggH125);
    analysis->AddModule(&totsignal125);
  }
  analysis->AddModule(&plotter);
  if(!dataonly) analysis->AddModule(&summary);

  analysis->RunAnalysis();

  return 0;

}

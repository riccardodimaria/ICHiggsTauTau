#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/LightTreeRDM.h"
#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"
#include "UserCode/ICHiggsTauTau/interface/TriggerPath.hh"
#include "UserCode/ICHiggsTauTau/interface/TriggerObject.hh"
#include "UserCode/ICHiggsTauTau/interface/EventInfo.hh"
#include "UserCode/ICHiggsTauTau/interface/Vertex.hh"
#include "UserCode/ICHiggsTauTau/interface/city.h"
#include "UserCode/ICHiggsTauTau/interface/L1TObject.hh"
#include "TVector3.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/HTFromLHEParticles.h"


namespace ic {

  LightTreeRDM::LightTreeRDM(std::string const& name): ModuleBase(name){
    fs_ = NULL;
    debug_ = 0;
    met_label_ = "pfMetType1Collection";
    jet_label_ = "pfJetsPFlow";
    dijet_label_ = "jjCandidates";
    is_data_ = false;
    do_trigskim_ = false;
    do_promptskim_ = false;
    do_noskim_ = false;
    ignoreLeptons_=false;

    countZmumu_ = 0;
    countZee_ = 0;

    outputTree_ = 0;

    //set vectors, do not reset later!!!
    //Used for setting branch addresses...
    nJetsSave_ = 6;


    //initialise arrays fixed size
    jet_pt_ = new double[nJetsSave_];
    jet_uncorpt_ = new double[nJetsSave_];
    jet_E_ = new double[nJetsSave_];
    jet_eta_ = new double[nJetsSave_];
    jet_phi_ = new double[nJetsSave_];
    jet_csv_ = new double[nJetsSave_];
    jet_jetid_ = new double[nJetsSave_];
    jet_loosepuid_ = new double[nJetsSave_];
    jet_tightpuid_ = new double[nJetsSave_];
    jet_flavour_ = new int[nJetsSave_];
    jet_neutralhadfrac_ = new double[nJetsSave_];
    jet_neutralemfrac_ = new double[nJetsSave_];
    jet_chargedemfrac_ = new double[nJetsSave_];
    jet_chargedhadfrac_ = new double[nJetsSave_];
    jet_chargedmult_ = new unsigned[nJetsSave_];
    jet_neutralmult_ = new unsigned[nJetsSave_];

    jet_genjet_mindR_ = new double[nJetsSave_];
    jet_genid_ = new unsigned[nJetsSave_];
    jet_genpt_ = new double[nJetsSave_];
    jet_geneta_ = new double[nJetsSave_];
    jet_genphi_ = new double[nJetsSave_];
    /*
    for (unsigned ij(0); ij<nJetsSave_;++ij){
      jet_pt_.push_back(0);
      jet_E_.push_back(0);
      jet_eta_.push_back(-5);
      jet_phi_.push_back(-5);
      jet_csv_.push_back(-1);
      jet_genjet_mindR_.push_back(99);
      jet_flavour_.push_back(0);
      jet_jetid_.push_back(-1);
      jet_puid_.push_back(-1);
      jet_genid_.push_back(1000);
      jet_genpt_.push_back(0);
      jet_geneta_.push_back(-5);
      jet_genphi_.push_back(-5);
      }*/
    resetAllTreeVariables();
  }

  void LightTreeRDM::resetAllTreeVariables(){

    run_=0;
    lumi_=0;
    event_=0;
    weight_nolepnotrig_=1;
    for (unsigned iT(0); iT<7; ++iT){
      weight_trig_[iT]=1;
    }
    weight_lepveto_ = 1;
    weight_leptight_ = 1;
    total_weight_lepveto_ = 1;
    total_weight_leptight_ = 1;
    puweight_up_scale_=1;
    puweight_down_scale_=1;
    v_nlo_Reweight_=1;

    nJets_ = 0;
    nGenJets_ = 0;


    mynewloosephotons_.clear();
    mynewmediumphotons_.clear();
    mynewtightphotons_.clear();
    mynewvetoelectrons_.clear();
    mynewselelectrons_.clear();
    mynewvetomuons_.clear();
    mynewselmuons_.clear();




    for (unsigned ij(0); ij<nJetsSave_;++ij){
      jet_pt_[ij] = -1;
      jet_uncorpt_[ij] = -1;
      jet_E_[ij] = -1;
      jet_eta_[ij] = -5;
      jet_phi_[ij] = -5;
      jet_csv_[ij] = -1;
      jet_jetid_[ij] = -1;
      jet_loosepuid_[ij] = -1;
      jet_tightpuid_[ij] = -1;
      jet_flavour_[ij] = 0;
      jet_neutralhadfrac_[ij] = -1;
      jet_neutralemfrac_[ij] = -1;
      jet_chargedemfrac_[ij] = -1;
      jet_chargedhadfrac_[ij] = -1;
      jet_chargedmult_[ij] = 0;
      jet_neutralmult_[ij] = 0;
      jet_genjet_mindR_[ij] = 99;
      jet_genid_[ij] = 1000;
      jet_genpt_[ij] = -1;
      jet_geneta_[ij] = -5;
      jet_genphi_[ij] = -5;
    }

    n_jets_cjv_30_ = 0;
    n_jets_cjv_20EB_30EE_ = 0;
    n_jets_15_ = 0;
    n_jets_30_ = 0;
    n_jets_csv2medium_ = 0;
    cjvjetpt_=-1;

    dijet_M_ = 0;
    dijet_deta_ = 0;
    dijet_sumeta_ = 0;
    dijet_dphi_ = 0;
    forward_tag_eta_ = -5;
    central_tag_eta_ = -5;

    sumet_ = 0;

    l1met_ = 0;
    l1mht_ = 0;
    met_ = 0;
    calomet_ = 0;
    genmet_ = 0;
    genmetphi_ = 0;
    metnomuons_ =0;
    metnoelectrons_ =0;
    met_x_ = 0;
    met_y_ = 0;
    metnomu_x_ = 0;
    metnomu_y_ = 0;
    metnoel_x_ = 0;
    metnoel_y_ = 0;
    met_significance_ = 0;
    metnomu_significance_ = 0;
    metnoel_significance_ = 0;

    ht_ = 0;
    ht30_ = 0;
    mht_ = 0;
    sqrt_ht_ = 0;

    unclustered_et_ = 0;
    jetunclet_mindphi_ = 0;
    metnomuunclet_dphi_ = 0;
    metnoelunclet_dphi_ = 0;

    alljetsmetnomu_mindphi_ = 0;
    alljetsnotaumetnomu_mindphi_ = 0;
    jetmetnomu_mindphi_ = 0;

    alljetsmetnoel_mindphi_ = 0;
    alljetsnotaumetnoel_mindphi_ = 0;
    jetmetnoel_mindphi_ = 0;

    dijetmetnomu_scalarSum_pt_ = 0;
    dijetmetnomu_vectorialSum_pt_ = 0;
    dijetmetnomu_ptfraction_ = 0;
    jet1metnomu_scalarprod_ = 0;
    jet2metnomu_scalarprod_ = 0;

    dijetmetnoel_scalarSum_pt_ = 0;
    dijetmetnoel_vectorialSum_pt_ = 0;
    dijetmetnoel_ptfraction_ = 0;
    jet1metnoel_scalarprod_ = 0;
    jet2metnoel_scalarprod_ = 0;

    pass_muontrigger_ = -1;
    pass_photontrigger_ = -1;
    pass_sigtrigger_ = -1;
    pass_mettrigger_ = -1;
    pass_metmht90trigger_ = -1;
    pass_metmht100trigger_ = -1;
    pass_metmht110trigger_ = -1;
    pass_metmht120trigger_ = -1;
    pass_controltrigger_ = -1;
    pass_singleEltrigger_ = -1;

    nvetomuons_=0;
    nselmuons_=0;
    nvetoelectrons_=0;
    nselelectrons_=0;
    nvetotaus_=0;
    ntaus_=0;
    nloosephotons_=0;
    nmediumphotons_=0;
    ntightphotons_=0;

    m_mumu_=-1;
    pt_mumu_=-1;
    oppsign_mumu_=false;
    m_ee_=-1;
    pt_ee_=-1;
    oppsign_ee_=false;
    m_mumu_gen_=-1;
    m_ee_gen_=-1;
    lep_mt_=-1;

// Define gen_mu1_*

    gen_mu1_pt_=-1;
    gen_mu1_eta_=-5;
    gen_mu1_phi_=-5;
    gen_mu1_mindR_j1_=99;
    gen_mu1_mindR_j2_=99;

    mu1_pt_=-1;
    mu1_eta_=-5;
    mu1_phi_=-5;
    mu1_isTight_=false;
    mu1_genmindR_=99;
    mu1_genpt_=-1;
    mu1_geneta_=-5;
    mu1_genphi_=-5;
    mu2_pt_=-1;
    mu2_eta_=-5;
    mu2_phi_=-5;
    mu2_isTight_=false;
    mu2_genmindR_=99;
    mu2_genpt_=-1;
    mu2_geneta_=-5;
    mu2_genphi_=-5;

// Define gen_ele1_*

    gen_ele1_pt_=-1;
    gen_ele1_eta_=-5;
    gen_ele1_phi_=-5;
    gen_ele1_mindR_j1_=99;
    gen_ele1_mindR_j2_=99;


    ele1_pt_=-1;
    ele1_eta_=-5;
    ele1_phi_=-5;
    ele1_isTight_=false;
    ele1_genmindR_=99;
    ele1_genpt_=-1;
    ele1_geneta_=-5;
    ele1_genphi_=-5;

    ele2_pt_=-1;
    ele2_eta_=-5;
    ele2_phi_=-5;
    ele2_isTight_=false;
    ele2_genmindR_=99;
    ele2_genpt_=-1;
    ele2_geneta_=-5;
    ele2_genphi_=-5;

    tau1_pt_=-1;
    tau1_eta_=-5;
    tau1_phi_=-5;
    tau1_isTight_=false;
    tau1_genmindR_=99;
    tau1_genpt_=-1;
    tau1_geneta_=-5;
    tau1_genphi_=-5;

    vetotau1_pt_=-1;
    vetotau1_eta_=-5;
    vetotau1_phi_=-5;
    vetotau1_isTight_=false;
    vetotau1_genmindR_=99;
    vetotau1_genpt_=-1;
    vetotau1_geneta_=-5;
    vetotau1_genphi_=-5;

    gamma1_pt_=-1;
    gamma1_eta_=-5;
    gamma1_phi_=-5;
    gamma1_genmindR_=99;
    gamma1_genpt_=-1;
    gamma1_geneta_=-5;
    gamma1_genphi_=-5;

    n_vertices_=0;


    // lheParticles
    lheHT_=0;

    // Boson pt
    boson_pt_=-50;
  }

  LightTreeRDM::~LightTreeRDM(){
    ;
  }


  int  LightTreeRDM::PreAnalysis(){
    std::cout << "-----------------------------------------------" << std::endl
              << "------ PreAnalysis Info for LightTreeRDM ------" << std::endl
              << "-----------------------------------------------" << std::endl;
    if (fs_) {
      std::cout << "MET Label: " << met_label_ << std::endl;
      std::cout << "Jets Label: " << jet_label_ << std::endl;
      std::cout << "dijet Label: " << dijet_label_ << std::endl;
      if (is_data_) std::cout << "Processing set for data !" << std::endl;
      else  std::cout << "Processing set for MC !" << std::endl;
      std::cout << "ignoreLeptons = " << ignoreLeptons_ << std::endl
                << "do_trigskim = " << do_trigskim_ << std::endl
                << "do_promptskim = " << do_promptskim_ << std::endl
                << "do_noskim = "  << do_noskim_ << std::endl;
    }
    fs_->cd();
    outputTree_=fs_->make<TTree>("LightTree","Tree containing LightTreeAna input variables");
//  TTree *     outputTree_ = new TTree("LightTree","Tree containing LightTreeAna input variables"); 
    outputTree_->Branch("run",&run_);
    outputTree_->Branch("lumi",&lumi_);
    outputTree_->Branch("event",&event_);
    outputTree_->Branch("weight_nolepnotrig",&weight_nolepnotrig_);
    for (unsigned iT(0); iT<7; ++iT){
      std::ostringstream label;
      label << "weight_trig_" << iT;
      outputTree_->Branch(label.str().c_str(),&weight_trig_[iT]);
    }
    outputTree_->Branch("weight_lepveto",&weight_lepveto_);
    outputTree_->Branch("weight_leptight",&weight_leptight_);
    outputTree_->Branch("total_weight_lepveto",&total_weight_lepveto_);
    outputTree_->Branch("total_weight_leptight",&total_weight_leptight_);
    outputTree_->Branch("puweight_up_scale",&puweight_up_scale_);
    outputTree_->Branch("puweight_down_scale",&puweight_down_scale_);
    outputTree_->Branch("v_nlo_Reweight",&v_nlo_Reweight_);

    outputTree_->Branch("nGenJets",&nGenJets_);

    for (unsigned ij(0); ij<nJetsSave_;++ij){
      std::ostringstream label;
      label << "jet" << ij+1;
      outputTree_->Branch((label.str()+"_pt").c_str(),&jet_pt_[ij]);
      outputTree_->Branch((label.str()+"_uncorpt").c_str(),&jet_uncorpt_[ij]);
      outputTree_->Branch((label.str()+"_E").c_str(),&jet_E_[ij]);
      outputTree_->Branch((label.str()+"_eta").c_str(),&jet_eta_[ij]);
      outputTree_->Branch((label.str()+"_phi").c_str(),&jet_phi_[ij]);
      outputTree_->Branch((label.str()+"_csv").c_str(),&jet_csv_[ij]);
      outputTree_->Branch((label.str()+"_jetid").c_str(),&jet_jetid_[ij]);
      outputTree_->Branch((label.str()+"_loosepuid").c_str(),&jet_loosepuid_[ij]);
      outputTree_->Branch((label.str()+"_tightpuid").c_str(),&jet_tightpuid_[ij]);

      outputTree_->Branch((label.str()+"_neutralhadfrac").c_str(),&jet_neutralhadfrac_[ij]);
      outputTree_->Branch((label.str()+"_neutralemfrac").c_str(),&jet_neutralemfrac_[ij]);
      outputTree_->Branch((label.str()+"_chargedemfrac").c_str(),&jet_chargedemfrac_[ij]);
      outputTree_->Branch((label.str()+"_chargedhadfrac").c_str(),&jet_chargedhadfrac_[ij]);
      outputTree_->Branch((label.str()+"_chargedmult").c_str(),&jet_chargedmult_[ij]);
      outputTree_->Branch((label.str()+"_neutralmult").c_str(),&jet_neutralmult_[ij]);

      outputTree_->Branch((label.str()+"_genjet_mindR").c_str(),&jet_genjet_mindR_[ij]);
      outputTree_->Branch((label.str()+"_flavour").c_str(),&jet_flavour_[ij]);
      outputTree_->Branch((label.str()+"_genid").c_str(),&jet_genid_[ij]);
      outputTree_->Branch((label.str()+"_genpt").c_str(),&jet_genpt_[ij]);
      outputTree_->Branch((label.str()+"_geneta").c_str(),&jet_geneta_[ij]);
      outputTree_->Branch((label.str()+"_genphi").c_str(),&jet_genphi_[ij]);
    }

    outputTree_->Branch("n_jets_cjv_30",&n_jets_cjv_30_);
    outputTree_->Branch("n_jets_cjv_20EB_30EE",&n_jets_cjv_20EB_30EE_);
    outputTree_->Branch("n_jets_15",&n_jets_15_);
    outputTree_->Branch("n_jets_csv2medium",&n_jets_csv2medium_);
    outputTree_->Branch("n_jets_30",&n_jets_30_);
    outputTree_->Branch("cjvjetpt",&cjvjetpt_);

    outputTree_->Branch("dijet_M",&dijet_M_);
    outputTree_->Branch("dijet_deta",&dijet_deta_);
    outputTree_->Branch("dijet_sumeta",&dijet_sumeta_);
    outputTree_->Branch("dijet_dphi",&dijet_dphi_);
    outputTree_->Branch("forward_tag_eta",&forward_tag_eta_);
    outputTree_->Branch("central_tag_eta",&central_tag_eta_);

    outputTree_->Branch("sumet",&sumet_);

    outputTree_->Branch("l1met",&l1met_);
    outputTree_->Branch("l1mht",&l1mht_);
    outputTree_->Branch("genmet",&genmet_);
    outputTree_->Branch("genmetphi",&genmetphi_);
    outputTree_->Branch("met",&met_);
    outputTree_->Branch("calomet",&calomet_);
    outputTree_->Branch("met_x",&met_x_);
    outputTree_->Branch("met_y",&met_y_);
    outputTree_->Branch("met_significance",&met_significance_);
    outputTree_->Branch("metnomuons",&metnomuons_);
    outputTree_->Branch("metnomu_x",&metnomu_x_);
    outputTree_->Branch("metnomu_y",&metnomu_y_);
    outputTree_->Branch("metnomu_significance",&metnomu_significance_);
    outputTree_->Branch("metnoelectrons",&metnoelectrons_);
    outputTree_->Branch("metnoel_x",&metnoel_x_);
    outputTree_->Branch("metnoel_y",&metnoel_y_);
    outputTree_->Branch("metnoel_significance",&metnoel_significance_);

    outputTree_->Branch("ht",&ht_);
    outputTree_->Branch("ht30",&ht30_);
    outputTree_->Branch("mht",&mht_);
    outputTree_->Branch("sqrt_ht",&sqrt_ht_);

    outputTree_->Branch("unclustered_et",&unclustered_et_);
    outputTree_->Branch("jetunclet_mindphi",&jetunclet_mindphi_);
    outputTree_->Branch("metnomuunclet_dphi",&metnomuunclet_dphi_);
    outputTree_->Branch("metnoelunclet_dphi",&metnoelunclet_dphi_);

    outputTree_->Branch("jetmetnomu_mindphi",&jetmetnomu_mindphi_);
    outputTree_->Branch("alljetsmetnomu_mindphi",&alljetsmetnomu_mindphi_);
    outputTree_->Branch("alljetsnotaumetnomu_mindphi",&alljetsnotaumetnomu_mindphi_);
    outputTree_->Branch("dijetmetnomu_scalarSum_pt",&dijetmetnomu_scalarSum_pt_);
    outputTree_->Branch("dijetmetnomu_vectorialSum_pt",&dijetmetnomu_vectorialSum_pt_);
    outputTree_->Branch("dijetmetnomu_ptfraction",&dijetmetnomu_ptfraction_);
    outputTree_->Branch("jet1metnomu_scalarprod",&jet1metnomu_scalarprod_);
    outputTree_->Branch("jet2metnomu_scalarprod",&jet2metnomu_scalarprod_);

    outputTree_->Branch("jetmetnoel_mindphi",&jetmetnoel_mindphi_);
    outputTree_->Branch("alljetsmetnoel_mindphi",&alljetsmetnoel_mindphi_);
    outputTree_->Branch("alljetsnotaumetnoel_mindphi",&alljetsnotaumetnoel_mindphi_);
    outputTree_->Branch("dijetmetnoel_scalarSum_pt",&dijetmetnoel_scalarSum_pt_);
    outputTree_->Branch("dijetmetnoel_vectorialSum_pt",&dijetmetnoel_vectorialSum_pt_);
    outputTree_->Branch("dijetmetnoel_ptfraction",&dijetmetnoel_ptfraction_);
    outputTree_->Branch("jet1metnoel_scalarprod",&jet1metnoel_scalarprod_);
    outputTree_->Branch("jet2metnoel_scalarprod",&jet2metnoel_scalarprod_);

    outputTree_->Branch("pass_muontrigger",&pass_muontrigger_);
    outputTree_->Branch("pass_photontrigger",&pass_photontrigger_);
    outputTree_->Branch("pass_sigtrigger",&pass_sigtrigger_);
    outputTree_->Branch("pass_mettrigger",&pass_mettrigger_);
    outputTree_->Branch("pass_metmht90trigger",&pass_metmht90trigger_);
    outputTree_->Branch("pass_metmht100trigger",&pass_metmht100trigger_);
    outputTree_->Branch("pass_metmht110trigger",&pass_metmht110trigger_);
    outputTree_->Branch("pass_metmht120trigger",&pass_metmht120trigger_);
    outputTree_->Branch("pass_controltrigger",&pass_controltrigger_);
    outputTree_->Branch("pass_singleEltrigger",&pass_singleEltrigger_);

    outputTree_->Branch("nvetomuons",&nvetomuons_);
    outputTree_->Branch("nselmuons",&nselmuons_);
    outputTree_->Branch("nvetoelectrons",&nvetoelectrons_);
    outputTree_->Branch("nselelectrons",&nselelectrons_);
    outputTree_->Branch("nvetotaus",&nvetotaus_);
    outputTree_->Branch("ntaus",&ntaus_);
    outputTree_->Branch("nloosephotons",&nloosephotons_);
    outputTree_->Branch("nmediumphotons",&nmediumphotons_);
    outputTree_->Branch("ntightphotons",&ntightphotons_);

    outputTree_->Branch("m_mumu",&m_mumu_);
    outputTree_->Branch("pt_mumu",&pt_mumu_);
    outputTree_->Branch("oppsign_mumu",&oppsign_mumu_);
    outputTree_->Branch("m_ee",&m_ee_);
    outputTree_->Branch("pt_ee",&pt_ee_);
    outputTree_->Branch("oppsign_ee",&oppsign_ee_);
    outputTree_->Branch("m_mumu_gen",&m_mumu_gen_);
    outputTree_->Branch("m_ee_gen",&m_ee_gen_);
    outputTree_->Branch("lep_mt",&lep_mt_);

    outputTree_->Branch("mu1_pt",&mu1_pt_);
    outputTree_->Branch("mu1_eta",&mu1_eta_);
    outputTree_->Branch("mu1_phi",&mu1_phi_);
    outputTree_->Branch("mu1_isTight",&mu1_isTight_);
    outputTree_->Branch("mu1_genmindR",&mu1_genmindR_);
    outputTree_->Branch("mu1_genpt",&mu1_genpt_);
    outputTree_->Branch("mu1_geneta",&mu1_geneta_);
    outputTree_->Branch("mu1_genphi",&mu1_genphi_);

    outputTree_->Branch("mu2_pt",&mu2_pt_);
    outputTree_->Branch("mu2_eta",&mu2_eta_);
    outputTree_->Branch("mu2_phi",&mu2_phi_);
    outputTree_->Branch("mu2_isTight",&mu2_isTight_);
    outputTree_->Branch("mu2_genmindR",&mu2_genmindR_);
    outputTree_->Branch("mu2_genpt",&mu2_genpt_);
    outputTree_->Branch("mu2_geneta",&mu2_geneta_);
    outputTree_->Branch("mu2_genphi",&mu2_genphi_);

    outputTree_->Branch("ele1_pt",&ele1_pt_);
    outputTree_->Branch("ele1_eta",&ele1_eta_);
    outputTree_->Branch("ele1_phi",&ele1_phi_);
    outputTree_->Branch("ele1_isTight",&ele1_isTight_);
    outputTree_->Branch("ele1_genmindR",&ele1_genmindR_);
    outputTree_->Branch("ele1_genpt",&ele1_genpt_);
    outputTree_->Branch("ele1_geneta",&ele1_geneta_);
    outputTree_->Branch("ele1_genphi",&ele1_genphi_);

    outputTree_->Branch("ele2_pt",&ele2_pt_);
    outputTree_->Branch("ele2_eta",&ele2_eta_);
    outputTree_->Branch("ele2_phi",&ele2_phi_);
    outputTree_->Branch("ele2_isTight",&ele2_isTight_);
    outputTree_->Branch("ele2_genmindR",&ele2_genmindR_);
    outputTree_->Branch("ele2_genpt",&ele2_genpt_);
    outputTree_->Branch("ele2_geneta",&ele2_geneta_);
    outputTree_->Branch("ele2_genphi",&ele2_genphi_);

    outputTree_->Branch("tau1_pt",&tau1_pt_);
    outputTree_->Branch("tau1_eta",&tau1_eta_);
    outputTree_->Branch("tau1_phi",&tau1_phi_);
    outputTree_->Branch("tau1_isTight",&tau1_isTight_);
    outputTree_->Branch("tau1_genmindR",&tau1_genmindR_);
    outputTree_->Branch("tau1_genpt",&tau1_genpt_);
    outputTree_->Branch("tau1_geneta",&tau1_geneta_);
    outputTree_->Branch("tau1_genphi",&tau1_genphi_);

    outputTree_->Branch("vetotau1_pt",&vetotau1_pt_);
    outputTree_->Branch("vetotau1_eta",&vetotau1_eta_);
    outputTree_->Branch("vetotau1_phi",&vetotau1_phi_);
    outputTree_->Branch("vetotau1_isTight",&vetotau1_isTight_);
    outputTree_->Branch("vetotau1_genmindR",&vetotau1_genmindR_);
    outputTree_->Branch("vetotau1_genpt",&vetotau1_genpt_);
    outputTree_->Branch("vetotau1_geneta",&vetotau1_geneta_);
    outputTree_->Branch("vetotau1_genphi",&vetotau1_genphi_);

    outputTree_->Branch("gamma1_pt",&gamma1_pt_);
    outputTree_->Branch("gamma1_eta",&gamma1_eta_);
    outputTree_->Branch("gamma1_phi",&gamma1_phi_);
    outputTree_->Branch("gamma1_genmindR",&gamma1_genmindR_);
    outputTree_->Branch("gamma1_genpt",&gamma1_genpt_);
    outputTree_->Branch("gamma1_geneta",&gamma1_geneta_);
    outputTree_->Branch("gamma1_genphi",&gamma1_genphi_);

    outputTree_->Branch("n_vertices",&n_vertices_);


//Branches (new) for Photons Electrons and Muons

    outputTree_->Branch("loosePhotons", &mynewloosephotons_);
    outputTree_->Branch("mediumPhotons",&mynewmediumphotons_);
    outputTree_->Branch("tightPhotons", &mynewtightphotons_);
    outputTree_->Branch("vetoElectrons",&mynewvetoelectrons_);
    outputTree_->Branch("selElectrons", &mynewselelectrons_);
    outputTree_->Branch("vetoMuons",    &mynewvetomuons_);
    outputTree_->Branch("selMuons",     &mynewselmuons_);

//Branches for gen_ Electrons and Muons

    outputTree_->Branch("gen_mu1_pt",       &gen_mu1_pt_);
    outputTree_->Branch("gen_mu1_eta",      &gen_mu1_eta_);
    outputTree_->Branch("gen_mu1_phi",      &gen_mu1_phi_);
    outputTree_->Branch("gen_mu1_mindR_j1", &gen_mu1_mindR_j1_);
    outputTree_->Branch("gen_mu1_mindR_j2", &gen_mu1_mindR_j2_);

    outputTree_->Branch("gen_ele1_pt",      &gen_ele1_pt_);
    outputTree_->Branch("gen_ele1_eta",     &gen_ele1_eta_);
    outputTree_->Branch("gen_ele1_phi",     &gen_ele1_phi_);
    outputTree_->Branch("gen_ele1_mindR_j1",&gen_ele1_mindR_j1_);
    outputTree_->Branch("gen_ele1_mindR_j2",&gen_ele1_mindR_j2_);


    // lheParticles
    outputTree_->Branch("lheHT",&lheHT_);

    // Boson pt
    outputTree_->Branch("boson_pt",&boson_pt_);

    return 0;
  }

  int  LightTreeRDM::Execute(TreeEvent *event){

    static unsigned processed = 0;

    resetAllTreeVariables();

    ///////////////////////////////
    // Get basic event variables //
    ///////////////////////////////

    EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
    run_= eventInfo->run();
    lumi_= eventInfo->lumi_block();
    event_= eventInfo->event();
    n_vertices_=eventInfo->good_vertices();

    if (debug_) {
      std::cout << " Run " << run_ << " event " << event_ << " n_vtx = " << n_vertices_ << std::endl;
    }


    ///////////////////////////
    // Get trigger variables //
    ///////////////////////////

    //if (is_data_){
    try {
      auto const& triggerPathPtrVec = event->GetPtrVec<TriggerPath>("triggerPathPtrVec","triggerPaths");

      for (unsigned i = 0; i < triggerPathPtrVec.size(); ++i) {
        std::string name = triggerPathPtrVec[i]->name();
        double prescale = triggerPathPtrVec[i]->prescale();
        //if (name.find("Photon") != name.npos){
        //std::cout << " Trigger " << name 
        //<< " prescale " << prescale
        //<< std::endl;
        //}
        if ( prescale==1 && ( ( name.find("HLT_Photon") != name.npos && (name.find("R9Id90_HE10_Iso") != name.npos || name.find("_PFMET") != name.npos) ) || name.find("HLT_Photon175_v") != name.npos ) ) pass_photontrigger_ = 1;
        //HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_v2
        //HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40
        //HLT_Photon36_R9Id90_HE10_IsoM
        ///HLT_Photon175_v*
        //HLT_Photon135_PFMET100

        //if (name.find("HLT_IsoMu20_") != name.npos) pass_muontrigger_ = prescale;
        if (name.find("HLT_IsoMu22_") != name.npos) pass_muontrigger_ = prescale;
        if (name.find("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140") != name.npos) pass_sigtrigger_ = prescale;
        if (name.find("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80") != name.npos) pass_controltrigger_ = prescale;
        if (name.find("HLT_PFMET170_") != name.npos) pass_mettrigger_ = prescale;
        if (name.find("HLT_PFMETNoMu90") != name.npos && name.find("PFMHTNoMu90") != name.npos) pass_metmht90trigger_ = prescale;
        if (name.find("HLT_PFMETNoMu100") != name.npos && name.find("PFMHTNoMu100") != name.npos) pass_metmht100trigger_ = prescale;
        if (name.find("HLT_PFMETNoMu110") != name.npos && name.find("PFMHTNoMu110") != name.npos) pass_metmht110trigger_ = prescale;
        if (name.find("HLT_PFMETNoMu120") != name.npos && name.find("PFMHTNoMu120") != name.npos) pass_metmht120trigger_ = prescale;
        if (name.find("HLT_Ele23_WPLoose_Gsf_v") != name.npos ||
            name.find("HLT_Ele27_WPLoose_Gsf_v") != name.npos ||
            name.find("HLT_Ele27_eta2p1_WPLoose_Gsf_v") != name.npos ||
            name.find("HLT_Ele27_WPTight_Gsf_v") != name.npos ||
            name.find("HLT_Ele35_WPLoose_Gsf_v") != name.npos ) pass_singleEltrigger_ = prescale;
      }
      if(do_trigskim_){
        if(!(pass_muontrigger_==1     ||
             pass_sigtrigger_>0       ||
             pass_controltrigger_>0   ||
             pass_mettrigger_>0       ||
             pass_metmht90trigger_>0  ||
             pass_metmht100trigger_>0 ||
             pass_metmht110trigger_>0 ||
             pass_metmht120trigger_>0 ||
             pass_photontrigger_>0    ||
             pass_singleEltrigger_>0    )){
          return 1;
        }
      }

      if (debug_) {
        std::cout << " -- Pass Sig         " << pass_sigtrigger_       << std::endl;
        std::cout << " -- Pass Muon        " << pass_muontrigger_      << std::endl;
        std::cout << " -- Pass Electron    " << pass_singleEltrigger_  << std::endl;
        std::cout << " -- Pass MET-MHT 90  " << pass_metmht90trigger_  << std::endl;
        std::cout << " -- Pass MET-MHT 100 " << pass_metmht100trigger_ << std::endl;
        std::cout << " -- Pass MET-MHT 110 " << pass_metmht110trigger_ << std::endl;
        std::cout << " -- Pass MET-MHT 120 " << pass_metmht120trigger_ << std::endl;
      }
    } catch (...) {
      //static bool trgbits = true;
      //if (trgbits) std::cout << " Trigger bits not filled ... " << std::endl;
      std::cout << " Trigger bits not filled ... " << std::endl;
      //trgbits = false; 
    }

    /////////////////
    // Get weights //
    /////////////////

    double wt = eventInfo->total_weight();
    //eventInfo->Print();

    double pileupwt=1;
    double pileupwtup=1;
    double pileupwtdown=1;

    if(!is_data_){
      v_nlo_Reweight_= eventInfo->weight_defined("v_nlo_Reweighting")?eventInfo->weight("v_nlo_Reweighting"):1;
      weight_lepveto_= eventInfo->weight("idisoVeto");
      weight_leptight_ = eventInfo->weight("idisoTight");
      pileupwt=eventInfo->weight("pileup");
      pileupwtup=eventInfo->weight("pileup_up");
      pileupwtdown=eventInfo->weight("pileup_down");

      //std::string label[7] = {"","_v0Up","_v0Down","_v1Up","_v1Down","_v2Up","_v2Down"};
      std::string label[7] = {"","_Up","_Down"};
      for (unsigned iT(0); iT<3; ++iT){
        std::string thislabel = "trig_2dbinned1d"+label[iT];
        weight_trig_[iT]=eventInfo->weight_defined(thislabel.c_str())?eventInfo->weight(thislabel.c_str()):0;
      }
    }
    if (pileupwt!=0) {
      puweight_up_scale_=pileupwtup/pileupwt;
      puweight_down_scale_=pileupwtdown/pileupwt;
    }
    else {
      puweight_up_scale_=1;
      puweight_down_scale_=1;
    }
    // Just for info: weight_trig_[0] = eventInfo->weight(trig_2dbinned1d)
    weight_nolepnotrig_ = wt;
    total_weight_lepveto_ =wt*weight_lepveto_*weight_trig_[0];
    total_weight_leptight_=wt*weight_leptight_*weight_trig_[0];

    if(total_weight_lepveto_!=total_weight_lepveto_){
      std::cout << " -- NAN lepveto weight: " << total_weight_lepveto_ << " " << wt << " " << weight_lepveto_ << " " << weight_trig_[0] << std::endl;//!!
    }

    //std::cout << " -- eventInfo->total_weight() : " << eventInfo->total_weight() << std::endl;
    //std::cout << " -- weight_nolepnotrig_ = eventInfo->weight(pileup)*eventInfo->weight(lumixs)*eventInfo->weight(wt_mc_sign) : " << eventInfo->weight("pileup")*eventInfo->weight("lumixs")*eventInfo->weight("wt_mc_sign") << std::endl;
    //std::cout << " -- total_weight_leptight_ = wt*weight_leptight_ : " << total_weight_leptight_ << std::endl;
    //std::cout << " -- wt (=eventInfo->total_weight())  : " << wt << std::endl;
    //std::cout << " -- weight_leptight_ : " << weight_leptight_ << std::endl;
    //std::cout << " -- trig_2dbinned1d : " << eventInfo->weight("trig_2dbinned1d") << std::endl;
    //std::cout << " -- What should be total_weight decomposed = eventInfo->weight(pileup)*eventInfo->weight(lumixs)*eventInfo->weight(wt_mc_sign)*eventInfo->weight(trig_2dbinned1d) : " << eventInfo->weight("pileup")*eventInfo->weight("lumixs")*eventInfo->weight("wt_mc_sign")*eventInfo->weight("trig_2dbinned1d") << std::endl;
    //std::cout << " -- What should be weight_trig_0 : " << weight_trig_[0] << std::endl;

    if (debug_) std::cout << " Event weight = " << wt << " pu up " << puweight_up_scale_ << " " << puweight_down_scale_ << std::endl;

    /////////////////////////
    // Get MET collections //
    /////////////////////////

    Met *met = 0;
    try {
      std::vector<Met*> & metCol = event->GetPtrVec<Met>(met_label_);
      met = metCol[0];
    } catch (...){
      //std::cout << " Met vec not found..." << std::endl;
      met = event->GetPtr<Met>(met_label_);
      if (!met) {
        std::cerr << " -- Found no MET " << met_label_ << " in event! Exiting..." << std::endl;
        exit(1);
      }
    }
    Met const* metnomuons = event->GetPtr<Met>("metNoMuons");
    Met const* metnoelectrons = event->GetPtr<Met>("metNoElectrons");

    ROOT::Math::PtEtaPhiEVector metvec = met->vector();
    ROOT::Math::PtEtaPhiEVector metnomuvec = metnomuons->vector();
    ROOT::Math::PtEtaPhiEVector metnoelvec = metnoelectrons->vector();


    met_ = met->pt();
    const Met::BasicMet & caloMet = met->GetCorrectedMet("RawCalo");
    calomet_ = caloMet.pt();
    //std::cout << " Check: met = " << met_ << " calomet = " << calomet_ << std::endl;

    met_x_ = metvec.Px();
    met_y_ = metvec.Py();
    met_significance_ = met->et_sig();
    sumet_ = met->sum_et();
    //if (debug_) std::cout << "l1metsize = " << l1met.size() << std::endl;
    //data have 5 bx saved, 2 before, 2 after. 
    //MC has only 1 but 3 objects (not sure why...)
    //if(!is_data_)
    //l1met_ = l1met[0]->pt();
    //else l1met_ = l1met[2]->pt();
    //}
    metnomuons_ = metnomuons->pt();
    metnomu_x_ = metnomuvec.Px();
    metnomu_y_ = metnomuvec.Py();
    if (met_>0) metnomu_significance_ = met_significance_/met_*metnomuons_;
    else metnomu_significance_ = met_significance_;


    metnoelectrons_ = metnoelectrons->pt();
    metnoel_x_ = metnoelvec.Px();
    metnoel_y_ = metnoelvec.Py();
    if (met_>0) metnoel_significance_ = met_significance_/met_*metnoelectrons_;
    else metnoel_significance_ = met_significance_;

    ////////////////
    // L1 objects //
    ////////////////
    try{
      std::vector<L1TObject*> const& l1sum = event->GetPtrVec<L1TObject>("l1tEtSum");
      //if(debug_ && l1met.size()!=1)std::cout<<"There seem to be "<<l1met.size()<<" l1mets!!"<<std::endl;
      //std::cout << " -- new event: l1sum size = " << l1sum.size() << " PFMET = " << met_ << std::endl;
      unsigned counterMet = 0;
      unsigned counterMht = 0;
      for (unsigned iL1(0); iL1<l1sum.size(); ++iL1){
        L1TObject* l1obj = l1sum[iL1];
        //std::cout << " -- sum " << iL1 << " bx " << l1obj->bx() << " type " << l1obj->sumType() << " pt " << l1obj->pt() << " " << l1obj->phi() << std::endl;
        if (l1obj->bx()!=0) continue;
        if (l1obj->sumType()==2) {
          //if (counterMet==2) 
          l1met_ = l1obj->pt();
          counterMet++;
        }
        if (l1obj->sumType()==3) {
          //if (counterMht==2) 
          l1mht_ = l1obj->pt();
          counterMht++;
        }
        //if (counterMet>2 && counterMht>2) break;
      }
      if (counterMet!=1 || counterMht != 1){
        std::cout << " *** Warning, found " << counterMet << " L1MET and " << counterMht << " L1MHT." << std::endl;
      }

      if (debug_) std::cout << " Met = " << met_ << " metnomu = " << metnomuons_ << " metnoel = " << metnoelectrons_ << " l1met = " << l1met_ << std::endl;

    } catch (...) {
      static bool firsttime = true;
      if (firsttime) std::cout << "No L1T trigger info in file, skipping.... " << std::endl;
      firsttime = false; 
    }

    //////////////////////////////
    // Get GenLevel collections //
    //////////////////////////////

    std::vector<Candidate*> genElecs;
    std::vector<Candidate*> genMus;
    std::vector<Candidate*> genTaus;
    std::vector<Candidate*> genPhotons;
    std::vector<GenJet *> genvec;

    if(!is_data_){

      // lheParticles and HT
      try{
        std::vector<GenParticle *> const& lheParticles = event->GetPtrVec<GenParticle>("lheParticles");
        lheHT_ = HTFromLHEParticles(lheParticles);
      } catch (...) {
        //std::cout << " -- lheParticles branch not found, go ahead ... " << std::endl;
      }

      genvec= event->GetPtrVec<GenJet>("genJets");
      std::sort(genvec.begin(), genvec.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
      nGenJets_ = genvec.size();

      std::vector<Met*> const& GenMet = event->GetPtrVec<Met>("genMetCollection");
      genmet_ = GenMet[0]->pt();
      genmetphi_ = GenMet[0]->phi();

      std::vector<GenParticle*> const& parts = event->GetPtrVec<GenParticle>("genParticles");
      std::vector<GenParticle*> Zmuons;
      std::vector<GenParticle*> Zelecs;

      //std::cout << " -- Event: " << event_ << " print of hardProcess particles..." << std::endl;
      for (unsigned iGenPart = 0; iGenPart < parts.size(); ++iGenPart) {//Loop over gen particles
        int id = parts[iGenPart]->pdgid();
        std::vector<bool> flags=parts[iGenPart]->statusFlags();
        if (flags[GenStatusBits::IsHardProcess] && flags[GenStatusBits::FromHardProcess] && flags[GenStatusBits::IsFirstCopy]){

          if (abs(id)==24 || abs(id)==23){// W+- || Z
            boson_pt_ = parts[iGenPart]->pt();
          }
          if (abs(id)==15){
            genTaus.push_back(parts[iGenPart]);
          }
          else if (parts[iGenPart]->mothers().size()==0 || (parts[iGenPart]->mothers().size()==1 && parts[parts[iGenPart]->mothers()[0]]->pdgid()==23)){
            if (abs(id)==11) Zelecs.push_back(parts[iGenPart]);
            else if (abs(id)==13) Zmuons.push_back(parts[iGenPart]);
          }
        }
        if (parts[iGenPart]->status()!=1) continue;
        if (abs(id)==11) genElecs.push_back(parts[iGenPart]);
        else if (abs(id)==13) genMus.push_back(parts[iGenPart]);
        else if (abs(id)==22) genPhotons.push_back(parts[iGenPart]);
      }//endof Loop over gen particles

      if (Zmuons.size()==2) {
        m_mumu_gen_ = (Zmuons[0]->vector()+Zmuons[1]->vector()).M();
        countZmumu_++;
      }
      else if (Zelecs.size()==2) {
        m_ee_gen_ = (Zelecs[0]->vector()+Zelecs[1]->vector()).M();
        countZee_++;
      }


      std::sort(genElecs.begin(), genElecs.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
      std::sort(genMus.begin(), genMus.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
      std::sort(genTaus.begin(), genTaus.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
      std::sort(genPhotons.begin(), genPhotons.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));


      //FILLING HERE THE NEW VARIABLES
      if(genElecs.size() != 0){
        gen_ele1_pt_       = genElecs[0]->pt();
        gen_ele1_eta_      = genElecs[0]->eta();
        gen_ele1_phi_      = genElecs[0]->phi();
      }

      if(genMus.size()   != 0){
        gen_mu1_pt_        = genMus[0]->pt();
        gen_mu1_eta_       = genMus[0]->eta();
        gen_mu1_phi_       = genMus[0]->phi();
      }


      if (debug_) std::cout << " nElecs = " << genElecs.size() << " nMu = " << genMus.size() << std::endl;
    }//for MC



    //////////////////////////
    // Get Jets collections //
    //////////////////////////
    //jet pair
    ROOT::Math::PtEtaPhiEVector jet1vec;
    ROOT::Math::PtEtaPhiEVector jet2vec;
    std::vector<CompositeCandidate *> const& dijet_vec = event->GetPtrVec<CompositeCandidate>(dijet_label_);

    if (dijet_vec.size() != 0) {

      CompositeCandidate const* dijet = dijet_vec.at(0);

      Candidate const* jet1 = dijet->GetCandidate("jet1");
      Candidate const* jet2 = dijet->GetCandidate("jet2");
      jet1vec = jet1->vector();
      jet2vec = jet2->vector();

      if(fabs(jet1->eta())>fabs(jet2->eta())){
        forward_tag_eta_=jet1->eta();
        central_tag_eta_=jet2->eta();
      }
      else{
        forward_tag_eta_=jet2->eta();
        central_tag_eta_=jet1->eta();
      }


      //FILLING HERE THE NEW VARIABLES
      if(genElecs.size() != 0){
        gen_ele1_mindR_j1_ = ROOT::Math::VectorUtil::DeltaR(genElecs[0]->vector(),jet1vec);
        gen_ele1_mindR_j2_ = ROOT::Math::VectorUtil::DeltaR(genElecs[0]->vector(),jet2vec);
      }

      if(genMus.size()   != 0){
        gen_mu1_mindR_j1_  = ROOT::Math::VectorUtil::DeltaR(genMus[0]->vector(),jet1vec);
        gen_mu1_mindR_j2_  = ROOT::Math::VectorUtil::DeltaR(genMus[0]->vector(),jet2vec);
      }


      dijet_M_ = dijet->M();
      dijet_deta_ = fabs(jet1->eta() - jet2->eta());
      dijet_sumeta_ = jet1->eta() + jet2->eta();
      dijet_dphi_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(jet1vec,jet2vec));

      double nomudphi1 = fabs(ROOT::Math::VectorUtil::DeltaPhi(jet1vec,metnomuvec));
      double nomudphi2 = fabs(ROOT::Math::VectorUtil::DeltaPhi(jet2vec,metnomuvec));
      jetmetnomu_mindphi_ = std::min(nomudphi1,nomudphi2);
      dijetmetnomu_scalarSum_pt_ = jet1->pt()+jet2->pt()+metnomuons->pt();
      dijetmetnomu_vectorialSum_pt_ = (jet1vec+jet2vec+metnomuvec).Pt();
      dijetmetnomu_ptfraction_ = dijet->pt()/(dijet->pt()+metnomuons->pt());
      jet1metnomu_scalarprod_ = (jet1vec.Px()*metnomu_x_+jet1vec.Py()*metnomu_y_)/met_;
      jet2metnomu_scalarprod_ = (jet2vec.Px()*metnomu_x_+jet2vec.Py()*metnomu_y_)/met_;

      double noeldphi1 = fabs(ROOT::Math::VectorUtil::DeltaPhi(jet1vec,metnoelvec));
      double noeldphi2 = fabs(ROOT::Math::VectorUtil::DeltaPhi(jet2vec,metnoelvec));
      jetmetnoel_mindphi_ = std::min(noeldphi1,noeldphi2);
      dijetmetnoel_scalarSum_pt_ = jet1->pt()+jet2->pt()+metnoelectrons->pt();
      dijetmetnoel_vectorialSum_pt_ = (jet1vec+jet2vec+metnoelvec).Pt();
      dijetmetnoel_ptfraction_ = dijet->pt()/(dijet->pt()+metnoelectrons->pt());
      jet1metnoel_scalarprod_ = (jet1vec.Px()*metnoel_x_+jet1vec.Py()*metnoel_y_)/met_;
      jet2metnoel_scalarprod_ = (jet2vec.Px()*metnoel_x_+jet2vec.Py()*metnoel_y_)/met_;


      if (debug_>1){
        std::cout << " Jet1 = ";
        jet1->Print();
        std::cout << " Jet2 = ";
        jet2->Print();
      }
    }//jet pair exists

    std::vector<PFJet*> & jets = event->GetPtrVec<PFJet>(jet_label_);
    std::sort(jets.begin(), jets.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));

    std::vector<std::pair<unsigned,bool> > recotogenmatch;
    recotogenmatch.resize(jets.size(),std::pair<unsigned,bool>(1000,false));
    if (!is_data_) getGenRecoMatches<PFJet,GenJet>(jets,genvec,recotogenmatch);

    ROOT::Math::PtEtaPhiEVector mhtVec(0,0,0,0);
    alljetsmetnomu_mindphi_=jetmetnomu_mindphi_;
    alljetsnotaumetnomu_mindphi_=1000;

    alljetsmetnoel_mindphi_=jetmetnoel_mindphi_;
    alljetsnotaumetnoel_mindphi_=1000;

    std::vector<Tau*> & taus=event->GetPtrVec<Tau>("taus");
    ntaus_=taus.size();
    std::vector<Tau*> & vetotaus=event->GetPtrVec<Tau>("vetoTaus");
    nvetotaus_=vetotaus.size();



    for (unsigned i = 0; i < jets.size(); ++i) {//Loop on jets
      ROOT::Math::PtEtaPhiEVector jetvec = jets[i]->vector();
      //just for first nJetsSave_ jets
      if (nJets_<nJetsSave_){
        jet_pt_[nJets_]          = jets[i]->pt();
        jet_uncorpt_[nJets_]     = jets[i]->uncorrected_energy()/jets[i]->energy()*jets[i]->pt();
        jet_neutralhadfrac_[nJets_] = jets[i]->neutral_had_energy()/jets[i]->uncorrected_energy();
        jet_chargedhadfrac_[nJets_] = jets[i]->charged_had_energy_frac();
	jet_neutralemfrac_[nJets_] = jets[i]->neutral_em_energy_frac();
	jet_chargedemfrac_[nJets_] = jets[i]->charged_em_energy_frac();
	jet_chargedmult_[nJets_] = jets[i]->charged_multiplicity();
	jet_neutralmult_[nJets_] = jets[i]->neutral_multiplicity();

        jet_E_[nJets_]   =jets[i]->energy();
        jet_eta_[nJets_] =jets[i]->eta();
        jet_phi_[nJets_] =jets[i]->phi();
	
	/*std::map<std::size_t, float> btaglist = jets[i]->b_discriminators();
	std::map<std::size_t, float>::iterator btagiter = btaglist.begin();
	for (unsigned ib=0; btagiter != btaglist.end();++ib,++btagiter){
	  std::cout << " tagger " << ib << " " << (btagiter->first) << " " << btagiter->second << std::endl;
	  }*/

        //jet_csv_[nJets_]=jets[i]->GetBDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        //jet_csv_[nJets_]=jets[i]->GetBDiscriminator("combinedSecondaryVertexBJetTags");
        //jet_jetid_[nJets_]=PFJetID2015(jets[i]);
        jet_csv_[nJets_]       = jets[i]->GetBDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet_jetid_[nJets_]     = PFJetID2016(jets[i]);
        jet_loosepuid_[nJets_] = PileupJetID(jets[i],4,false);
        jet_tightpuid_[nJets_] = PileupJetID(jets[i],4,true);
        jet_flavour_[nJets_]   = jets[i]->parton_flavour();

        if (!is_data_ && recotogenmatch[i].second){
          unsigned genid = recotogenmatch[i].first;
          jet_genjet_mindR_[nJets_]=ROOT::Math::VectorUtil::DeltaR(genvec[genid]->vector(),jetvec);
          jet_genid_[nJets_]=genid;
          jet_genpt_[nJets_]=genvec[genid]->pt();
          jet_geneta_[nJets_]=genvec[genid]->eta();
          jet_genphi_[nJets_]=genvec[genid]->phi();
        }
        nJets_++;
        if (debug_>1) {
          std::cout << " Jet " << i ;
          jets[i]->Print();
        }
      }
      //do stuff for all jets
      ht_+=jetvec.Et();
      if(jets[i]->pt()>30)	ht30_+=jetvec.Et();
      mhtVec += jetvec;

      if(jets[i]->pt()>15) n_jets_15_++;
      if(jets[i]->pt()>30) n_jets_30_++;
      if (fabs(jets[i]->eta())<2.4 && jets[i]->GetBDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.8) n_jets_csv2medium_++;
      //3rd jet
      if (i > 1) {
        double eta_high = -5;
        double eta_low = -5;

        if (jets[0]->eta() > jets[1]->eta()){
          eta_high = jets[0]->eta();
          eta_low = jets[1]->eta();
        }
        else {
          eta_high = jets[1]->eta();
          eta_low = jets[0]->eta();
        }

        bool isInCentralGap = fabs(jets[i]->eta())<4.7 && jets[i]->eta() > eta_low && jets[i]->eta() < eta_high;
        double tmppt=jets[i]->pt();
        if(isInCentralGap&&(tmppt>cjvjetpt_)){
          cjvjetpt_=tmppt;
        }
        if (jets[i]->pt() > 30.0 && isInCentralGap){
          ++n_jets_cjv_30_;
        }
        if ( ((jets[i]->eta()<2.4 && jets[i]->pt() > 20.0) ||
              (jets[i]->eta()>=2.4 && jets[i]->pt() > 30.0)) && 
            isInCentralGap){
          ++n_jets_cjv_20EB_30EE_;
        }
        if(jets[i]->pt()>30.0 && n_jets_30_<=5){
          double thisjetmetnomudphi = fabs(ROOT::Math::VectorUtil::DeltaPhi(jetvec,metnomuvec));
          if(thisjetmetnomudphi<alljetsmetnomu_mindphi_)alljetsmetnomu_mindphi_=thisjetmetnomudphi;
          //check matching with taus
          double thisjetmetnoeldphi = fabs(ROOT::Math::VectorUtil::DeltaPhi(jetvec,metnoelvec));
          if(thisjetmetnoeldphi<alljetsmetnoel_mindphi_)alljetsmetnoel_mindphi_=thisjetmetnoeldphi;
          bool noTauMatch=true;
          for (unsigned itau = 0; itau < ntaus_; ++itau) {//loop on taus
            ROOT::Math::PtEtaPhiEVector tauvec = taus[itau]->vector();
            if (ROOT::Math::VectorUtil::DeltaR(tauvec,jetvec)<0.4){
              noTauMatch=false;
              break;
            }
          }//loop on taus
          if(noTauMatch && thisjetmetnomudphi<alljetsnotaumetnomu_mindphi_)alljetsnotaumetnomu_mindphi_=thisjetmetnomudphi;
          if(noTauMatch && thisjetmetnoeldphi<alljetsnotaumetnoel_mindphi_)alljetsnotaumetnoel_mindphi_=thisjetmetnoeldphi;
        }
      }


      //if (nJets_>=nJetsSave_) {
        //std::cout << " -- Warning! maximum number of jets reached ! nJets = " << jets.size() << ". Saving only first " << nJetsSave_ << " jets." << std::endl;
        //break;
      //}
    }//endof Loop on jets

    ROOT::Math::PtEtaPhiEVector unclVec = mhtVec + metvec;

    mht_ = mhtVec.Et();
    sqrt_ht_ = sqrt(ht_);

    unclustered_et_ = unclVec.Et();
    metnomuunclet_dphi_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(unclVec,metnomuvec));
    metnoelunclet_dphi_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(unclVec,metnoelvec));

    //if jet1vec and jet2vec are properly 
    //set by the selected jet pair
    if (dijet_vec.size() != 0){
      jetunclet_mindphi_ = std::min(fabs(ROOT::Math::VectorUtil::DeltaPhi(jet1vec,unclVec)),
      fabs(ROOT::Math::VectorUtil::DeltaPhi(jet2vec,unclVec)));
    }

    if (debug_) std::cout << " mht = " << mht_ << std::endl;

   //////////////////////////////
    // Get Leptons collections //
    /////////////////////////////

    std::vector<Muon*> & vetomuons=event->GetPtrVec<Muon>("vetoMuons");
    std::vector<Muon*> & selmuons=event->GetPtrVec<Muon>("selMuons");
    std::vector<Electron*> & vetoelectrons=event->GetPtrVec<Electron>("vetoElectrons");
    std::vector<Electron*> & selelectrons=event->GetPtrVec<Electron>("selElectrons");

    for (unsigned vetoelectronselements = 0; vetoelectronselements<vetoelectrons.size(); ++vetoelectronselements){
    mynewvetoelectrons_.push_back(*vetoelectrons[vetoelectronselements]);
    }

    for (unsigned selelectronselements = 0; selelectronselements<selelectrons.size(); ++selelectronselements){
    mynewselelectrons_.push_back(*selelectrons[selelectronselements]);
    }

    for (unsigned vetomuonselements = 0; vetomuonselements<vetomuons.size(); ++vetomuonselements){
    mynewvetomuons_.push_back(*vetomuons[vetomuonselements]);
    }

    for (unsigned selmuonselements = 0; selmuonselements<selmuons.size(); ++selmuonselements){
    mynewselmuons_.push_back(*selmuons[selmuonselements]);
    }
    /*
    for (unsigned elements = 0; elements<vetoelectrons.size(); ++elements){
    mynewvetoelectrons_.push_back(*vetoelectrons[elements]);
    }
    for (unsigned elements = 0; elements<vetoelectrons.size(); ++elements){
    mynewvetoelectrons_.push_back(*vetoelectrons[elements]);
    }
    for (unsigned elements = 0; elements<vetoelectrons.size(); ++elements){
    mynewvetoelectrons_.push_back(*vetoelectrons[elements]);
    }*/
    if(!ignoreLeptons_){
      nvetomuons_=vetomuons.size();
      nvetoelectrons_=vetoelectrons.size();
    }
    else {
      nvetomuons_=0;
      nvetoelectrons_=0;
    }
    nselmuons_=selmuons.size();
    nselelectrons_=selelectrons.size();
    //fill lepton info using veto objects: veto=tight when requiring ==1 tight....
    //allows to use one loose one tight for mumu selection

    std::sort(vetomuons.begin(), vetomuons.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
    std::sort(vetoelectrons.begin(), vetoelectrons.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
    std::sort(taus.begin(), taus.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
    std::sort(vetotaus.begin(), vetotaus.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));

    std::vector<std::pair<unsigned,bool> > recotogen_elecs;
    recotogen_elecs.resize(vetoelectrons.size(),std::pair<unsigned,bool>(1000,false));
    if (!is_data_) getGenRecoMatches<Electron,Candidate>(vetoelectrons,genElecs,recotogen_elecs);

    std::vector<std::pair<unsigned,bool> > recotogen_muons;
    recotogen_muons.resize(vetomuons.size(),std::pair<unsigned,bool>(1000,false));
    if (!is_data_) getGenRecoMatches<Muon,Candidate>(vetomuons,genMus,recotogen_muons);

    std::vector<std::pair<unsigned,bool> > recotogen_taus;
    recotogen_taus.resize(taus.size(),std::pair<unsigned,bool>(1000,false));
    if (!is_data_) getGenRecoMatches<Tau,Candidate>(taus,genTaus,recotogen_taus);

    std::vector<std::pair<unsigned,bool> > recotogen_vetotaus;
    recotogen_vetotaus.resize(vetotaus.size(),std::pair<unsigned,bool>(1000,false));
    if (!is_data_) getGenRecoMatches<Tau,Candidate>(vetotaus,genTaus,recotogen_vetotaus);

    double lep_pt=0;
    double lep_phi=0;
    if(nvetomuons_>=1){
      mu1_pt_=vetomuons[0]->pt();
      mu1_eta_=vetomuons[0]->eta();
      mu1_phi_=vetomuons[0]->phi();
      mu1_isTight_ = isTightMuon(vetomuons[0],selmuons);
      if (nvetoelectrons_==0&&ntaus_==0){
        lep_pt=mu1_pt_;
        lep_phi=mu1_phi_;
      }
      if (!is_data_ && recotogen_muons[0].second){
        unsigned genid = recotogen_muons[0].first;
        mu1_genmindR_ = ROOT::Math::VectorUtil::DeltaR(genMus[genid]->vector(),vetomuons[0]->vector());
        mu1_genpt_=genMus[genid]->pt();
        mu1_geneta_=genMus[genid]->eta();
        mu1_genphi_=genMus[genid]->phi();
      }
      if(nvetomuons_>=2){
        mu2_pt_=vetomuons[1]->pt();
        mu2_eta_=vetomuons[1]->eta();
        mu2_phi_=vetomuons[1]->phi();
        mu2_isTight_ = isTightMuon(vetomuons[1],selmuons);
        m_mumu_=((vetomuons.at(0)->vector())+(vetomuons.at(1)->vector())).M();
        pt_mumu_=((vetomuons.at(0)->vector())+(vetomuons.at(1)->vector())).Pt();
        oppsign_mumu_ = vetomuons.at(0)->charge() != vetomuons.at(1)->charge();
        if (!is_data_ && recotogen_muons[1].second){
          unsigned genid = recotogen_muons[1].first;
          mu2_genmindR_ = ROOT::Math::VectorUtil::DeltaR(genMus[genid]->vector(),vetomuons[1]->vector());
          mu2_genpt_=genMus[genid]->pt();
          mu2_geneta_=genMus[genid]->eta();
          mu2_genphi_=genMus[genid]->phi();
        }
      }
    }
    if(nvetoelectrons_>=1){
      ele1_pt_=vetoelectrons[0]->pt();
      ele1_eta_=vetoelectrons[0]->eta();
      ele1_phi_=vetoelectrons[0]->phi();
      ele1_isTight_ = isTightElectron(vetoelectrons[0],selelectrons);
      if (nvetomuons_==0&&ntaus_==0){
        lep_pt=ele1_pt_;
        lep_phi=ele1_phi_;
      }
      if (!is_data_ && recotogen_elecs[0].second){
        unsigned genid = recotogen_elecs[0].first;
        ele1_genmindR_ = ROOT::Math::VectorUtil::DeltaR(genElecs[genid]->vector(),vetoelectrons[0]->vector());
        ele1_genpt_=genElecs[genid]->pt();
        ele1_geneta_=genElecs[genid]->eta();
        ele1_genphi_=genElecs[genid]->phi();
      }
      if(nvetoelectrons_>=2){
	ele2_pt_=vetoelectrons[1]->pt();
	ele2_eta_=vetoelectrons[1]->eta();
	ele2_phi_=vetoelectrons[1]->phi();
	ele2_isTight_ = isTightElectron(vetoelectrons[1],selelectrons);
	if (!is_data_ && recotogen_elecs[1].second){
	  unsigned genid = recotogen_elecs[1].first;
	  ele2_genmindR_ = ROOT::Math::VectorUtil::DeltaR(genElecs[genid]->vector(),vetoelectrons[1]->vector());
	  ele2_genpt_=genElecs[genid]->pt();
	  ele2_geneta_=genElecs[genid]->eta();
	  ele2_genphi_=genElecs[genid]->phi();
	}
        m_ee_=((vetoelectrons.at(0)->vector())+(vetoelectrons.at(1)->vector())).M();
        pt_ee_=((vetoelectrons.at(0)->vector())+(vetoelectrons.at(1)->vector())).Pt();
	oppsign_ee_ = vetoelectrons.at(0)->charge() != vetoelectrons.at(1)->charge();
      }
    }
    if(ntaus_>=1){
      tau1_pt_=taus[0]->pt();
      tau1_eta_=taus[0]->eta();
      tau1_phi_=taus[0]->phi();
      tau1_isTight_ = true;
      if (nvetoelectrons_==0&&nvetomuons_==0){
        lep_pt=tau1_pt_;
        lep_phi=tau1_phi_;
      }
      if (!is_data_ && recotogen_taus[0].second){
        unsigned genid = recotogen_taus[0].first;
        tau1_genmindR_ = ROOT::Math::VectorUtil::DeltaR(genTaus[genid]->vector(),taus[0]->vector());
        tau1_genpt_=genTaus[genid]->pt();
        tau1_geneta_=genTaus[genid]->eta();
        tau1_genphi_=genTaus[genid]->phi();
      }
    }
    if(nvetotaus_>=1){
      vetotau1_pt_=vetotaus[0]->pt();
      vetotau1_eta_=vetotaus[0]->eta();
      vetotau1_phi_=vetotaus[0]->phi();
      vetotau1_isTight_ = isTightTau(vetotaus[0],taus);
      if (!is_data_ && recotogen_vetotaus[0].second){
        unsigned genid = recotogen_vetotaus[0].first;
        vetotau1_genmindR_ = ROOT::Math::VectorUtil::DeltaR(genTaus[genid]->vector(),vetotaus[0]->vector());
        vetotau1_genpt_=genTaus[genid]->pt();
        vetotau1_geneta_=genTaus[genid]->eta();
        vetotau1_genphi_=genTaus[genid]->phi();
      }
    }


    lep_mt_ =sqrt(2*lep_pt*met->pt()*(1-cos(lep_phi-met->phi())));

    if (debug_) std::cout << " nelecs = " << nvetoelectrons_ << " nvetomuons = " << nvetomuons_ << std::endl;

    ////////////////////////////
    // Get Photon collections //
    ////////////////////////////
    try {
      std::vector<Photon*>  const& loosephotons=event->GetPtrVec<Photon>("loosePhotons");
      std::vector<Photon*>  const& mediumphotons=event->GetPtrVec<Photon>("mediumPhotons");
      std::vector<Photon*>  & tightphotons=event->GetPtrVec<Photon>("tightPhotons");
      nloosephotons_=loosephotons.size();
      nmediumphotons_=mediumphotons.size();
      ntightphotons_=tightphotons.size();
      std::sort(tightphotons.begin(), tightphotons.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
      
      std::vector<std::pair<unsigned,bool> > recotogen_photons;
      recotogen_photons.resize(tightphotons.size(),std::pair<unsigned,bool>(1000,false));
      if (!is_data_) getGenRecoMatches<Photon,Candidate>(tightphotons,genPhotons,recotogen_photons);
      
      for (unsigned loosephotonselements = 0; loosephotonselements<loosephotons.size(); ++loosephotonselements){
	mynewloosephotons_.push_back(*loosephotons[loosephotonselements]);
      }
      
      for (unsigned mediumphotonselements = 0; mediumphotonselements<mediumphotons.size(); ++mediumphotonselements){
	mynewmediumphotons_.push_back(*mediumphotons[mediumphotonselements]);
      }
      
      for (unsigned tightphotonselements = 0; tightphotonselements<tightphotons.size(); ++tightphotonselements){
	mynewtightphotons_.push_back(*tightphotons[tightphotonselements]);
      }
      
      if(ntightphotons_>=1){
	gamma1_pt_=tightphotons[0]->pt();
	gamma1_eta_=tightphotons[0]->eta();
	gamma1_phi_=tightphotons[0]->phi();
	if (!is_data_ && recotogen_photons[0].second){
	  unsigned genid = recotogen_photons[0].first;
	  gamma1_genmindR_ = ROOT::Math::VectorUtil::DeltaR(genPhotons[genid]->vector(),tightphotons[0]->vector());
	  gamma1_genpt_=genPhotons[genid]->pt();
	  gamma1_geneta_=genPhotons[genid]->eta();
	  gamma1_genphi_=genPhotons[genid]->phi();
	}
      }
      
      
      if (debug_) std::cout << " nPhotons = " << nloosephotons_ << " " << nmediumphotons_ << " " << ntightphotons_ << std::endl;
    } catch (...) {
      //if photon collections not there, ignore....
    }

    //IF PASSES CUTS FILL TREE    
    //if(do_noskim_){
    //but at least one interesting offline object....
    //if (n_jets_30_>1 || 
    //nselmuons_>0 || nselelectrons_>0 || ntaus_ > 0 ||
    //nloosephotons_ > 0 || nmediumphotons_ > 0 || ntightphotons_ > 0){
    outputTree_->Fill();
    //}
    //++processed;
    //}
    /*else if(!ignoreLeptons_){
      if(!do_promptskim_){
	if (jet_pt_[1]>40&& dijet_M_ > 600 &&  dijet_deta_>3.6){
	  //if (dijet_M_>1000 &&  dijet_deta_>3.6 && metnomuons_>100 && jet1_pt_>50){//for prompt presel
	  outputTree_->Fill();
	  //++processed;
	}
      }
      else{
	if ((pass_sigtrigger_>0||pass_controltrigger_>0)&&dijet_deta_>3.6&&metnomuons_>90&&jet_pt_[0]>50){
	  //if (dijet_M_>1000 &&  dijet_deta_>3.6 && metnomuons_>100 && jet1_pt_>50){//for prompt presel
	  outputTree_->Fill();
	  //++processed;
	}
      }
    }//Not ignoreLeptons
    else{
      if(!do_promptskim_){
	if (metnomu_significance_ > 3.0 &&  dijet_deta_>3.6 &&m_mumu_gen_>80&&m_mumu_gen_<100){
	  //if (dijet_M_>1000 &&  dijet_deta_>3.6 && metnomuons_>100 && jet1_pt_>50){//for prompt presel
	  outputTree_->Fill();
	  //++processed;
	}
      }
      else{
	if ((pass_sigtrigger_>0||pass_controltrigger_>0)&&dijet_deta_>3.6&&metnomuons_>90&&jet_pt_[0]>50){
	  //if (dijet_M_>1000 &&  dijet_deta_>3.6 && metnomuons_>100 && jet1_pt_>50){//for prompt presel
	  outputTree_->Fill();
	  //++processed;
	}
      }	
    }
    */
    
    ++processed;
    //if (processed == 500) outputTree_->OptimizeBaskets();
    if ((processed%100000) == 0) outputTree_->AutoSave();

    return 0;


  }//execute method

  int  LightTreeRDM::PostAnalysis(){

    fs_->cd();
    outputTree_->Write();

    std::cout << "----------------------------------------" << std::endl
              << "-- PostAnalysis Info for LightTreeRDM --" << std::endl
              << "----------------------------------------" << std::endl
              << " -- Number of Zmumu at ME level: " << countZmumu_ << std::endl
              << " -- Number of Zee at ME level: " << countZee_ << std::endl;

    return 0;
  }

  void  LightTreeRDM::PrintInfo(){
    ;
  }

}//namespace


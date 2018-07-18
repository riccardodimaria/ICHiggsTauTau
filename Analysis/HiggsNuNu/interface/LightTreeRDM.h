#ifndef ICHiggsTauTau_Module_LightTreeRDM_h
#define ICHiggsTauTau_Module_LightTreeRDM_h

#include <string>

#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "UserCode/ICHiggsTauTau/interface/Candidate.hh"
#include "UserCode/ICHiggsTauTau/interface/Photon.hh"
#include "UserCode/ICHiggsTauTau/interface/Electron.hh"
#include "UserCode/ICHiggsTauTau/interface/Muon.hh"
#include "UserCode/ICHiggsTauTau/interface/Tau.hh"
#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "UserCode/ICHiggsTauTau/interface/GenJet.hh"

#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/BTagCalibrationStandalone.h"


#include "TTree.h"
#include "TVector3.h"
#include "Math/VectorUtil.h"


namespace ic {
  class LightTreeRDM : public ModuleBase {

  private:

    CLASS_MEMBER(LightTreeRDM,fwlite::TFileService*, fs);
    CLASS_MEMBER(LightTreeRDM,unsigned, debug);
    CLASS_MEMBER(LightTreeRDM,std::string, met_label);
    CLASS_MEMBER(LightTreeRDM,std::string, pfmet_label);
    CLASS_MEMBER(LightTreeRDM,std::string, jet_label);
    CLASS_MEMBER(LightTreeRDM,std::string, dijet_label);
    CLASS_MEMBER(LightTreeRDM,bool, is_data);
    CLASS_MEMBER(LightTreeRDM,bool, ignoreLeptons);
    CLASS_MEMBER(LightTreeRDM,bool, do_jetht);
    CLASS_MEMBER(LightTreeRDM,bool, do_trigskim);
    CLASS_MEMBER(LightTreeRDM,bool, do_promptskim);
    CLASS_MEMBER(LightTreeRDM,bool, do_noskim);

    unsigned countZmumu_;
    unsigned countZee_;
    unsigned countQuarks_;
    unsigned countGenjets_;

    unsigned *Eht_;
    unsigned *Ejet_;

    TTree *outputTree_;

    unsigned run_;
    unsigned lumi_;
    unsigned event_;
    double weight_nolepnotrig_;
    double weight_trig_[2];
    double weight_lepveto_;
    double weight_leptight_;
    double total_weight_lepveto_;
    double total_weight_leptight_;
    double puweight_up_scale_;
    double puweight_down_scale_;
    double weight_pileup_;
    double weight_xsection_;

    double v_nlo_Reweight_;
    double ewk_v_nlo_Reweight_;
    double prefiring_Reweight_;
    double prefiring_SingleMuon_Reweight_;
    double prefiring_JetHT_Reweight_;

    double weight_eleVeto_[7];
    double weight_muVeto_[7];
    double weight_eleTight_[3];
    double weight_gsfTight_[3];
    double weight_muidTight_[3];
    double weight_muisoTight_[3];
    double weight_mutkTight_[3];

    double weight_0muloose_;
    double weight_0muloose_idup_;
    double weight_0muloose_iddown_;
    double weight_0muloose_isoup_;
    double weight_0muloose_isodown_;
    double weight_0muloose_tkup_;
    double weight_0muloose_tkdown_;

    double weight_eletrigEff_[3];

    double weight_0b_alljets_;
    double weight_0bup_alljets_;
    double weight_0bdown_alljets_;
    double weight_0b_extrajets_;
    double weight_0bup_extrajets_;
    double weight_0bdown_extrajets_;

    unsigned nJetsSave_;
    unsigned nJets_;
    unsigned nGenJets_;

    /*std::vector<double> jet_pt_;
    std::vector<double> jet_E_;
    std::vector<double> jet_eta_;
    std::vector<double> jet_phi_;
    std::vector<double> jet_csv_;
    std::vector<double> jet_jetid_;
    std::vector<double> jet_puid_;
    std::vector<int> jet_flavour_;

    std::vector<double> jet_genjet_mindR_;
    std::vector<unsigned> jet_genid_;
    std::vector<double> jet_genpt_;
    std::vector<double> jet_geneta_;
    std::vector<double> jet_genphi_;*/

// Define mynew* -- Electron.hh also included

    std::vector<Photon> mynewloosephotons_;
    std::vector<Photon> mynewmediumphotons_;
    std::vector<Photon> mynewtightphotons_;
    std::vector<Electron> mynewvetoelectrons_;
    std::vector<Electron> mynewselelectrons_;
    std::vector<Muon> mynewvetomuons_;
    std::vector<Muon> mynewselmuons_;


    double* jet_pt_;
    double* jet_uncorpt_;
    double* jet_E_;
    double* jet_eta_;
    double* jet_phi_;
    double* jet_csv_;
    double* jet_jetid_;
    double* jet_loosepuid_;
    double* jet_tightpuid_;
    int* jet_flavour_;
    double* jet_neutralhadfrac_;
    double* jet_neutralemfrac_;
    double* jet_chargedemfrac_;
    double* jet_chargedhadfrac_;
    unsigned* jet_chargedmult_;
    unsigned* jet_neutralmult_;


    double* jet_genjet_mindR_;
    unsigned* jet_genid_;
    double* jet_genpt_;
    double* jet_geneta_;
    double* jet_genphi_;

    unsigned n_jets_cjv_30_;
    unsigned n_jets_cjv_20EB_30EE_;
    unsigned n_jets_15_;
    unsigned n_jets_30_;
    unsigned n_jets_csv2medium_;
    unsigned n_extrajets_csv2medium_;
    double cjvjetpt_;

    double dijet_M_;
    double dijet_deta_;
    double dijet_sumeta_;
    double dijet_dphi_;
    double forward_tag_eta_;
    double central_tag_eta_;

    double sumet_;

    double l1met_;
    double l1mht_;
    double met_;
    double met_phi_;
    double calomet_;
    double calomet_phi_;
    double genmet_;
    double genmetphi_;
    double metnomuons_;
    double metnomuons_phi_;
    double metnoelectrons_;
    double metnoelectrons_phi_;
    double met_x_;
    double met_y_;
    double metnomu_x_;
    double metnomu_y_;
    double metnoel_x_;
    double metnoel_y_;
    double met_significance_;
    double metnomu_significance_;
    double metnoel_significance_;

    double ht_;
    double ht30_;
    double mht_;
    double sqrt_ht_;

    double ht_eta3d0_;
    double ht30_eta3d0_;
    double mht30_eta3d0_;
    double ht_eta2d4_;
    double ht30_eta2d4_;
    double mht30_eta2d4_;

    double unclustered_et_;
    double jetunclet_mindphi_;
    double metnomuunclet_dphi_;
    double metnoelunclet_dphi_;

    double threejetsmetnomu_mindphi_;
    double threejetsnotaumetnomu_mindphi_;
    double fourjetsmetnomu_mindphi_;
    double fourjetsnotaumetnomu_mindphi_;
    double fivejetsmetnomu_mindphi_;
    double fivejetsnotaumetnomu_mindphi_;
    double sixjetsmetnomu_mindphi_;
    double sixjetsnotaumetnomu_mindphi_;
    double alljetsmetnomu_mindphi_;
    double alljetsnotaumetnomu_mindphi_;
    double jetmetnomu_mindphi_;

    double threejetsmetnoel_mindphi_;
    double threejetsnotaumetnoel_mindphi_;
    double fourjetsmetnoel_mindphi_;
    double fourjetsnotaumetnoel_mindphi_;
    double fivejetsmetnoel_mindphi_;
    double fivejetsnotaumetnoel_mindphi_;
    double sixjetsmetnoel_mindphi_;
    double sixjetsnotaumetnoel_mindphi_;
    double alljetsmetnoel_mindphi_;
    double alljetsnotaumetnoel_mindphi_;
    double jetmetnoel_mindphi_;

    double dijetmetnomu_scalarSum_pt_;
    double dijetmetnomu_vectorialSum_pt_;
    double dijetmetnomu_ptfraction_;
    double jet1metnomu_scalarprod_;
    double jet2metnomu_scalarprod_;

    double dijetmetnoel_scalarSum_pt_;
    double dijetmetnoel_vectorialSum_pt_;
    double dijetmetnoel_ptfraction_;
    double jet1metnoel_scalarprod_;
    double jet2metnoel_scalarprod_;

    int pass_muontrigger_;
    int pass_photontrigger_;
    int pass_sigtrigger_;
    int pass_mettrigger_;
    int pass_metmht90trigger_;
    int pass_metmht100trigger_;
    int pass_metmht110trigger_;
    int pass_metmht120trigger_;
    int pass_controltrigger_;
    int pass_singleEltrigger_;

    int pass_httrigger_[11];
    int pass_jettrigger_[10];
    int pass_alljethttrigger_;

    unsigned nvetomuons_;
    unsigned nselmuons_;
    unsigned nvetoelectrons_;
    unsigned nselelectrons_;
    unsigned nvetotaus_;
    unsigned ntaus_;
    unsigned nloosephotons_;
    unsigned nmediumphotons_;
    unsigned ntightphotons_;

    double m_mumu_;
    double pt_mumu_;
    bool   oppsign_mumu_;
    double m_ee_;
    double pt_ee_;
    bool   oppsign_ee_;
    double m_mumu_gen_;
    double m_ee_gen_;
    double lep_mt_;

//  Define gen_mu1_*

    double gen_mu1_pt_;
    double gen_mu1_eta_;
    double gen_mu1_phi_;
    double gen_mu1_mindR_j1_;
    double gen_mu1_mindR_j2_;
    double gen_mu2_pt_;
    double gen_mu2_eta_;


    double mu1_pt_;
    double mu1_eta_;
    double mu1_phi_;
    bool   mu1_isTight_;
    double mu1_genmindR_;
    double mu1_genpt_;
    double mu1_geneta_;
    double mu1_genphi_;

    double mu2_pt_;
    double mu2_eta_;
    double mu2_phi_;
    bool   mu2_isTight_;
    double mu2_genmindR_;
    double mu2_genpt_;
    double mu2_geneta_;
    double mu2_genphi_;

//  Define gen_ele1_*

    double gen_ele1_pt_;
    double gen_ele1_eta_;
    double gen_ele1_phi_;
    double gen_ele1_mindR_j1_;
    double gen_ele1_mindR_j2_;
    double gen_ele2_pt_;
    double gen_ele2_eta_;


    double ele1_pt_;
    double ele1_eta_;
    double ele1_phi_;
    bool   ele1_isTight_;
    double ele1_genmindR_;
    double ele1_genpt_;
    double ele1_geneta_;
    double ele1_genphi_;
    double ele1_relISOwithEA_;

    double ele2_pt_;
    double ele2_eta_;
    double ele2_phi_;
    bool   ele2_isTight_;
    double ele2_genmindR_;
    double ele2_genpt_;
    double ele2_geneta_;
    double ele2_genphi_;
    double ele2_relISOwithEA_;


//  Define gen_tau1_*

    double gen_tau1_pt_;
    double gen_tau1_eta_;
    double gen_tau1_phi_;
    double gen_tau1_mindR_j1_;
    double gen_tau1_mindR_j2_;
    double gen_tau2_pt_;
    double gen_tau2_eta_;


    double tau1_pt_;
    double tau1_eta_;
    double tau1_phi_;
    bool   tau1_isTight_;
    double tau1_genmindR_;
    double tau1_genpt_;
    double tau1_geneta_;
    double tau1_genphi_;

    double vetotau1_pt_;
    double vetotau1_eta_;
    double vetotau1_phi_;
    bool   vetotau1_isTight_;
    double vetotau1_genmindR_;
    double vetotau1_genpt_;
    double vetotau1_geneta_;
    double vetotau1_genphi_;

    double gamma1_pt_;
    double gamma1_eta_;
    double gamma1_phi_;
    double gamma1_genmindR_;
    double gamma1_genpt_;
    double gamma1_geneta_;
    double gamma1_genphi_;

    unsigned n_vertices_;
    float n_true_int_;

    // lheParticles
    double lheHT_;

    // Boson pt
    double boson_pt_;

    //dijet mass from genjets from VBF quarks
    double vbf_digenjet_m_;
    double vbf_diquark_m_;

    const BTagCalibration *calib;
    BTagCalibrationReader* reader_comb;

  public:
    LightTreeRDM(std::string const& name);
    virtual ~LightTreeRDM();

    void resetAllTreeVariables();

    virtual int PreAnalysis();
    virtual int Execute(TreeEvent *event);
    virtual int PostAnalysis();
    virtual void PrintInfo();

   };


}


#endif

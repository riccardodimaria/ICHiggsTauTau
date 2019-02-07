#ifndef ICHiggsTauTau_Module_HinvWeights_h
#define ICHiggsTauTau_Module_HinvWeights_h

#include "UserCode/ICHiggsTauTau/interface/GenParticle.hh"
#include "UserCode/ICHiggsTauTau/interface/TH2DAsymErr.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/BTagWeight.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvConfig.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/HTFromLHEParticles.h"
#include <string>
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "TH3F.h"

namespace ic {

class HinvWeights : public ModuleBase {
 private:
  CLASS_MEMBER(HinvWeights, fwlite::TFileService*, fs);
  CLASS_MEMBER(HinvWeights, ic::mc, mc)
  CLASS_MEMBER(HinvWeights, ic::era, era)
  CLASS_MEMBER(HinvWeights, bool, save_weights)
  CLASS_MEMBER(HinvWeights, bool, save_lumixs_weights)  
  CLASS_MEMBER(HinvWeights, bool, do_top_reweighting)
  CLASS_MEMBER(HinvWeights, bool, do_trg_weights)

  CLASS_MEMBER(HinvWeights, bool, trg_applied_in_mc)
  CLASS_MEMBER(HinvWeights, bool, do_idiso_tight_weights)
  CLASS_MEMBER(HinvWeights, bool, do_idiso_veto_weights)
  CLASS_MEMBER(HinvWeights, bool, do_w_soup)
  CLASS_MEMBER(HinvWeights, bool, do_w_reweighting)
  CLASS_MEMBER(HinvWeights, bool, do_ewk_w_reweighting)
  CLASS_MEMBER(HinvWeights, bool, do_dy_soup)
  CLASS_MEMBER(HinvWeights, bool, do_dy_soup_htbinned)
  CLASS_MEMBER(HinvWeights, bool, do_dy_reweighting)
  CLASS_MEMBER(HinvWeights, bool, do_ewk_dy_reweighting)
  CLASS_MEMBER(HinvWeights, std::string, input_met)
  CLASS_MEMBER(HinvWeights, std::string, input_jet)
  CLASS_MEMBER(HinvWeights, bool, do_lumixs_weights)
  CLASS_MEMBER(HinvWeights, std::string, input_params)
  CLASS_MEMBER(HinvWeights, std::string, sample_name)
  CLASS_MEMBER(HinvWeights, std::string, mettrg_weight_file)
  CLASS_MEMBER(HinvWeights, std::string, mettrg_zmm_weight_file)

  // For v_nlo_Reweighting (kfactors.root file in input/scalefactors from MIT group)
  CLASS_MEMBER(HinvWeights, std::string, kfactors_wjets_file)
  CLASS_MEMBER(HinvWeights, std::string, kfactors_zjets_file)
  TFile *kfactors_wjets_;
  TH1F *hist_kfactors_qcdewk_W;
  TH1F *hist_kfactors_vbf_cnc_W;
  TFile *kfactors_zjets_;
  TH1F *hist_kfactors_qcdewk_Z;
  TH1F *hist_kfactors_vbf_cnc_Z;

  CLASS_MEMBER(HinvWeights, std::string, kFactor_ZToNuNu_pT_Mjj_file)
  CLASS_MEMBER(HinvWeights, std::string, kFactor_WToLNu_pT_Mjj_file)
  TFile *kFactor_ZToNuNu_pT_Mjj_;
  TFile *kFactor_WToLNu_pT_Mjj_;
  TH2F *hist_kFactors_ewk_Z;
  TH2F *hist_kFactors_ewk_W;

  TFile *mettrigSF_;
  TFile *mettrigZmmSF_;
  std::vector<TF1*> func_mettrigSF_;
  std::vector<TF1*> func_mettrigZmmSF_;
  std::vector<unsigned> trigMjjBins_;

  TH1F *tighteleweight;
  TH1F *tightmuweight;
  TH1F *vetoeleweight;
  TH1F *vetomuweight;

  double lumixsweight;

  //[3]=central-up-down
  std::vector<double> eTight_idisoSF_[3];
  std::vector<double> e_trigDataEff_[3];
  std::vector<double> eVeto_idisoDataEff_[3];
  std::vector<double> eVeto_idisoMCEff_[3];
  std::vector<double> e_gsfidSF_[3];
  std::vector<double> e_gsfidDataEff_[3];
  std::vector<double> e_gsfidMCEff_[3];
  std::vector<double> muTight_idSF_[3];
  std::vector<double> muTight_isoSF_[3];
  std::vector<double> muVeto_idSF_[3];
  std::vector<double> muVeto_isoSF_[3];
  std::vector<double> muVeto_idDataEff_[3];
  std::vector<double> muVeto_isoDataEff_[3];
  std::vector<double> muVeto_idMCEff_[3];
  std::vector<double> muVeto_isoMCEff_[3];
  std::vector<double> mu_tkSF_[3];

  std::vector<double> gsf_etabin_;
  std::vector<double> gsf_ptbin_;
  std::vector<double> e_etabin_;
  std::vector<double> e_ptbin_;
  std::vector<double> e_etatrig_;
  std::vector<double> e_pttrig_;

  std::vector<double> tk_etabin_;
  std::vector<double> tk_ptbin_;
  std::vector<double> mu_etabin_;
  std::vector<double> mu_ptbin_;

  double f0_,f1_,f2_,f3_,f4_,f5_,f6_,f7_,n_inc_,n1_,n2_,n3_,n4_,n5_,n6_,n7_,w0_,w1_,w2_,w3_,w4_,w5_,w6_,w7_;
  double zf0_,zf1_,zf2_,zf3_,zf4_,zn_inc_,zn1_,zn2_,zn3_,zn4_,zw0_,zw1_,zw2_,zw3_,zw4_;

  unsigned eventsWithGenElectron_;
  unsigned eventsWithGenMuon_;
  unsigned eventsWithGenTau_;
  unsigned eventsWithGenElectronFromTau_;
  unsigned eventsWithGenMuonFromTau_;
  unsigned eventsWithGenElectronInAcc_;
  unsigned eventsWithGenMuonInAcc_;
  unsigned eventsWithGenElectronFromTauInAcc_;
  unsigned eventsWithGenMuonFromTauInAcc_;

 public:
  HinvWeights(std::string const& name);
  virtual ~HinvWeights();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
  double Efficiency(double m, double m0, double sigma, double alpha, double n, double norm);
  void SetWTargetFractions(double f0, double f1, double f2, double f3, double f4);
  void SetWTargetFractions(double f0, double f1, double f2, double f3, double f4, double f5, double f6, double f7);
  void SetWInputYields(double n_inc, double n1, double n2, double n3, double n4);
  void SetWInputYields(double n_inc, double n1, double n2, double n3, double n4, double n5, double n6, double n7);
  void SetDYTargetFractions(double zf0, double zf1, double zf2, double zf3, double zf4);
  void SetDYInputYields(double zn_inc, double zn1, double zn2, double zn3, double zn4);

  unsigned getPartonNumber(std::vector<GenParticle*> const& parts);

  unsigned findPtEtaBin(const double & pt, const double & eta, const std::vector<double> & ptbin, const std::vector<double> & etabin);

  void fillVector(const std::string & aFileName, 
                  const unsigned nPtBins,
                  const unsigned nEtaBins,
                  std::vector<double> aVector[3],
                  std::vector<double> & ptbin,
                  std::vector<double> & etabin,
		  const bool protect);

  double nloReweighting(const double & aMjj, const double & aYstar);


};

}
#undef MEMBER_NP
#endif

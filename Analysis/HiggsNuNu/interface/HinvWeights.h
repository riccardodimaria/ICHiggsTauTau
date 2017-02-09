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
  CLASS_MEMBER(HinvWeights, bool, do_1dparkedtrg_weights)
  CLASS_MEMBER(HinvWeights, bool, do_fitted1dparkedtrg_weights)
  CLASS_MEMBER(HinvWeights, bool, do_3dtrg_weights)
  CLASS_MEMBER(HinvWeights, bool, do_binnedin2d1dfittedtrg_weights)
  CLASS_MEMBER(HinvWeights, std::vector<std::string>, binnedin2d1dfitweightvarorder)//bin in first two fit in 3rd
  CLASS_MEMBER(HinvWeights, std::vector<double>, binnedin2d1dfitweightvar1binning)//binning of first var
  CLASS_MEMBER(HinvWeights, std::vector<double>, binnedin2d1dfitweightvar2binning)//binning of second var
  CLASS_MEMBER(HinvWeights, bool, do_run2)
  CLASS_MEMBER(HinvWeights, bool, do_metmht)


  CLASS_MEMBER(HinvWeights, bool, trg_applied_in_mc)
  CLASS_MEMBER(HinvWeights, bool, do_idiso_tight_weights)
  CLASS_MEMBER(HinvWeights, bool, do_idiso_veto_weights)
  CLASS_MEMBER(HinvWeights, bool, do_w_soup)
  CLASS_MEMBER(HinvWeights, bool, do_w_reweighting)
  CLASS_MEMBER(HinvWeights, bool, do_dy_soup)
  CLASS_MEMBER(HinvWeights, bool, do_dy_soup_htbinned)
  CLASS_MEMBER(HinvWeights, bool, do_dy_reweighting)
  CLASS_MEMBER(HinvWeights, std::string, input_met)
  CLASS_MEMBER(HinvWeights, std::string, input_jet)
  CLASS_MEMBER(HinvWeights, bool, do_idiso_err)
  CLASS_MEMBER(HinvWeights, bool, do_idiso_errmuore)
  CLASS_MEMBER(HinvWeights, bool, do_idiso_errupordown)
  CLASS_MEMBER(HinvWeights, bool, do_lumixs_weights)
  CLASS_MEMBER(HinvWeights, std::string, input_params)
  CLASS_MEMBER(HinvWeights, std::string, sample_name)
  CLASS_MEMBER(HinvWeights, std::string, trg_weight_file)
  CLASS_MEMBER(HinvWeights, double, Alumi)
  CLASS_MEMBER(HinvWeights, double, BClumi)
  CLASS_MEMBER(HinvWeights, double, Dlumi)

  // For v_nlo_Reweighting (kfactors.root file in input/scalefactors from MIT group)
  CLASS_MEMBER(HinvWeights, std::string, kfactors_file)
  TFile *kfactors_;
  TH1F *hist_kfactors_EWKcorr_W;
  TH1F *hist_kfactors_WJets_012j_NLO;
  TH1F *hist_kfactors_EWKcorr_Z;
  TH1F *hist_kfactors_ZJets_012j_NLO;

  CLASS_MEMBER(HinvWeights, std::string, kfactor_VBF_zjets_v2_file)
  CLASS_MEMBER(HinvWeights, std::string, kfactor_VBF_wjets_v2_file)
  TFile *kfactor_VBF_zjets_v2_;
  TFile *kfactor_VBF_wjets_v2_;
  TH1F *hist_kfactors_N_W;
  TH1F *hist_kfactors_D_W;
  TH1F *hist_kfactors_N_Z;
  TH1F *hist_kfactors_D_Z;



  TFile *triggerSF_;
  std::vector<TH1F*> hist_trigSF_METL1vec;
  std::vector<TH1F*> hist_trigSF_METHLTvec;
  std::vector<TH1F*> hist_trigSF_MjjHLTvec;
  std::vector<TH1F*> hist_trigSF_JetHLTvec;
  std::vector<TF1*> func_trigSF_METL1vec;
  std::vector<TF1*> func_trigSF_METHLTvec;
  std::vector<TF1*> func_trigSF_MjjHLTvec;
  std::vector<TF1*> func_trigSF_JetHLTvec;

  std::vector<std::vector<std::vector<TF1*> > > func_trigSF_binnedin2d[7];
  std::vector<std::string> errLabel;

  TH1F *hist_trigSF_METL1;
  TH1F *hist_trigSF_METHLT;
  TH1F *hist_trigSF_MjjHLT;
  TH1F *hist_trigSF_JetHLT;

  std::vector<TH3F*> hist_trigSF_3D;

  TH1F *tighteleweight;
  TH1F *tightmuweight;
  TH1F *vetoeleweight;
  TH1F *vetomuweight;

  double lumixsweight;

  std::vector<double> eTight_idisoSF_;
  std::vector<double> eVeto_idisoDataEff_;
  std::vector<double> eVeto_idisoMCEff_;
  std::vector<double> e_gsfidSF_;
  std::vector<double> e_gsfidDataEff_;
  std::vector<double> e_gsfidMCEff_;
  std::vector<double> muTight_idSF_;
  std::vector<double> muTight_isoSF_;
  std::vector<double> muVeto_idDataEff_;
  std::vector<double> muVeto_isoDataEff_;
  std::vector<double> muVeto_idMCEff_;
  std::vector<double> muVeto_isoMCEff_;
  std::vector<double> muTight_idisoSF_;
  std::vector<double> muVeto_idisoDataEff_;
  std::vector<double> muVeto_idisoMCEff_;
  std::vector<double> mu_tkSF_;
  std::vector<double> mu_tkDataEff_;

  std::vector<double> gsf_etabin_;
  std::vector<double> gsf_ptbin_;
  std::vector<double> e_etabin_;
  std::vector<double> e_ptbin_;

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
                  std::vector<double> & aVector,
                  std::vector<double> & ptbin,
                  std::vector<double> & etabin);

  void fillVectorError(const std::string & aFileName, std::vector<double> & aVector, bool upordown);

  double nloReweighting(const double & aMjj, const double & aYstar);


};

}
#undef MEMBER_NP
#endif
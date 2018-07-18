#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvWeights.h"
#include "UserCode/ICHiggsTauTau/interface/TriggerPath.hh"
#include "UserCode/ICHiggsTauTau/interface/TriggerObject.hh"
#include "UserCode/ICHiggsTauTau/interface/EventInfo.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/SimpleParamParser.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/HTFromLHEParticles.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TFile.h"
#include "TROOT.h"
#include "TKey.h"
#include "TCollection.h"
#include "TClass.h"
#include "iostream"
#include "fstream"
#include "sstream"

namespace ic {//namespace

  HinvWeights::HinvWeights(std::string const& name) : ModuleBase(name),
    mc_(mc::summer12_53X),
    era_(era::data_2012_moriond) {
    save_weights_                     = true;
    do_top_reweighting_               = false;
    do_trg_weights_                   = false; //Store as part of total weight

//     //binning in Mjj
//     trigMjjBins_.push_back(0);
//     trigMjjBins_.push_back(800);
//     trigMjjBins_.push_back(1200);
//     trigMjjBins_.push_back(1700);
//     trigMjjBins_.push_back(3000);
// 
//     //binning in Detajj
//     trigDetajjBins_cc_.push_back(1.0);
//     trigDetajjBins_cc_.push_back(1.5);
//     trigDetajjBins_cc_.push_back(2.0);
//     trigDetajjBins_cc_.push_back(2.5);
//     trigDetajjBins_cc_.push_back(3.0);
//     trigDetajjBins_cc_.push_back(3.5);
//     trigDetajjBins_cc_.push_back(4.0);
//     trigDetajjBins_cc_.push_back(5.0);
//     trigDetajjBins_cc_.push_back(10.0);
// 
//     trigDetajjBins_cf_.push_back(3.0);
//     trigDetajjBins_cf_.push_back(3.5);
//     trigDetajjBins_cf_.push_back(4.0);
//     trigDetajjBins_cf_.push_back(5.0);
//     trigDetajjBins_cf_.push_back(6.0);
//     trigDetajjBins_cf_.push_back(10.0);
// 
//     trigDetajjBins_Zmm_cc_.push_back(1.0);
//     trigDetajjBins_Zmm_cc_.push_back(2.0);
//     trigDetajjBins_Zmm_cc_.push_back(3.0);
//     trigDetajjBins_Zmm_cc_.push_back(4.0);
//     trigDetajjBins_Zmm_cc_.push_back(10.0);
// 
//     trigDetajjBins_Zmm_cf_.push_back(3.0);
//     trigDetajjBins_Zmm_cf_.push_back(4.0);
//     trigDetajjBins_Zmm_cf_.push_back(5.0);
//     trigDetajjBins_Zmm_cf_.push_back(10.0);

    trg_applied_in_mc_      = false;
    do_idiso_tight_weights_ = false;
    do_idiso_veto_weights_  = false;
    do_w_soup_              = false;
    do_w_reweighting_       = false;
    do_ewk_w_reweighting_   = false;
    do_dy_soup_             = false;
    do_dy_soup_htbinned_    = false;
    do_dy_reweighting_      = false;
    do_z_reweighting_       = false;
    do_ewk_dy_reweighting_  = false;
    do_prefiring_reweighting_ = false;
    do_lumixs_weights_      = false;
    save_lumixs_weights_    = true;
    input_params_ = "";
    sample_name_= "test";
    input_met_ = "metNoMuons";
    input_jet_ = "pfJetsPFlow";
    mettrg_weight_file_="";
    mettrg_zmm_weight_file_="";

    // For v_nlo_Reweighting (in input/nlo_factors from gitlab)
    kfactors_file_        = "input/nlo_factors/kfactor_24bins.root";
    kfactors_dyjets_file_ = "input/nlo_factors/kfactor_VBF_zll.root";
    kfactors_zjets_file_  = "input/nlo_factors/kfactor_VBF_znn.root";
    kfactors_wjets_file_  = "input/nlo_factors/kfactor_VBF_wjet.root";

    kFactor_ZToNuNu_pT_Mjj_file_="input/nlo_factors/kFactor_ZToNuNu_pT_Mjj.root";
    kFactor_WToLNu_pT_Mjj_file_ ="input/nlo_factors/kFactor_WToLNu_pT_Mjj.root";

    prefiring_etaPt_file_="input/trig_factors/l1v3.root";
  }

  HinvWeights::~HinvWeights() {
    ;
  }

  int HinvWeights::PreAnalysis() {
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "-- PreAnalysis Info for HinvWeights --" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "-- Era: " << Era2String(era_) << std::endl;
    std::cout << "-- MC: " << MC2String(mc_) << std::endl;
    std::cout << "-- Save weights: " << save_weights_ << std::endl;
    std::cout << "-- Do Trg Weights?: \t\t" << do_trg_weights_ << std::endl;
    if(trg_applied_in_mc_) std::cout << "-- Trg Sel Is Applied In MC make sure you use DataMC scale factors for trigger weightst" << std::endl;
    else std::cout << "-- Trg Sel Not Applied In MC make sure you use raw Data efficiencies for trigger weightst" << std::endl;
    //std::cout << "Trg Sel Applied?: \t\t" << trg_applied_in_mc_ << std::endl;
    std::cout << "-- Do ID & iso weights for Tight leptons ?: \t\t" << do_idiso_tight_weights_ << std::endl;
    std::cout << "-- Do ID & iso weights for veto leptons ?: \t\t" << do_idiso_veto_weights_ << std::endl;

    std::cout << "-- Input MET for MET HLT:  \t\t" << input_met_ << std::endl;
    std::cout << "-- Note: Input MET for MET L1 is always metNoMuons." << std::endl;

    //Making output histograms
    TFileDirectory const& dir = fs_->mkdir("leptoneffweights");
    tighteleweight = dir.make<TH1F>("tighteleweight","tighteleweight",2000,0.,2.);
    tightmuweight  = dir.make<TH1F>("tightmuweight","tightmuweight",2000,0.,2.);
    vetoeleweight = dir.make<TH1F>("vetoeleweight","vetoeleweight",2000,0.,2.);
    vetomuweight  = dir.make<TH1F>("vetomuweight","vetomuweight",2000,0.,2.);

    if(do_lumixs_weights_){
      //Get Sample name file
      std::string suffix=".root";
      std::size_t found = sample_name_.find(suffix);
      if(found==std::string::npos){
        lumixsweight=1;
        std::cout<<"Non-standard sample name format not doing lumixs weight"<<std::endl;
      }
      else{
        //std::cout << sample_name_ << " " << sample_name_.find("split") <<std::endl;
        if (sample_name_.find("split")!=sample_name_.npos){
          sample_name_.erase(found-8,13);
        }
        else sample_name_.erase(found,5);
        std::cout << "Sample Name: "<<sample_name_<<std::endl;

        //Get lumi xs and events from params file
        SimpleParamParser parser;
        std::cout << "** Parsing parameter file... **" << input_params_ << std::endl;
        parser.ParseFile(input_params_);
        std::cout<<"-- Parsed"<<std::endl;
        double xs=parser.GetParam<double>("XS_"+sample_name_);
        std::cout<<"-- Got xs"<<std::endl;
        double events=parser.GetParam<double>("EVT_"+sample_name_);
        std::cout<<"-- Got events"<<std::endl;
        double lumi=parser.GetParam<double>("LUMI_DATA");
        std::cout<<"-- Got lumi"<<std::endl;

        std::cout<<"-- XS is: "<<xs<<"pb"<<std::endl;
        std::cout<<"-- EVT is: "<<events<<std::endl;
        std::cout<<"-- LUMI is: "<<lumi<<"pb^-1"<<std::endl;

        lumixsweight=xs*lumi/events;
        std::cout<<"-- LUMIXSWEIGHT is: "<<lumixsweight<<std::endl;
      }
    }

    if (do_w_soup_) {
      if (mc_ == mc::spring16_80X || mc_ == mc::summer16_80X){
        std::cout << "-- Making W Soup:" << std::endl;
        std::cout << "nInc = " << n_inc_ << std::endl;
        w1_ = (n_inc_*f1_) / ( (n_inc_*f1_) + n1_ );
        w2_ = (n_inc_*f2_) / ( (n_inc_*f2_) + n2_ );
        w3_ = (n_inc_*f3_) / ( (n_inc_*f3_) + n3_ );
        w4_ = (n_inc_*f4_) / ( (n_inc_*f4_) + n4_ );
        w5_ = (n_inc_*f5_) / ( (n_inc_*f5_) + n5_ );
        w6_ = (n_inc_*f6_) / ( (n_inc_*f6_) + n6_ );
        w7_ = (n_inc_*f7_) / ( (n_inc_*f7_) + n7_ );
        std::cout << "f1 = " << f1_ << "\t" << "n1 = " << n1_ << "\t" << "w1 = " << w1_ << std::endl;
        std::cout << "f2 = " << f2_ << "\t" << "n2 = " << n2_ << "\t" << "w2 = " << w2_ << std::endl;
        std::cout << "f3 = " << f3_ << "\t" << "n3 = " << n3_ << "\t" << "w3 = " << w3_ << std::endl;
        std::cout << "f4 = " << f4_ << "\t" << "n4 = " << n4_ << "\t" << "w4 = " << w4_ << std::endl;
        std::cout << "f5 = " << f5_ << "\t" << "n5 = " << n5_ << "\t" << "w5 = " << w5_ << std::endl;
        std::cout << "f6 = " << f6_ << "\t" << "n6 = " << n6_ << "\t" << "w6 = " << w6_ << std::endl;
        std::cout << "f7 = " << f7_ << "\t" << "n7 = " << n7_ << "\t" << "w7 = " << w7_ << std::endl;
      }
      else{
      std::cout << "Making W Soup:" << std::endl;
      std::cout << "nInc = " << n_inc_ << std::endl;
      w1_ = (n_inc_*f1_) / ( (n_inc_*f1_) + n1_ );
      w2_ = (n_inc_*f2_) / ( (n_inc_*f2_) + n2_ );
      w3_ = (n_inc_*f3_) / ( (n_inc_*f3_) + n3_ );
      w4_ = (n_inc_*f4_) / ( (n_inc_*f4_) + n4_ );
      std::cout << "f1 = " << f1_ << "\t" << "n1 = " << n1_ << "\t" << "w1 = " << w1_ << std::endl;
      std::cout << "f2 = " << f2_ << "\t" << "n2 = " << n2_ << "\t" << "w2 = " << w2_ << std::endl;
      std::cout << "f3 = " << f3_ << "\t" << "n3 = " << n3_ << "\t" << "w3 = " << w3_ << std::endl;
      std::cout << "f4 = " << f4_ << "\t" << "n4 = " << n4_ << "\t" << "w4 = " << w4_ << std::endl;
      }
    }
    if (do_dy_soup_ || do_dy_soup_htbinned_) {
      std::cout << "-- Making DY Soup:" << std::endl;
      std::cout << "nInc = " << zn_inc_ << std::endl;
      zw1_ = (zn_inc_*zf1_) / ( (zn_inc_*zf1_) + zn1_ );
      zw2_ = (zn_inc_*zf2_) / ( (zn_inc_*zf2_) + zn2_ );
      zw3_ = (zn_inc_*zf3_) / ( (zn_inc_*zf3_) + zn3_ );
      zw4_ = (zn_inc_*zf4_) / ( (zn_inc_*zf4_) + zn4_ );
      std::cout << "f1 = " << zf1_ << "\t" << "n1 = " << zn1_ << "\t" << "w1 = " << zw1_ << std::endl;
      std::cout << "f2 = " << zf2_ << "\t" << "n2 = " << zn2_ << "\t" << "w2 = " << zw2_ << std::endl;
      std::cout << "f3 = " << zf3_ << "\t" << "n3 = " << zn3_ << "\t" << "w3 = " << zw3_ << std::endl;
      std::cout << "f4 = " << zf4_ << "\t" << "n4 = " << zn4_ << "\t" << "w4 = " << zw4_ << std::endl;
    }
    if (do_w_reweighting_ || do_dy_reweighting_ || do_z_reweighting_) { // For v_nlo_Reweighting (kfactors.root file in input/scalefactors from MIT group)

      kfactors_wjets_  = TFile::Open(kfactors_wjets_file_.c_str());
      kfactors_zjets_  = TFile::Open(kfactors_zjets_file_.c_str());
      kfactors_dyjets_ = TFile::Open(kfactors_dyjets_file_.c_str());
      kfactors_        = TFile::Open(kfactors_file_.c_str());


      if (do_w_reweighting_) {
        std::cout << " -- Applying reweighting of W events to NLO from MIT (Raffaele)." << std::endl;
        hist_kfactors_EWKcorr_W      = (TH1F*)kfactors_->Get("EWKcorr/W");
        hist_kfactors_WJets_012j_NLO = (TH1F*)kfactors_->Get("WJets_012j_NLO/nominal");
        hist_kfactors_vbf_cnc_W      = (TH1F*)kfactors_wjets_->Get("kfactors_cc/kfactor_vbf");
      }
      if (do_dy_reweighting_) {
        std::cout << " -- Applying reweighting of DY events to NLO from MIT (Raffaele)." << std::endl;
        hist_kfactors_EWKcorr_Z      = (TH1F*)kfactors_->Get("EWKcorr/Z");
        hist_kfactors_ZJets_012j_NLO = (TH1F*)kfactors_->Get("ZJets_012j_NLO/nominal");
        hist_kfactors_vbf_cnc_DY     = (TH1F*)kfactors_dyjets_->Get("kfactors_cc/kfactor_vbf");
      }
      if (do_z_reweighting_) {
        std::cout << " -- Applying reweighting of Zvv events to NLO from MIT (Raffaele)." << std::endl;
        hist_kfactors_EWKcorr_Z      = (TH1F*)kfactors_->Get("EWKcorr/Z");
        hist_kfactors_ZJets_012j_NLO = (TH1F*)kfactors_->Get("ZJets_012j_NLO/nominal");
        hist_kfactors_vbf_cnc_Z      = (TH1F*)kfactors_zjets_->Get("kfactors_cc/kfactor_vbf");
      }
    }

    if (do_ewk_w_reweighting_ || do_ewk_dy_reweighting_) {
      kFactor_ZToNuNu_pT_Mjj_ = TFile::Open(kFactor_ZToNuNu_pT_Mjj_file_.c_str());
      kFactor_WToLNu_pT_Mjj_ = TFile::Open(kFactor_WToLNu_pT_Mjj_file_.c_str());

      if (!kFactor_ZToNuNu_pT_Mjj_) return 1;
      if (!kFactor_WToLNu_pT_Mjj_) return 1;

      if (do_ewk_w_reweighting_) {
        std::cout << " -- Applying reweighting of W events to NLO from MIT (Raffaele)." << std::endl;
        hist_kFactors_ewk_W = (TH2F*)kFactor_WToLNu_pT_Mjj_->Get("TH2F_kFactor");
      }
      if (do_ewk_dy_reweighting_) {
        std::cout << " -- Applying reweighting of DY events to NLO from MIT (Raffaele)." << std::endl;
        hist_kFactors_ewk_Z = (TH2F*)kFactor_ZToNuNu_pT_Mjj_->Get("TH2F_kFactor");
      }
    }

    if (do_prefiring_reweighting_) {
      prefiring_etaPt_ = TFile::Open(prefiring_etaPt_file_.c_str());

      if (!prefiring_etaPt_) return 1;

      std::cout << " -- Applying reweighting signal events (VBF and ggH) to account for pre-firing." << std::endl;
      hist_prefiring_etaPt_SingleMuon = (TH2F*)prefiring_etaPt_->Get("h_SingleMuon_spike_finor_pteta");
      hist_prefiring_etaPt_JetHT      = (TH2F*)prefiring_etaPt_->Get("h_JetHT_spike_finor_pteta");
    }

    if (save_weights_ && do_trg_weights_){
      // do weights even if not applied, to fill histo with weight for comparison !

      //get trigger scale factor histograms from file
      mettrigSF_ = new TFile(mettrg_weight_file_.c_str());
      if (!mettrigSF_) return 1;
      mettrigZmmSF_ = new TFile(mettrg_zmm_weight_file_.c_str());
      if (!mettrigZmmSF_) return 1;
      std::cout<<"Getting trigger efficiency functions"<<std::endl;

      mettrigSF_->cd();
      TIter next(mettrigSF_->GetListOfKeys());
      TKey *key;
      TF1 *ftmp = 0;
      while ((key = (TKey*)next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TF1")) continue;
        TF1 *tmp = (TF1*)key->ReadObj();
//         //binning in Mjj
//         for(unsigned iVar1=0;iVar1<(trigMjjBins_.size()-1);iVar1++){
//           std::ostringstream flabel;
//           flabel << "fitfunc_vbf_mjj_"
//           << trigMjjBins_[iVar1]
//           << ".000000_"
//           << trigMjjBins_[iVar1+1]
//           << ".000000";
// 
//           std::string str1 = flabel.str();
//           std::string str2 = tmp->GetName();
//           if (str1==str2){
//             ftmp[iVar1] = tmp;
//             std::cout << " -- trigger Function label " << flabel.str() << " bin " << iVar1 << std::endl;
//             break;
//           }
//         }

//         //binning in Detajj cc
//         for(unsigned iVar1=0;iVar1<(trigDetajjBins_cc_.size()-1);iVar1++){
//           std::ostringstream flabel_cc;
//           flabel_cc << "fitfunc_recoil_cc_vs_detajj_"
//                     << std::setprecision(1) << std::fixed << trigDetajjBins_cc_[iVar1]
//                     << "_"
//                     << trigDetajjBins_cc_[iVar1+1];
// 
//           std::string str1 = flabel_cc.str();
//           std::string str2 = tmp->GetName();
//           if (str1==str2){
//             ftmp[0][iVar1] = tmp;
//             std::cout << " -- trigger Function label " << flabel_cc.str() << " bin " << iVar1 << std::endl;
//             break;
//           }
//         }
//         //binning in Detajj cf
//         for(unsigned iVar1=0;iVar1<(trigDetajjBins_cf_.size()-1);iVar1++){
//           std::ostringstream flabel_cf;
//           flabel_cf << "fitfunc_recoil_cf_vs_detajj_"
//                     << std::setprecision(1) << std::fixed << trigDetajjBins_cf_[iVar1]
//                     << "_"
//                     << trigDetajjBins_cf_[iVar1+1];
// 
//           std::string str1 = flabel_cf.str();
//           std::string str2 = tmp->GetName();
//           if (str1==str2){
//             ftmp[1][iVar1] = tmp;
//             std::cout << " -- trigger Function label " << flabel_cf.str() << " bin " << iVar1 << std::endl;
//             break;
//           }
//         }
        std::ostringstream flabel;
        flabel << "sigmoidFit";

        std::string str1 = flabel.str();
        std::string str2 = tmp->GetName();
        if (str1==str2){
          ftmp = tmp;
          std::cout << " -- trigger Function label " << flabel.str() << std::endl;
          break;
        }
      }


      mettrigZmmSF_->cd();
      TIter nextZmm(mettrigZmmSF_->GetListOfKeys());
      TF1 *ftmp2 = 0;

      while ((key = (TKey*)nextZmm())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TF1")) continue;
        TF1 *tmp = (TF1*)key->ReadObj();
//         //binning in Mjj
//         for(unsigned iVar1=0;iVar1<(trigMjjBins_.size()-1);iVar1++){
//           std::ostringstream flabel;
//           flabel << "fitfunc_vbf_mjj_"
//           << trigMjjBins_[iVar1]
//           << ".000000_"
//           << trigMjjBins_[iVar1+1]
//           << ".000000";
// 
//           std::string str1 = flabel.str();
//           std::string str2 = tmp->GetName();
//           if (str1==str2){
//             ftmp2[iVar1] = tmp;
//             std::cout << " -- Zmm trigger Function label " << flabel.str()  << " bin " << iVar1 << std::endl;
//             break;
//           }
//         }

//         //binning in Detajj cc
//         for( unsigned iVar1=0; iVar1<(trigDetajjBins_Zmm_cc_.size()-1); iVar1++ ){
//           std::ostringstream flabel_cc;
//           flabel_cc << "fitfunc_recoil_cc_vs_detajj_"
//                     << std::setprecision(1) << std::fixed << trigDetajjBins_Zmm_cc_[iVar1]
//                     << "_"
//                     << trigDetajjBins_Zmm_cc_[iVar1+1];
// 
//           std::string str1 = flabel_cc.str();
//           std::string str2 = tmp->GetName();
//           if (str1==str2){
//             ftmp2[0][iVar1] = tmp;
//             std::cout << " -- Zmm trigger Function label " << flabel_cc.str() << " bin " << iVar1 << std::endl;
//             break;
//           }
//         }
//         //binning in Detajj cf
//         for(unsigned iVar1=0;iVar1<(trigDetajjBins_Zmm_cf_.size()-1);iVar1++){
//           std::ostringstream flabel_cf;
//           flabel_cf << "fitfunc_recoil_cf_vs_detajj_"
//                     << std::setprecision(1) << std::fixed << trigDetajjBins_Zmm_cf_[iVar1]
//                     << "_"
//                     << trigDetajjBins_Zmm_cf_[iVar1+1];
// 
//           std::string str1 = flabel_cf.str();
//           std::string str2 = tmp->GetName();
//           if (str1==str2){
//             ftmp2[1][iVar1] = tmp;
//             std::cout << " -- Zmm trigger Function label " << flabel_cf.str() << " bin " << iVar1 << std::endl;
//             break;
//           }
//         }
        std::ostringstream flabel;
        flabel << "sigmoidFit";

        std::string str1 = flabel.str();
        std::string str2 = tmp->GetName();
        if (str1==str2){
          ftmp2 = tmp;
          std::cout << " -- Zmm trigger Function label " << flabel.str() << std::endl;
          break;
        }
      }
//       //binning in Mjj
//       for(unsigned iVar1=0;iVar1<(trigMjjBins_.size()-1);iVar1++){
//         func_mettrigSF_.push_back(ftmp[iVar1]);
//         func_mettrigZmmSF_.push_back(ftmp2[iVar1]);
//       }
//       //binning in Detajj cc
//       for( unsigned iVar1=0; iVar1<(trigDetajjBins_cc_.size()-1); iVar1++ ){
//         func_mettrigSF_[0].push_back(ftmp[0][iVar1]);
//       }
//       for( unsigned iVar1=0; iVar1<(trigDetajjBins_Zmm_cc_.size()-1); iVar1++ ){
//         func_mettrigZmmSF_[0].push_back(ftmp2[0][iVar1]);
//       }
//       //binning in Detajj cf
//       for( unsigned iVar1=0; iVar1<(trigDetajjBins_cf_.size()-1); iVar1++ ){
//         func_mettrigSF_[1].push_back(ftmp[1][iVar1]);
//       }
//       for( unsigned iVar1=0; iVar1<(trigDetajjBins_Zmm_cf_.size()-1); iVar1++ ){
//         func_mettrigZmmSF_[1].push_back(ftmp2[1][iVar1]);
//       }
      func_mettrigSF_.push_back(ftmp);
      func_mettrigZmmSF_.push_back(ftmp2);
      std::cout<<"-- Done!"<<std::endl;
    }
    if (save_weights_){
      std::vector<double> dummypt;
      std::vector<double> dummyeta;
      //last bool: protect against values > 1, just for efficiencies not for SF
      fillVector("input/scale_factors/Summer16_80X_ele_tight_id_SF.txt",6,10,eTight_idisoSF_,e_ptbin_,e_etabin_,false);
      fillVector("input/scale_factors/Summer16_80X_ele_trig_data_eff.txt",15,6,e_trigDataEff_,e_pttrig_,e_etatrig_,true);
      fillVector("input/scale_factors/Summer16_80X_gsf_id_SF.txt",3,30,e_gsfidSF_,gsf_ptbin_,gsf_etabin_,false);
      fillVector("input/scale_factors/Summer16_80X_mu_tight_id_SF.txt",6,8,muTight_idSF_,mu_ptbin_,mu_etabin_,false);
      fillVector("input/scale_factors/Summer16_80X_mu_trackingSF_smooth.txt",1,23,mu_tkSF_,tk_ptbin_,tk_etabin_,true);

      fillVector("input/scale_factors/Summer16_80X_ele_veto_id_data_eff.txt",6,10,eVeto_idisoDataEff_,dummypt,dummyeta,true);
      fillVector("input/scale_factors/Summer16_80X_ele_veto_id_mc_eff.txt",6,10,eVeto_idisoMCEff_,dummypt,dummyeta,true);
      fillVector("input/scale_factors/Summer16_80X_gsf_id_data_eff.txt",3,30,e_gsfidDataEff_,dummypt,dummyeta,true);
      fillVector("input/scale_factors/Summer16_80X_gsf_id_mc_eff.txt",3,30,e_gsfidMCEff_,dummypt,dummyeta,true);

      fillVector("input/scale_factors/Summer16_80X_mu_tight_iso_SF.txt",6,8,muTight_isoSF_,dummypt,dummyeta,false);
      fillVector("input/scale_factors/Summer16_80X_mu_loose_id_SF.txt",6,8,muVeto_idSF_,dummypt,dummyeta,false);
      fillVector("input/scale_factors/Summer16_80X_mu_loose_iso_SF.txt",6,8,muVeto_isoSF_,dummypt,dummyeta,false);
      fillVector("input/scale_factors/Summer16_80X_mu_loose_id_data_eff.txt",6,8,muVeto_idDataEff_,dummypt,dummyeta,true);
      fillVector("input/scale_factors/Summer16_80X_mu_loose_iso_data_eff.txt",6,8,muVeto_isoDataEff_,dummypt,dummyeta,true);
      fillVector("input/scale_factors/Summer16_80X_mu_loose_id_mc_eff.txt",6,8,muVeto_idMCEff_,dummypt,dummyeta,true);
      fillVector("input/scale_factors/Summer16_80X_mu_loose_iso_mc_eff.txt",6,8,muVeto_isoMCEff_,dummypt,dummyeta,true);

      /*for (unsigned iBin(0); iBin<muTight_idSF_[0].size();++iBin){
        muTight_idisoSF_[0].push_back(muTight_idSF_[0][iBin]*muTight_isoSF_[0][iBin]);
        muTight_idisoSF_[1].push_back(muTight_idSF_[1][iBin]*muTight_isoSF_[1][iBin]);
        muTight_idisoSF_[2].push_back(muTight_idSF_[2][iBin]*muTight_isoSF_[2][iBin]);
        muVeto_idisoDataEff_[0].push_back(muVeto_idDataEff_[0][iBin]*muVeto_isoDataEff_[0][iBin]);
        muVeto_idisoMCEff_[0].push_back(muVeto_idMCEff_[0][iBin]*muVeto_isoMCEff_[0][iBin]);
        muVeto_idisoDataEff_[1].push_back(muVeto_idDataEff_[1][iBin]*muVeto_isoDataEff_[1][iBin]);
        muVeto_idisoMCEff_[1].push_back(muVeto_idMCEff_[1][iBin]*muVeto_isoMCEff_[1][iBin]);
        muVeto_idisoDataEff_[2].push_back(muVeto_idDataEff_[2][iBin]*muVeto_isoDataEff_[2][iBin]);
        muVeto_idisoMCEff_[2].push_back(muVeto_idMCEff_[2][iBin]*muVeto_isoMCEff_[2][iBin]);
        //std::cout<<muVeto_idisoMCEff_.back()<<" "<<muVeto_idisoDataEff_.back()<<" "<<muTight_idisoSF_.back()<<std::endl;//!!
	}*/

      eventsWithGenElectron_ = 0;
      eventsWithGenElectronFromTau_ = 0;
      eventsWithGenMuon_ = 0;
      eventsWithGenMuonFromTau_ = 0;
      eventsWithGenElectronInAcc_ = 0;
      eventsWithGenElectronFromTauInAcc_ = 0;
      eventsWithGenMuonInAcc_ = 0;
      eventsWithGenMuonFromTauInAcc_ = 0;
      eventsWithGenTau_ = 0;
    }

    return 0;
  }

  int HinvWeights::Execute(TreeEvent *event) {
    EventInfo * eventInfo = event->GetPtr<EventInfo>("eventInfo");

    if(do_lumixs_weights_){
      //std::cout<<"weight before lumixs: "<<eventInfo->total_weight()<<std::endl;
      //std::cout<<"intended lumixsweight: "<<lumixsweight<<std::endl;
      eventInfo->set_weight("lumixs",lumixsweight);
      //std::cout<<"weight after lumixs: "<<eventInfo->total_weight()<<std::endl;
      //std::cout<<"lumixsweight info "<<eventInfo->weight_is_enabled("lumixs") << " " << eventInfo->weight_defined("lumixs") << " " << eventInfo->weight("lumixs") << std::endl;
    }
    //else eventInfo->set_weight("!lumixs",lumixsweight);

    if (save_weights_){//Save weights

          //double weight = 1.0;

 //    if (do_btag_weight_) {
//       std::vector<PFJet*> jets = event->GetPtrVec<PFJet>(input_jet_); // Make a copy of the jet collection
//       ic::erase_if(jets,!boost::bind(MinPtMaxEta, _1, 20.0, 2.4));
//       //double no_btag_weight = btag_weight.GetWeight(jets, "CSVM", 0, 0, is_2012_);
//       //double inclusive_btag_weight = btag_weight.GetWeight(jets, "CSVM", 1, 99, is_2012_);
//       double no_btag_weight = 1.0;
//       double inclusive_btag_weight = 1.0;
//       btag_weight.ReTag(jets, mc_ == mc::summer12_53X);
//       event->Add("no_btag_weight", no_btag_weight);
//       event->Add("inclusive_btag_weight", inclusive_btag_weight);
//     }

//     //get METnoMuons:
//     Met const* metHLT = event->GetPtr<Met>(input_met_);
//     Met const* metL1 = event->GetPtr<Met>("metNoMuons");
//     double l1met = metL1->pt();
//     double hltmet = metHLT->pt();

    std::vector<PFJet*> & jets = event->GetPtrVec<PFJet>("pfJetsPFlow");
    std::sort(jets.begin(), jets.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
    ROOT::Math::PtEtaPhiEVector mhtVec30_eta3d0(0,0,0,0);
    double mht30_eta3d0 = -1;

    for (unsigned i = 0; i < jets.size(); ++i) {//Loop on jets
      ROOT::Math::PtEtaPhiEVector jetvec = jets[i]->vector();
      if(jets[i]->pt()>30 && fabs(jets[i]->eta())<3.0) mhtVec30_eta3d0 += jetvec;
    }//endof Loop on jets
    mht30_eta3d0 = mhtVec30_eta3d0.Et();
//     std::cout << " ***---*** mht30_eta3d0 : " << mht30_eta3d0 << std::endl;
//     if(mht30_eta3d0<=0){
//       std::cout << " ---*** Error: mht30_eta3d0 has not been filled: " << mht30_eta3d0 << std::endl;
//       return 1;
//     }

    if (do_trg_weights_){//do_trg_weights_
//       double mjj=0.;
//       double detajj=-1;
//       double jet1pt=0.;
//       double jet2pt=0.;
//       bool hasJetsInHF = false;

//       //get 2 leading jets
//       std::vector<CompositeCandidate *> const& dijet_vec = event->GetPtrVec<CompositeCandidate>("jjLeadingCandidates");
//       if (dijet_vec.size() > 0) {//if dijets

//         CompositeCandidate const* dijet = dijet_vec.at(0);
// 
//         Candidate const* jet1 = dijet->GetCandidate("jet1");
//         Candidate const* jet2 = dijet->GetCandidate("jet2");
// 
//         mjj = dijet->M();
//         detajj = fabs(jet1->eta() - jet2->eta());
//         jet1pt = jet1->pt();
//         jet2pt = jet2->pt();
//         hasJetsInHF = fabs(jet1->eta())>=3 || fabs(jet2->eta())>=3 ;
//         std::cout<<"mjj "<<mjj<<" j2pt "<<jet2pt<<" metl1 "<<l1met<<" hltmet "<<hltmet<<std::endl;
//         if(l1met!=hltmet){
//           std::cout<<"Error: you must use metnomuons for both l1met and hltmet"<<std::endl;
//           return 1;
//         }

//         //FIND WHICH BIN YOU'RE IN
//         int var1bin=-10;
//         if(mjj<trigMjjBins_[0])var1bin=-1;
//         else{
//           for(unsigned iBin=0;iBin<(trigMjjBins_.size()-1);iBin++){
//             if(mjj<trigMjjBins_[iBin+1]){
//               var1bin=iBin;
//               break;
//             }
//           }
//           if(var1bin==-10)var1bin=trigMjjBins_.size()-2;
//         }
// 
// 
//         double trgweight = 0;
//         double trgweightzmm = 0;
//         double xmin,xmax;
//         if(var1bin!=-1){
// 
//           func_mettrigSF_[var1bin]->GetRange(xmin,xmax);
// 
//           if(hltmet<=xmax){
//             if(hltmet>=xmin){
//               trgweight = func_mettrigSF_[var1bin]->Eval(hltmet);
//               trgweightzmm = func_mettrigZmmSF_[var1bin]->Eval(hltmet);
//             }
//           }
//           else {
//             trgweight = func_mettrigSF_[var1bin]->Eval(xmax);
//             trgweightzmm = func_mettrigZmmSF_[var1bin]->Eval(xmax);
//           }
//         }
//         std::cout << " -- Mjj " << mjj << " bin " << var1bin << " " << trigMjjBins_[var1bin] << "-" << trigMjjBins_[var1bin+1]
//             << " weights " << func_mettrigSF_[var1bin]->Eval(260.) << " " << func_mettrigZmmSF_[var1bin]->Eval(260.)
//             << std::endl;

//         double trgweight = 0;
//         double trgweightzmm = 0;

//         if (fabs(jet1->eta())<3 && fabs(jet2->eta())<3) {
//           //FIND WHICH BIN YOU'RE IN
//           int var1bin=-10;
//           if(detajj<trigDetajjBins_cc_[0]){
//             var1bin=-1;
//           }
//           else{
//             for(unsigned iBin=0;iBin<(trigDetajjBins_cc_.size()-1);iBin++){
//               if(detajj<trigDetajjBins_cc_[iBin+1]){
//                 var1bin=iBin;
//                 break;
//               }
//             }
//             if(var1bin==-10)var1bin=trigDetajjBins_cc_.size()-2;
//           }
// 
//           int var1binzmm=-10;
//           if(detajj<trigDetajjBins_Zmm_cc_[0]){
//             var1binzmm=-1;
//           }
//           else{
//             for(unsigned iBin=0;iBin<(trigDetajjBins_Zmm_cc_.size()-1);iBin++){
//               if(detajj<trigDetajjBins_Zmm_cc_[iBin+1]){
//                 var1binzmm=iBin;
//                 break;
//               }
//             }
//             if(var1binzmm==-10)var1binzmm=trigDetajjBins_Zmm_cc_.size()-2;
//           }
// 
//           double xmin,xmax;
//           if(var1bin!=-1){
// 
//             func_mettrigSF_[0].at(var1bin)->GetRange(xmin,xmax);
// 
//             if(hltmet<=xmax){
//               if(hltmet>=xmin){
//                   trgweight = func_mettrigSF_[0].at(var1bin)->Eval(hltmet);
//               }
//             }
//             else {
//                 trgweight = func_mettrigSF_[0].at(var1bin)->Eval(xmax);
//             }
//           }
// 
//           double xminzmm,xmaxzmm;
//           if(var1binzmm!=-1){
// 
//             func_mettrigZmmSF_[0][var1binzmm]->GetRange(xminzmm,xmaxzmm);
// 
//             if(hltmet<=xmaxzmm){
//               if(hltmet>=xminzmm){
//                 trgweightzmm = func_mettrigZmmSF_[0][var1binzmm]->Eval(hltmet);
//               }
//             }
//             else {
//               trgweightzmm = func_mettrigZmmSF_[0][var1binzmm]->Eval(xmaxzmm);
//             }
//           }
//         } else if ((fabs(jet1->eta())<3 && fabs(jet2->eta())>=3)||(fabs(jet1->eta())>=3 && fabs(jet2->eta())<3)) {
//           //FIND WHICH BIN YOU'RE IN
//           int var1bin=-10;
//           if(detajj<trigDetajjBins_cf_[0]){
//             var1bin=-1;
//           }
//           else{
//             for(unsigned iBin=0;iBin<(trigDetajjBins_cf_.size()-1);iBin++){
//               if(detajj<trigDetajjBins_cf_[iBin+1]){
//                 var1bin=iBin;
//                 break;
//               }
//             }
//             if(var1bin==-10)var1bin=trigDetajjBins_cf_.size()-2;
//           }
// 
//           int var1binzmm=-10;
//           if(detajj<trigDetajjBins_Zmm_cf_[0]) {
//             var1binzmm=-1;
//           }
//           else{
//             for(unsigned iBin=0;iBin<(trigDetajjBins_Zmm_cf_.size()-1);iBin++){
//               if(detajj<trigDetajjBins_Zmm_cf_[iBin+1]){
//                 var1binzmm=iBin;
//                 break;
//               }
//             }
//             if(var1binzmm==-10)var1binzmm=trigDetajjBins_Zmm_cf_.size()-2;
//           }
// 
//           double xmin,xmax;
//           if(var1bin!=-1){
// 
//             func_mettrigSF_[1][var1bin]->GetRange(xmin,xmax);
// 
//             if(hltmet<=xmax){
//               if(hltmet>=xmin){
//                 trgweight = func_mettrigSF_[1][var1bin]->Eval(hltmet);
//               }
//             }
//             else {
//               trgweight = func_mettrigSF_[1][var1bin]->Eval(xmax);
//             }
//           }
// 
//           double xminzmm,xmaxzmm;
//           if(var1binzmm!=-1){
// 
//             func_mettrigZmmSF_[1][var1binzmm]->GetRange(xminzmm,xmaxzmm);
// 
//             if(hltmet<=xmaxzmm){
//               if(hltmet>=xminzmm){
//                 trgweightzmm = func_mettrigZmmSF_[1][var1binzmm]->Eval(hltmet);
//               }
//             }
//             else {
//               trgweightzmm = func_mettrigZmmSF_[1][var1binzmm]->Eval(xmaxzmm);
//             }
//           }
//         }

      if (mht30_eta3d0>0.) {
        double trgweight = 0;
        double trgweightzmm = 0;

        double xmin,xmax;
        func_mettrigSF_[0]->GetRange(xmin,xmax);
        if(mht30_eta3d0<=xmax){
          if(mht30_eta3d0>=xmin){
            trgweight = func_mettrigSF_[0]->Eval(mht30_eta3d0);
          }
        } else {
          trgweight = func_mettrigSF_[0]->Eval(xmax);
        }

        double xminzmm,xmaxzmm;
        func_mettrigZmmSF_[0]->GetRange(xminzmm,xmaxzmm);
        if(mht30_eta3d0<=xmaxzmm){
          if(mht30_eta3d0>=xminzmm){
            trgweightzmm = func_mettrigZmmSF_[0]->Eval(mht30_eta3d0);
          }
        } else {
          trgweightzmm = func_mettrigZmmSF_[0]->Eval(xmaxzmm);
        }
        //SET TRIGGER WEIGHT
        eventInfo->set_weight("!mettrigSF",trgweight);
        eventInfo->set_weight("!mettrigZmmSF",trgweightzmm);
      }

//       }//endof if dijets //dijet pair
    }//do trig weights
    //else {
    //std::cout << " skipping trigger stuff" << std::endl;
    //}
    //eventInfo->set_weight("!trigger",metl1weight*methltweight*mjjhltweight*jet1hltweight*jet2hltweight);

    //eventInfo->set_weight("lepton", weight);
    if (do_top_reweighting_) {
      double top_wt = 1.0;
      double top_wt_up = 1.0;
      double top_wt_down = 1.0;
      std::vector<GenParticle *> const& parts = event->GetPtrVec<GenParticle>("genParticles");
      for (unsigned i = 0; i < parts.size(); ++i) {
        if (parts[i]->status() == 3 && abs(parts[i]->pdgid()) == 6) {
          double pt = parts[i]->pt();
          pt = std::min(pt, 400.);
          top_wt *= std::exp(0.156-0.00137*pt);
        }
      }
      top_wt = std::sqrt(top_wt);
      top_wt_up = top_wt * top_wt;
      top_wt_down = 1.0;
      eventInfo->set_weight("!tquark_weight_up", top_wt_up);
      eventInfo->set_weight("!tquark_weight_down", top_wt_down);
      eventInfo->set_weight("tquark_weight", top_wt);

      //madgraph ttbar reweighting
      std::size_t foundttbar = sample_name_.find("TTJets");
      if(foundttbar!=std::string::npos){
        int ngenLeptonsStatus3=0;
        double genWeight=1;
        for (unsigned i = 0; i < parts.size(); ++i) {
          if (parts[i]->status() == 3 && ((abs(parts[i]->pdgid()) == 11||(abs(parts[i]->pdgid()) == 13)||(abs(parts[i]->pdgid()) == 15)))) {
            ngenLeptonsStatus3++;
          }
        }
        if(ngenLeptonsStatus3==2) { genWeight=pow(0.1086/(1./9.),2); }
        else if(ngenLeptonsStatus3==1) { genWeight=(0.1086/(1./9.))*(0.6741/(2./3.)); }
        else { genWeight=pow(0.6741/(2./3.),2); }
        eventInfo->set_weight("madgraph_ttbarbr_weight",genWeight);
      }
    }//do top
    //else {
    //std::cout << " skipping top stuff" << std::endl;
    //}

    //ID+iso tight leptons
    std::vector<Electron*> const& elecs = event->GetPtrVec<Electron>("selElectrons");
    double ele_weight[3] = {1.0,1.0,1.0};
    double gsf_weight[3] = {1.0,1.0,1.0};
    //record first two electrons
    double trigW[3][2];
    for (unsigned err(0); err<3;++err){
      trigW[err][0] = 1.;
      trigW[err][1] = 1.;
    }
    for (unsigned iEle(0); iEle<elecs.size();++iEle){
      for (unsigned err(0); err<3;++err){
	ele_weight[err] *= eTight_idisoSF_[err][findPtEtaBin(elecs[iEle]->pt(),elecs[iEle]->eta(),e_ptbin_,e_etabin_)];
	//ele_weight[err] *= e_gsfidSF_[err][findPtEtaBin(elecs[iEle]->pt(),elecs[iEle]->eta(),gsf_ptbin_,gsf_etabin_)];
	gsf_weight[err] *= e_gsfidSF_[err][findPtEtaBin(elecs[iEle]->pt(),elecs[iEle]->eta(),gsf_ptbin_,gsf_etabin_)];
	if (iEle<2) trigW[err][iEle] = e_trigDataEff_[err][findPtEtaBin(elecs[iEle]->pt(),fabs(elecs[iEle]->eta()),e_pttrig_,e_etatrig_)];
      }
    }

    //calculate electron trigger weight
    double eleTrigW[3] = {1,1,1};
    for (unsigned err(0); err<3;++err){
      if (elecs.size()>1) eleTrigW[err] = trigW[err][0]+trigW[err][1]-trigW[err][0]*trigW[err][1];
      else eleTrigW[err] = trigW[err][0]; 
    }
    eventInfo->set_weight("!ele_trigEff",eleTrigW[0]);
    eventInfo->set_weight("!ele_trigEff_up",eleTrigW[1]);
    eventInfo->set_weight("!ele_trigEff_down",eleTrigW[2]);

    //add veto which are not tight
    std::vector<Electron*> const& loose = event->GetPtrVec<Electron>("vetoElectrons");
    for (unsigned iEle(0); iEle<loose.size();++iEle){
      //check overlap with tight
      if (isTightElectron(loose[iEle],elecs)) continue;
      unsigned lBin = findPtEtaBin(loose[iEle]->pt(),loose[iEle]->eta(),e_ptbin_,e_etabin_);
      for (unsigned err(0); err<3;++err){
	ele_weight[err] *= eVeto_idisoDataEff_[err][lBin]/eVeto_idisoMCEff_[0][lBin];
	//ele_weight[err] *= e_gsfidSF_[err][findPtEtaBin(loose[iEle]->pt(),loose[iEle]->eta(),gsf_ptbin_,gsf_etabin_)];
	gsf_weight[err] *= e_gsfidSF_[err][findPtEtaBin(loose[iEle]->pt(),loose[iEle]->eta(),gsf_ptbin_,gsf_etabin_)];
      }
    }
    eventInfo->set_weight("!eleTight_idisoSF",ele_weight[0]);
    eventInfo->set_weight("!eleTight_idisoSF_up",ele_weight[1]);
    eventInfo->set_weight("!eleTight_idisoSF_down",ele_weight[2]);
    eventInfo->set_weight("!eleTight_gsfSF",gsf_weight[0]);
    eventInfo->set_weight("!eleTight_gsfSF_up",gsf_weight[1]);
    eventInfo->set_weight("!eleTight_gsfSF_down",gsf_weight[2]);
    tighteleweight->Fill(ele_weight[0]*gsf_weight[0]);

    //std::cout << " ele OK" << std::endl;

    std::vector<Muon*> const& mus = event->GetPtrVec<Muon>("selMuons");
    double mu_id_weight[3] = {1.0,1.0,1.0};
    double mu_iso_weight[3] = {1.0,1.0,1.0};
    double mu_tk_weight[3] = {1.0,1.0,1.0};
    for (unsigned iEle(0); iEle<mus.size();++iEle){
      unsigned lBin = findPtEtaBin(mus[iEle]->pt(),mus[iEle]->eta(),mu_ptbin_,mu_etabin_);
      unsigned mBin = findPtEtaBin(mus[iEle]->pt(),mus[iEle]->eta(),tk_ptbin_,tk_etabin_);
      for (unsigned err(0); err<3;++err){
	mu_id_weight[err] *= muTight_idSF_[err][lBin];
	mu_iso_weight[err] *= muTight_isoSF_[err][lBin];
	mu_tk_weight[err] *= mu_tkSF_[err][mBin];
      }
    }

    //add veto which are not tight
    std::vector<Muon*> const& loosemus = event->GetPtrVec<Muon>("vetoMuons");
    double muloose_0mu_weight[7] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    for (unsigned iEle(0); iEle<loosemus.size();++iEle){
      unsigned lBin = findPtEtaBin(loosemus[iEle]->pt(),loosemus[iEle]->eta(),mu_ptbin_,mu_etabin_);
      unsigned mBin = findPtEtaBin(loosemus[iEle]->pt(),loosemus[iEle]->eta(),tk_ptbin_,tk_etabin_);      
      //fill loose mu weights for veto...
      for (unsigned err(0); err<3;++err){
	muloose_0mu_weight[err] *= (1-(muVeto_idSF_[err][lBin]*muVeto_isoSF_[0][lBin]*mu_tkSF_[0][mBin]));
	//mu_veto_weight[err] *= (1-(muVeto_idDataEff_[err][lBin]*muVeto_isoDataEff_[0][lBin]))/(1-(muVeto_idMCEff_[0][lBin]*muVeto_isoMCEff_[0][lBin]));
      }
      muloose_0mu_weight[3] *= (1-(muVeto_idSF_[0][lBin]*muVeto_isoSF_[1][lBin]*mu_tkSF_[0][mBin]));
      muloose_0mu_weight[4] *= (1-(muVeto_idSF_[0][lBin]*muVeto_isoSF_[2][lBin]*mu_tkSF_[0][mBin]));
      muloose_0mu_weight[5] *= (1-(muVeto_idSF_[0][lBin]*muVeto_isoSF_[0][lBin]*mu_tkSF_[1][mBin]));
      muloose_0mu_weight[6] *= (1-(muVeto_idSF_[0][lBin]*muVeto_isoSF_[0][lBin]*mu_tkSF_[2][mBin]));
      //for (unsigned err(0); err<7;++err){
	//if (muloose_0mu_weight[err]<0) muloose_0mu_weight[err]=0;
      //}
      
      //check overlap with tight
      if (isTightMuon(loosemus[iEle],mus)) continue;
      for (unsigned err(0); err<3;++err){
	mu_id_weight[err] *= muVeto_idDataEff_[err][lBin]/muVeto_idMCEff_[0][lBin];
	mu_iso_weight[err] *= muVeto_isoDataEff_[err][lBin]/muVeto_isoMCEff_[0][lBin];
	mu_tk_weight[err] *= mu_tkSF_[err][mBin];
      }
    }
    eventInfo->set_weight("!muTight_idSF",mu_id_weight[0]);
    eventInfo->set_weight("!muTight_idSF_up",mu_id_weight[1]);
    eventInfo->set_weight("!muTight_idSF_down",mu_id_weight[2]);
    eventInfo->set_weight("!muTight_isoSF",mu_iso_weight[0]);
    eventInfo->set_weight("!muTight_isoSF_up",mu_iso_weight[1]);
    eventInfo->set_weight("!muTight_isoSF_down",mu_iso_weight[2]);
    eventInfo->set_weight("!muTight_tkSF",mu_tk_weight[0]);
    eventInfo->set_weight("!muTight_tkSF_up",mu_tk_weight[1]);
    eventInfo->set_weight("!muTight_tkSF_down",mu_tk_weight[2]);
    tightmuweight->Fill(mu_id_weight[0]*mu_iso_weight[0]*mu_tk_weight[0]);

    eventInfo->set_weight("!muLoose_0muSF",muloose_0mu_weight[0]);
    eventInfo->set_weight("!muLoose_0muSF_idup",muloose_0mu_weight[1]);
    eventInfo->set_weight("!muLoose_0muSF_iddown",muloose_0mu_weight[2]);
    eventInfo->set_weight("!muLoose_0muSF_isoup",muloose_0mu_weight[3]);
    eventInfo->set_weight("!muLoose_0muSF_isodown",muloose_0mu_weight[4]);
    eventInfo->set_weight("!muLoose_0muSF_tkup",muloose_0mu_weight[5]);
    eventInfo->set_weight("!muLoose_0muSF_tkdown",muloose_0mu_weight[6]);


    //std::cout << " mu OK" << std::endl;

    if(do_idiso_tight_weights_){
      eventInfo->set_weight("idisoTight",ele_weight[0]*gsf_weight[0]*mu_id_weight[0]*mu_iso_weight[0]*mu_tk_weight[0]);
    }
    else{
      eventInfo->set_weight("!idisoTight",ele_weight[0]*gsf_weight[0]*mu_id_weight[0]*mu_iso_weight[0]*mu_tk_weight[0]);
    }

    //std::cout << " IDISO tight done." << std::endl;

    //TO DO: id+iso veto leptons
    //first try: take leptons from W in pT,eta acceptance
    std::vector<GenParticle*> const& genParts = event->GetPtrVec<GenParticle>("genParticles");
    double ele_veto_weight[7] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    double mu_veto_weight[7] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0};

    for (unsigned iEle(0); iEle<genParts.size(); ++iEle){//Loop on genparticles

      //if (genParts[iEle]->status()!=3) continue;
      unsigned id = abs(genParts[iEle]->pdgid());
      std::vector<bool> flags=genParts[iEle]->statusFlags();
      if ( !( (flags[GenStatusBits::IsPrompt] && flags[GenStatusBits::FromHardProcess] && flags[GenStatusBits::IsFirstCopy]) ||
           (flags[GenStatusBits::IsDirectPromptTauDecayProduct]) ) ) continue;

      bool isTau = flags[GenStatusBits::IsDirectPromptTauDecayProduct];
      if (id==15) eventsWithGenTau_++;
      //do Electrons
      if (id==11){
        //std::cout << "Electron, status = " << genParts[iEle]->status() << std::endl;
        if (isTau) eventsWithGenElectronFromTau_++;
        else eventsWithGenElectron_++;
        if (genParts[iEle]->pt() > 10 && fabs(genParts[iEle]->eta()) < 2.4) {
          unsigned lBin = findPtEtaBin(genParts[iEle]->pt(),genParts[iEle]->eta(),e_ptbin_,e_etabin_);
          unsigned lBinGsf = findPtEtaBin(genParts[iEle]->pt(),genParts[iEle]->eta(),gsf_ptbin_,gsf_etabin_);
	  for (unsigned err(0); err<3;++err){
	    ele_veto_weight[err] *= (1-(eVeto_idisoDataEff_[err][lBin]*e_gsfidDataEff_[0][lBinGsf]))/(1-(eVeto_idisoMCEff_[0][lBin]*e_gsfidMCEff_[0][lBinGsf]));
	  }
	  ele_veto_weight[3] *= (1-(eVeto_idisoDataEff_[0][lBin]*e_gsfidDataEff_[1][lBinGsf]))/(1-(eVeto_idisoMCEff_[0][lBin]*e_gsfidMCEff_[0][lBinGsf]));
	  ele_veto_weight[4] *= (1-(eVeto_idisoDataEff_[0][lBin]*e_gsfidDataEff_[2][lBinGsf]))/(1-(eVeto_idisoMCEff_[0][lBin]*e_gsfidMCEff_[0][lBinGsf]));
	  double sumsqunc = sqrt(pow(e_gsfidDataEff_[0][lBinGsf]*(eVeto_idisoDataEff_[1][lBin]-eVeto_idisoDataEff_[0][lBin]),2)+pow(eVeto_idisoDataEff_[0][lBin]*(e_gsfidDataEff_[1][lBinGsf]-e_gsfidDataEff_[0][lBinGsf]),2));
	  ele_veto_weight[5] *= (1-(eVeto_idisoDataEff_[0][lBin]*e_gsfidDataEff_[0][lBinGsf]+sumsqunc))/(1-(eVeto_idisoMCEff_[0][lBin]*e_gsfidMCEff_[0][lBinGsf]));
	  sumsqunc = sqrt(pow(e_gsfidDataEff_[0][lBinGsf]*(eVeto_idisoDataEff_[2][lBin]-eVeto_idisoDataEff_[0][lBin]),2)+pow(eVeto_idisoDataEff_[0][lBin]*(e_gsfidDataEff_[2][lBinGsf]-e_gsfidDataEff_[0][lBinGsf]),2));
	  ele_veto_weight[6] *= (1-(eVeto_idisoDataEff_[0][lBin]*e_gsfidDataEff_[0][lBinGsf]-sumsqunc))/(1-(eVeto_idisoMCEff_[0][lBin]*e_gsfidMCEff_[0][lBinGsf]));
	  //for (unsigned err(0); err<7;++err){
	  //if (ele_veto_weight[err]<0) ele_veto_weight[err]=0;
	  //}

          if (isTau) eventsWithGenElectronFromTauInAcc_++;
          else eventsWithGenElectronInAcc_++;
        }
      }

      //doMuons
      if (abs(genParts[iEle]->pdgid())==13){
        //std::cout << "Muon, status = " << genParts[iEle]->status() << std::endl;
        if (isTau) eventsWithGenMuonFromTau_++;
        else eventsWithGenMuon_++;
        if (genParts[iEle]->pt() > 10 && fabs(genParts[iEle]->eta()) < 2.4) {
          unsigned lBin = findPtEtaBin(genParts[iEle]->pt(),genParts[iEle]->eta(),mu_ptbin_,mu_etabin_);
          unsigned lBinTk = findPtEtaBin(genParts[iEle]->pt(),genParts[iEle]->eta(),tk_ptbin_,tk_etabin_);
	  for (unsigned err(0); err<3;++err){
	    mu_veto_weight[err] *= (1-(muVeto_idSF_[err][lBin]*muVeto_idMCEff_[0][lBin]*muVeto_isoSF_[0][lBin]*muVeto_isoMCEff_[0][lBin]*mu_tkSF_[0][lBinTk]))/(1-(muVeto_idMCEff_[0][lBin]*muVeto_isoMCEff_[0][lBin]));
	    //mu_veto_weight[err] *= (1-(muVeto_idDataEff_[err][lBin]*muVeto_isoDataEff_[0][lBin]))/(1-(muVeto_idMCEff_[0][lBin]*muVeto_isoMCEff_[0][lBin]));
	  }
	  mu_veto_weight[3] *= (1-(muVeto_idSF_[0][lBin]*muVeto_idMCEff_[0][lBin]*muVeto_isoSF_[1][lBin]*muVeto_isoMCEff_[0][lBin]*mu_tkSF_[0][lBinTk]))/(1-(muVeto_idMCEff_[0][lBin]*muVeto_isoMCEff_[0][lBin]));
	  mu_veto_weight[4] *= (1-(muVeto_idSF_[0][lBin]*muVeto_idMCEff_[0][lBin]*muVeto_isoSF_[2][lBin]*muVeto_isoMCEff_[0][lBin]*mu_tkSF_[0][lBinTk]))/(1-(muVeto_idMCEff_[0][lBin]*muVeto_isoMCEff_[0][lBin]));
	  mu_veto_weight[5] *= (1-(muVeto_idSF_[0][lBin]*muVeto_idMCEff_[0][lBin]*muVeto_isoSF_[0][lBin]*muVeto_isoMCEff_[0][lBin]*mu_tkSF_[1][lBinTk]))/(1-(muVeto_idMCEff_[0][lBin]*muVeto_isoMCEff_[0][lBin]));
	  mu_veto_weight[6] *= (1-(muVeto_idSF_[0][lBin]*muVeto_idMCEff_[0][lBin]*muVeto_isoSF_[0][lBin]*muVeto_isoMCEff_[0][lBin]*mu_tkSF_[2][lBinTk]))/(1-(muVeto_idMCEff_[0][lBin]*muVeto_isoMCEff_[0][lBin]));

	  //for (unsigned err(0); err<7;++err){
	  //if (mu_veto_weight[err]<0) mu_veto_weight[err]=0;
	  //}
	  //double sumsqunc = sqrt(pow(muVeto_idDataEff_[0][lBin]*muVeto_isoDataEff_[0][lBin]*(mu_tkSF_[1][lBinTk]-mu_tkSF_[0][lBinTk]),2)
	  //			 +pow(muVeto_idDataEff_[0][lBin]*mu_tkSF_[0][lBinTk]*(muVeto_isoDataEff_[1][lBin]-muVeto_isoDataEff_[0][lBin]),2)+
	  //			 pow(muVeto_isoDataEff_[0][lBin]*mu_tkSF_[0][lBinTk]*(muVeto_idDataEff_[1][lBin]-muVeto_idDataEff_[0][lBin]),2));
	//mu_veto_weight[7] *= (1-(muVeto_idDataEff_[0][lBin]*muVeto_isoDataEff_[0][lBin]*mu_tkSF_[0][lBinTk]+sumsqunc))/(1-(muVeto_idMCEff_[0][lBin]*muVeto_isoMCEff_[0][lBin]));
	//sumsqunc = sqrt(pow(muVeto_idDataEff_[0][lBin]*muVeto_isoDataEff_[0][lBin]*(mu_tkSF_[2][lBinTk]-mu_tkSF_[0][lBinTk]),2)
	//		  +pow(muVeto_idDataEff_[0][lBin]*mu_tkSF_[0][lBinTk]*(muVeto_isoDataEff_[2][lBin]-muVeto_isoDataEff_[0][lBin]),2)+
	//		  pow(muVeto_isoDataEff_[0][lBin]*mu_tkSF_[0][lBinTk]*(muVeto_idDataEff_[2][lBin]-muVeto_idDataEff_[0][lBin]),2));
      //mu_veto_weight[8] *= (1-(muVeto_idDataEff_[0][lBin]*muVeto_isoDataEff_[0][lBin]*mu_tkSF_[0][lBinTk]-sumsqunc))/(1-(muVeto_idMCEff_[0][lBin]*muVeto_isoMCEff_[0][lBin]));

	  //mu_veto_weight[3] *= (1-(muVeto_idDataEff_[0][lBin]*muVeto_isoDataEff_[1][lBin]))/(1-(muVeto_idMCEff_[0][lBin]*muVeto_isoMCEff_[0][lBin]));
	  //mu_veto_weight[4] *= (1-(muVeto_idDataEff_[0][lBin]*muVeto_isoDataEff_[2][lBin]))/(1-(muVeto_idMCEff_[0][lBin]*muVeto_isoMCEff_[0][lBin]));
	  //mu_veto_weight[5] *= (1-(muVeto_idDataEff_[0][lBin]*muVeto_isoDataEff_[0][lBin]))/(1-(muVeto_idMCEff_[0][lBin]*muVeto_isoMCEff_[0][lBin]));
	  //mu_veto_weight[6] *= (1-(muVeto_idDataEff_[0][lBin]*muVeto_isoDataEff_[0][lBin]))/(1-(muVeto_idMCEff_[0][lBin]*muVeto_isoMCEff_[0][lBin]));

          //if(mu_veto_weight<0)std::cout<<"Below zero weight:"<<(1-muVeto_idisoDataEff_[lBin])/(1-muVeto_idisoMCEff_[lBin])<<" "<<muVeto_idisoDataEff_[lBin]<<" "<<muVeto_idisoMCEff_[lBin]<<std::endl;//!!
          //if(mu_veto_weight>10000)std::cout<<"Very high weight:"<<(1-muVeto_idisoDataEff_[lBin])/(1-muVeto_idisoMCEff_[lBin])<<" "<<muVeto_idisoDataEff_[lBin]<<" "<<muVeto_idisoMCEff_[lBin]<<" "<<genParts[iEle]->pt()<<" "<<genParts[iEle]->eta()<<std::endl;//!!
          if (isTau) eventsWithGenMuonFromTauInAcc_++;
          else eventsWithGenMuonInAcc_++;
        }
      }
    }//endof Loop on genparticles
    
    vetoeleweight->Fill(ele_veto_weight[0]);
    vetomuweight->Fill(mu_veto_weight[0]);
    if (do_idiso_veto_weights_) eventInfo->set_weight("idisoVeto",ele_veto_weight[0]*mu_veto_weight[0]);
    else eventInfo->set_weight("!idisoVeto",ele_veto_weight[0]*mu_veto_weight[0]);

    eventInfo->set_weight("!eleVeto_idisoSF",ele_veto_weight[0]);
    eventInfo->set_weight("!eleVeto_idisoSF_up",ele_veto_weight[1]);
    eventInfo->set_weight("!eleVeto_idisoSF_down",ele_veto_weight[2]);
    eventInfo->set_weight("!eleVeto_gsfSF_up",ele_veto_weight[3]);
    eventInfo->set_weight("!eleVeto_gsfSF_down",ele_veto_weight[4]);
    eventInfo->set_weight("!eleVeto_up",ele_veto_weight[5]);
    eventInfo->set_weight("!eleVeto_down",ele_veto_weight[6]);

    eventInfo->set_weight("!muVeto_idisotkSF",mu_veto_weight[0]);
    eventInfo->set_weight("!muVeto_idSF_up",mu_veto_weight[1]);
    eventInfo->set_weight("!muVeto_idSF_down",mu_veto_weight[2]);
    eventInfo->set_weight("!muVeto_isoSF_up",mu_veto_weight[3]);
    eventInfo->set_weight("!muVeto_isoSF_down",mu_veto_weight[4]);
    eventInfo->set_weight("!muVeto_tkSF_up",mu_veto_weight[5]);
    eventInfo->set_weight("!muVeto_tkSF_down",mu_veto_weight[6]);
    //eventInfo->set_weight("!muVeto_up",mu_veto_weight[7]);
    //eventInfo->set_weight("!muVeto_down",mu_veto_weight[8]);

    //std::cout << " IDISO veto done." << std::endl;

    }//endof Save weights


    bool zeroParton = false;

    if (do_w_soup_) {
      try{
        std::vector<GenParticle *> const& lheParticles = event->GetPtrVec<GenParticle>("lheParticles");
        double lheHT = HTFromLHEParticles(lheParticles);

      if (mc_ == mc::fall15_76X){
        //double gen_ht = eventInfo->gen_ht() ;
        double gen_ht = lheHT;
        if (100 <= gen_ht&&gen_ht <200) eventInfo->set_weight("wsoup", w1_);
        else if (200 <= gen_ht&&gen_ht <400) eventInfo->set_weight("wsoup", w2_);
        else if (400 <= gen_ht &&gen_ht<600) eventInfo->set_weight("wsoup", w3_);
        else if (gen_ht >= 600) eventInfo->set_weight("wsoup", w4_);
      }
      else if (mc_ == mc::spring16_80X || mc_ == mc::summer16_80X){
        double gen_ht = lheHT;
        if (100 <= gen_ht&&gen_ht <200) eventInfo->set_weight("wsoup", w1_);
        else if (200 <= gen_ht&&gen_ht <400) eventInfo->set_weight("wsoup", w2_);
        else if (400 <= gen_ht &&gen_ht<600) eventInfo->set_weight("wsoup", w3_);
        else if (600 <= gen_ht &&gen_ht<800) eventInfo->set_weight("wsoup", w4_);
        else if (800 <= gen_ht &&gen_ht<1200) eventInfo->set_weight("wsoup", w5_);
        else if (1200 <= gen_ht &&gen_ht<2500) eventInfo->set_weight("wsoup", w6_);
        else if (gen_ht >= 2500) eventInfo->set_weight("wsoup", w7_);
      }
      else {
        std::vector<GenParticle*> const& parts = event->GetPtrVec<GenParticle>("genParticles");
        unsigned partons = getPartonNumber(parts);//,24);

        if (partons == 1) eventInfo->set_weight("wsoup", w1_);
        if (partons == 2) eventInfo->set_weight("wsoup", w2_);
        if (partons == 3) eventInfo->set_weight("wsoup", w3_);
        if (partons >= 4) eventInfo->set_weight("wsoup", w4_);

        if (partons == 0) zeroParton = true;
      }
      } catch (...) {
        //std::cout << " -- lheParticles branch not found, go ahead ... " << std::endl;
      }
    }

    if (do_w_reweighting_ || do_dy_reweighting_ || do_z_reweighting_) { // For v_nlo_Reweighting (kfactors.root file in input/scalefactors from MIT group)
      double v_nlo_Reweight = 1.0;
      double v_pt = -50.0;
      double v_pt_oldBinning = -50.0;
      double boson_pt = -50.0;

      std::vector<GenParticle*> const& parts = event->GetPtrVec<GenParticle>("genParticles");

      TLorentzVector l1vec;
      TLorentzVector l2vec;
      bool boson_found = false;
      bool l1_found =false;
      bool l2_found =false;
      bool reco_boson_found = false;

      for (unsigned iGenPart = 0; iGenPart < parts.size(); ++iGenPart) {//Loop over gen particles
        int id = parts[iGenPart]->pdgid();
        std::vector<bool> flags=parts[iGenPart]->statusFlags();

        if ( (abs(id)==24 || abs(id)==23) && 
          parts[iGenPart]->daughters().size() > 1 && 
          abs(parts[parts[iGenPart]->daughters()[0]]->pdgid()) > 10 &&
          abs(parts[parts[iGenPart]->daughters()[0]]->pdgid()) < 17 ){// W+- || Z
            boson_pt = parts[iGenPart]->pt();
            boson_found = true;
          } else if ( (flags[GenStatusBits::IsPrompt] && parts[iGenPart]->status()==1) || 
            (flags[GenStatusBits::IsPrompt] && flags[GenStatusBits::IsDecayedLeptonHadron]) ) {
            if (id > 10 && id < 17) {
              l1vec.SetPtEtaPhiM(parts[iGenPart]->pt(), parts[iGenPart]->eta(), parts[iGenPart]->phi(), 0.);
              l1_found = true;
            }
            if (id < -10 && id > -17) {
              l2vec.SetPtEtaPhiM(parts[iGenPart]->pt(), parts[iGenPart]->eta(), parts[iGenPart]->phi(), 0.);
              l2_found = true;
            }
            if ( l1_found && l2_found ) {
              reco_boson_found = true;
            }
          }
      }//endof Loop over gen particles

      if (!boson_found && reco_boson_found) {
        TLorentzVector wzvec(l1vec);
        wzvec += l2vec;
        boson_pt = wzvec.Pt();
      }
      v_pt = boson_pt;
      v_pt_oldBinning = boson_pt;

      if (v_pt<0 || v_pt_oldBinning<0 || boson_pt<0) {
        std::cout << " SOMETHING IS GOING WRONG! " << std::endl;
        std::cout << boson_pt << std::endl;
        std::cout << v_pt << std::endl;
        std::cout << v_pt_oldBinning << std::endl;
      }

      if (v_pt<=150) {
        //std::cout << " -- Underflow! v_pt = "<< v_pt << " has been re-set to v_pt = 151.0" << std::endl;
        v_pt = 151.0;
      }
      if (v_pt>=1000) {
        //std::cout << " -- Overflow! v_pt = "<< v_pt << " has been re-set to v_pt = 1249.0" << std::endl;
        v_pt = 999.0;
      }
      if (v_pt_oldBinning<=150) {
        v_pt_oldBinning = 151.0;
      }
      if (v_pt_oldBinning>=1250) {
        v_pt_oldBinning = 1249.0;
      }

      if (do_w_reweighting_) {
        int idx0 = hist_kfactors_EWKcorr_W->FindBin(v_pt_oldBinning);
        int idx1 = hist_kfactors_WJets_012j_NLO->FindBin(v_pt_oldBinning);
        int idx2 = hist_kfactors_vbf_cnc_W->FindBin(v_pt);
        v_nlo_Reweight = (hist_kfactors_EWKcorr_W->GetBinContent(idx0)/hist_kfactors_WJets_012j_NLO->GetBinContent(idx1))*hist_kfactors_vbf_cnc_W->GetBinContent(idx2);
      } else if (do_dy_reweighting_) {
        int idx0 = hist_kfactors_EWKcorr_Z->FindBin(v_pt_oldBinning);
        int idx1 = hist_kfactors_ZJets_012j_NLO->FindBin(v_pt_oldBinning);
        int idx2 = hist_kfactors_vbf_cnc_DY->FindBin(v_pt);
        v_nlo_Reweight = (hist_kfactors_EWKcorr_Z->GetBinContent(idx0)/hist_kfactors_ZJets_012j_NLO->GetBinContent(idx1))*hist_kfactors_vbf_cnc_DY->GetBinContent(idx2);
      } else if (do_z_reweighting_) {
        int idx0 = hist_kfactors_EWKcorr_Z->FindBin(v_pt_oldBinning);
        int idx1 = hist_kfactors_ZJets_012j_NLO->FindBin(v_pt_oldBinning);
        int idx2 = hist_kfactors_vbf_cnc_Z->FindBin(v_pt);
        v_nlo_Reweight = (hist_kfactors_EWKcorr_Z->GetBinContent(idx0)/hist_kfactors_ZJets_012j_NLO->GetBinContent(idx1))*hist_kfactors_vbf_cnc_Z->GetBinContent(idx2);
      }

      eventInfo->set_weight("!v_nlo_Reweighting", v_nlo_Reweight);

    }

    if (do_ewk_w_reweighting_ || do_ewk_dy_reweighting_) { // For ewk_v_nlo_Reweighting (kFactor_V*_pT_Mjj.root file in input/scalefactors from MIT group)

      double ewk_v_nlo_Reweight = 1.0;

      double v_pt     = -50.0;
      double boson_pt = -50.0;

      //dijet mass from genjets from VBF quarks
      double vbf_digenjet_m = 0;
//       double vbf_diquark_m = 0;
      double mjj      = -50.0;

      std::vector<GenParticle*> const& parts = event->GetPtrVec<GenParticle>("genParticles");
      std::vector<GenParticle*> Quarks;

      std::vector<GenJet *> genvec;
      genvec = event->GetPtrVec<GenJet>("genJets");
      std::sort(genvec.begin(), genvec.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));

      unsigned countQuarks = 0;
      unsigned countGenjets = 0;

      TLorentzVector l1vec;
      TLorentzVector l2vec;
      bool boson_found = false;
      bool l1_found = false;
      bool l2_found = false;
      bool reco_boson_found = false;
      bool foundVBF = false;

//       std::vector<CompositeCandidate *> const& dijet_vec = event->GetPtrVec<CompositeCandidate>("jjLeadingCandidates");
//       if (dijet_vec.size() > 0) {//if dijets
//         CompositeCandidate const* dijet = dijet_vec.at(0);
//         mjj = dijet->M();
//       }

      for (unsigned iGenPart = 0; iGenPart < parts.size(); ++iGenPart) {//Loop over gen particles
        int id = parts[iGenPart]->pdgid();
        std::vector<bool> flags=parts[iGenPart]->statusFlags();

        if ( (abs(id)==24 || abs(id)==23) && 
          parts[iGenPart]->daughters().size() > 1 && 
          abs(parts[parts[iGenPart]->daughters()[0]]->pdgid()) > 10 &&
          abs(parts[parts[iGenPart]->daughters()[0]]->pdgid()) < 17 ){// W+- || Z
            boson_pt = parts[iGenPart]->pt();
            boson_found = true;
          } else if ( (flags[GenStatusBits::IsPrompt] && parts[iGenPart]->status()==1) || 
            (flags[GenStatusBits::IsPrompt] && flags[GenStatusBits::IsDecayedLeptonHadron]) ) {
            if (id > 10 && id < 17) {
              l1vec.SetPtEtaPhiM(parts[iGenPart]->pt(), parts[iGenPart]->eta(), parts[iGenPart]->phi(), 0.);
              l1_found = true;
            }
            if (id < -10 && id > -17) {
              l2vec.SetPtEtaPhiM(parts[iGenPart]->pt(), parts[iGenPart]->eta(), parts[iGenPart]->phi(), 0.);
              l2_found = true;
            }
            if ( l1_found && l2_found ) {
              reco_boson_found = true;
            }
          }

          if (abs(id)>0 && abs(id)<6 && 
              parts[iGenPart]->mothers().size()==2 && 
              parts[iGenPart]->mothers()[0]==2 && 
              parts[iGenPart]->mothers()[1]==3){
            Quarks.push_back(parts[iGenPart]);
            //std::cout << " Select " << std::endl;
          }
          //if(abs(id)>0 && abs(id)<6) parts[iGenPart]->Print();
      }//endof Loop over gen particles

      if (!boson_found && reco_boson_found) {
        TLorentzVector wzvec(l1vec);
        wzvec += l2vec;
        boson_pt = wzvec.Pt();
      }
      v_pt = boson_pt;

      if (Quarks.size()==2) foundVBF=true;
      if (!foundVBF){
        //loop again to find quarks from W or Z from WWZ and ZZZ vertices.
        for (unsigned iGenPart = 0; iGenPart < parts.size(); ++iGenPart) {//Loop over gen particles
          int id = parts[iGenPart]->pdgid();
          std::vector<bool> flags=parts[iGenPart]->statusFlags();

          if (abs(id)>0 && abs(id)<6 && 
              parts[iGenPart]->mothers().size()==1 && 
              (abs(parts[parts[iGenPart]->mothers()[0]]->pdgid())==24 || 
               abs(parts[parts[iGenPart]->mothers()[0]]->pdgid())==23)){
            Quarks.push_back(parts[iGenPart]);
          }
        }//endof Loop over gen particles
        //if (Quarks.size()==2) foundVBF=true;
      }


      if (Quarks.size()==2){
        countQuarks++;
//         vbf_diquark_m = (Quarks[0]->vector()+Quarks[1]->vector()).M();
        //get genlevel dijet mass
        unsigned genjet1 = 1000; 
        unsigned genjet2 = 1000; 
        double mindr1 = 1000;
        double mindr2 = 1000;
        for (unsigned ig(0); ig<genvec.size(); ++ig){
          double dr1 = ROOT::Math::VectorUtil::DeltaR(genvec[ig]->vector(),Quarks[0]->vector());
          if (dr1<mindr1){
            genjet1 = ig;
            mindr1 = dr1;
          }
          double dr2 = ROOT::Math::VectorUtil::DeltaR(genvec[ig]->vector(),Quarks[1]->vector());
          if (dr2<mindr2){
            genjet2 = ig;
            mindr2=dr2;
          }
        }
        //if (debug_>1) std::cout << " genjet1/2 " << genjet1 << " " << genjet2 << " mindr1=" << mindr1 << " mindr2=" << mindr2 << std::endl;
        if (genjet1<genvec.size() && genjet2<genvec.size()){
          vbf_digenjet_m = (genvec[genjet1]->vector()+genvec[genjet2]->vector()).M();
          countGenjets++;
        } else {
//           std::cout << " Warning, event genjet pair for VBF jets not found! Taking leading pair." << std::endl;
          if (genvec.size()>1) vbf_digenjet_m = (genvec[0]->vector()+genvec[1]->vector()).M();
//           std::cout << " Check mass: " << vbf_diquark_m << " " << vbf_digenjet_m << std::endl;
        }
      } else {
        std::cout << " Problem event found " << Quarks.size() << " quarks" << std::endl;
      }
      mjj = vbf_digenjet_m;

      if (v_pt<0 || boson_pt<0 || mjj<0) {
        std::cout << " Warning: SOMETHING IS GOING WRONG! " << std::endl;
        std::cout << " -- Boson pT: " << boson_pt << std::endl;
        std::cout << " -- V pT: " << v_pt << std::endl;
        std::cout << " -- Mjj: " << mjj << std::endl;
      }

      if (v_pt<200) {
        //std::cout << " -- Underflow! v_pt = "<< v_pt << " has been re-set to v_pt = 151.0" << std::endl;
        v_pt = 201.0;
      }
      if (v_pt>=1000) {
        //std::cout << " -- Overflow! v_pt = "<< v_pt << " has been re-set to v_pt = 1249.0" << std::endl;
        v_pt = 999.0;
      }
      if (mjj<200) {
        mjj = 201.0;
      }
      if (mjj>=3000) {
        mjj = 2999.0;
      }

      if (do_ewk_w_reweighting_) {
//         ewk_v_nlo_Reweight = hist_kFactors_ewk_W->GetBinContent( hist_kFactors_ewk_W->FindBin(v_pt), hist_kFactors_ewk_W->FindBin(mjj) ); WRONG!
        ewk_v_nlo_Reweight = hist_kFactors_ewk_W->GetBinContent( hist_kFactors_ewk_W->FindBin(v_pt,mjj) );
      } else if (do_ewk_dy_reweighting_) {
//         ewk_v_nlo_Reweight = hist_kFactors_ewk_Z->GetBinContent( hist_kFactors_ewk_Z->FindBin(v_pt), hist_kFactors_ewk_Z->FindBin(mjj) ); WRONG!
        ewk_v_nlo_Reweight = hist_kFactors_ewk_Z->GetBinContent( hist_kFactors_ewk_Z->FindBin(v_pt,mjj) );
      }

      eventInfo->set_weight("!ewk_v_nlo_Reweighting", ewk_v_nlo_Reweight);

    }

    if (do_prefiring_reweighting_) { // For prefiring_Reweighting (l1.root file in input/trig_factors from MIT group)
      double prefiring_Reweight            = 1.0;
      double prefiring_SingleMuon_Reweight = 1.0;
      double prefiring_JetHT_Reweight      = 1.0;
      double dijet_jet1_eta = -50.;
      double dijet_jet2_eta = -50.;
      double dijet_jet1_pt  = -50.;
      double dijet_jet2_pt  = -50.;
      double jet1_SingleMuon_weight = 0.;
      double jet1_JetHT_weight      = 0.;
      double jet1_weight            = 0.;
      double jet2_SingleMuon_weight = 0.;
      double jet2_JetHT_weight      = 0.;
      double jet2_weight            = 0.;

      std::vector<CompositeCandidate *> const& dijet_vec = event->GetPtrVec<CompositeCandidate>("jjLeadingCandidates");
      if (dijet_vec.size() > 0) {//if dijets
        CompositeCandidate const* dijet = dijet_vec.at(0);
        Candidate const* jet1 = dijet->GetCandidate("jet1");
        Candidate const* jet2 = dijet->GetCandidate("jet2");
        dijet_jet1_eta = fabs(jet1->eta());
        dijet_jet2_eta = fabs(jet2->eta());
        dijet_jet1_pt  = jet1->pt();
        dijet_jet2_pt  = jet2->pt();
      }

      if ((dijet_jet1_eta>2.25 && dijet_jet1_eta<3. && dijet_jet1_pt>40.)||
          (dijet_jet2_eta>2.25 && dijet_jet2_eta<3. && dijet_jet2_pt>40.)) {
        if ( dijet_jet1_pt<0. || dijet_jet2_pt<0. ) {
          std::cout << " Warning: SOMETHING IS GOING WRONG! " << std::endl;
          std::cout << " -- Jet1 pT: " << dijet_jet1_pt << std::endl;
          std::cout << " -- Jet2 pT: " << dijet_jet2_pt << std::endl;
        }
        if (dijet_jet1_pt>600) dijet_jet1_pt = 599.0;
        if (dijet_jet2_pt>600) dijet_jet2_pt = 599.0;

        if (dijet_jet1_eta>2.25 && dijet_jet1_eta<3. && dijet_jet1_pt>40.) {
          if (dijet_jet1_pt<=300) {
            jet1_SingleMuon_weight = hist_prefiring_etaPt_SingleMuon->GetBinContent( hist_prefiring_etaPt_SingleMuon->FindBin(dijet_jet1_pt,dijet_jet1_eta) );
            jet1_JetHT_weight      = hist_prefiring_etaPt_JetHT->GetBinContent( hist_prefiring_etaPt_JetHT->FindBin(dijet_jet1_pt,dijet_jet1_eta) );
            jet1_weight            = hist_prefiring_etaPt_SingleMuon->GetBinContent( hist_prefiring_etaPt_SingleMuon->FindBin(dijet_jet1_pt,dijet_jet1_eta) );
          } else {
            jet1_SingleMuon_weight = hist_prefiring_etaPt_SingleMuon->GetBinContent( hist_prefiring_etaPt_SingleMuon->FindBin(dijet_jet1_pt,dijet_jet1_eta) );
            if (jet1_SingleMuon_weight<0.001) {
              jet1_SingleMuon_weight = 1.;
            }
            jet1_JetHT_weight      = hist_prefiring_etaPt_JetHT->GetBinContent( hist_prefiring_etaPt_JetHT->FindBin(dijet_jet1_pt,dijet_jet1_eta) );
            jet1_weight            = hist_prefiring_etaPt_JetHT->GetBinContent( hist_prefiring_etaPt_JetHT->FindBin(dijet_jet1_pt,dijet_jet1_eta) );
          }
        }

        if (dijet_jet2_eta>2.25 && dijet_jet2_eta<3. && dijet_jet2_pt>40.) {
          if (dijet_jet2_pt<=300) {
            jet2_SingleMuon_weight = hist_prefiring_etaPt_SingleMuon->GetBinContent( hist_prefiring_etaPt_SingleMuon->FindBin(dijet_jet2_pt,dijet_jet2_eta) );
            jet2_JetHT_weight      = hist_prefiring_etaPt_JetHT->GetBinContent( hist_prefiring_etaPt_JetHT->FindBin(dijet_jet2_pt,dijet_jet2_eta) );
            jet2_weight            = hist_prefiring_etaPt_SingleMuon->GetBinContent( hist_prefiring_etaPt_SingleMuon->FindBin(dijet_jet2_pt,dijet_jet2_eta) );
          } else {
            jet2_SingleMuon_weight = hist_prefiring_etaPt_SingleMuon->GetBinContent( hist_prefiring_etaPt_SingleMuon->FindBin(dijet_jet2_pt,dijet_jet2_eta) );
            if (jet2_SingleMuon_weight<0.001) {
              jet2_SingleMuon_weight = 1.;
            }
            jet2_JetHT_weight      = hist_prefiring_etaPt_JetHT->GetBinContent( hist_prefiring_etaPt_JetHT->FindBin(dijet_jet2_pt,dijet_jet2_eta) );
            jet2_weight            = hist_prefiring_etaPt_JetHT->GetBinContent( hist_prefiring_etaPt_JetHT->FindBin(dijet_jet2_pt,dijet_jet2_eta) );
          }
        }
        prefiring_SingleMuon_Reweight = (1.-jet1_SingleMuon_weight)*(1.-jet2_SingleMuon_weight);
        prefiring_JetHT_Reweight      = (1.-jet1_JetHT_weight)*(1.-jet2_JetHT_weight);
        prefiring_Reweight            = (1.-jet1_weight)*(1.-jet2_weight);
      }
      eventInfo->set_weight("!prefiring_SingleMuon_Reweighting", prefiring_SingleMuon_Reweight);
      eventInfo->set_weight("!prefiring_JetHT_Reweighting", prefiring_JetHT_Reweight);
      eventInfo->set_weight("!prefiring_Reweighting", prefiring_Reweight); 
    }


    if (do_dy_soup_) {
      std::vector<GenParticle*> const& parts = event->GetPtrVec<GenParticle>("genParticles");
      unsigned partons = getPartonNumber(parts);//,23);

      if (partons == 1) eventInfo->set_weight("dysoup", zw1_);
      if (partons == 2) eventInfo->set_weight("dysoup", zw2_);
      if (partons == 3) eventInfo->set_weight("dysoup", zw3_);
      if (partons == 4) eventInfo->set_weight("dysoup", zw4_);
      if (partons == 0) zeroParton = true;
    }

    if (do_dy_soup_htbinned_){
      try{
        std::vector<GenParticle *> const& lheParticles = event->GetPtrVec<GenParticle>("lheParticles");
        double lheHT = HTFromLHEParticles(lheParticles);

      double gen_ht = lheHT;
      if (100 <= gen_ht&&gen_ht <200) eventInfo->set_weight("dysoup", zw1_);
      if (200 <= gen_ht&&gen_ht <400) eventInfo->set_weight("dysoup", zw2_);
      if (400 <= gen_ht &&gen_ht<600) eventInfo->set_weight("dysoup", zw3_);
      if (gen_ht >= 600) eventInfo->set_weight("dysoup", zw4_);
      } catch (...) {
        //std::cout << " -- lheParticles branch not found, go ahead ... " << std::endl;
      }
    }

    if (!save_weights_) event->Add("NoParton",zeroParton);
    //std::cout<<"Final weight: "<<eventInfo->total_weight()<<std::endl;
    return 0;

  }//execute method

  int HinvWeights::PostAnalysis() {
    if (save_weights_) {
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "-- PostAnalysis Info for HinvWeights --" << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    std::cout << " -- eventsWithGenElectron_ = " << eventsWithGenElectron_ << std::endl;
    std::cout << " -- eventsWithGenElectronInAcc_ = " << eventsWithGenElectronInAcc_ << std::endl;
    std::cout << " -- eventsWithGenElectronFromTau_ = " << eventsWithGenElectronFromTau_ << std::endl;
    std::cout << " -- eventsWithGenElectronFromTauInAcc_ = " << eventsWithGenElectronFromTauInAcc_ << std::endl;
    std::cout << " -- eventsWithGenMuon_ = " << eventsWithGenMuon_ << std::endl;
    std::cout << " -- eventsWithGenMuonInAcc_ = " << eventsWithGenMuonInAcc_ << std::endl;
    std::cout << " -- eventsWithGenMuonFromTau_ = " << eventsWithGenMuonFromTau_ << std::endl;
    std::cout << " -- eventsWithGenMuonFromTauInAcc_ = " << eventsWithGenMuonFromTauInAcc_ << std::endl;
    std::cout << " -- eventsWithGenTau_ = " << eventsWithGenTau_ << std::endl;
    }

    return 0;
  }

  void HinvWeights::PrintInfo() {
    ;
  }

  void HinvWeights::SetWTargetFractions(double f0, double f1, double f2, double f3, double f4) {
    f0_ = f0;
    f1_ = f1;
    f2_ = f2;
    f3_ = f3;
    f4_ = f4;
  }
  void HinvWeights::SetWTargetFractions(double f0, double f1, double f2, double f3, double f4, double f5, double f6, double f7) {
    f0_ = f0;
    f1_ = f1;
    f2_ = f2;
    f3_ = f3;
    f4_ = f4;
    f5_ = f5;
    f6_ = f6;
    f7_ = f7;
  }
  void HinvWeights::SetWInputYields(double n_inc, double n1, double n2, double n3, double n4) {
    n_inc_ = n_inc;
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    n4_ = n4;
  }
  void HinvWeights::SetWInputYields(double n_inc, double n1, double n2, double n3, double n4, double n5, double n6, double n7) {
    n_inc_ = n_inc;
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    n4_ = n4;
    n5_ = n5;
    n6_ = n6;
    n7_ = n7;
  }
  void HinvWeights::SetDYTargetFractions(double zf0, double zf1, double zf2, double zf3, double zf4) {
    zf0_ = zf0;
    zf1_ = zf1;
    zf2_ = zf2;
    zf3_ = zf3;
    zf4_ = zf4;
  }
  void HinvWeights::SetDYInputYields(double zn_inc, double zn1, double zn2, double zn3, double zn4) {
    zn_inc_ = zn_inc;
    zn1_ = zn1;
    zn2_ = zn2;
    zn3_ = zn3;
    zn4_ = zn4;
  }

  unsigned HinvWeights::getPartonNumber(std::vector<GenParticle*> const& parts){

    //bool count_jets = false;
    unsigned partons = 0;
    unsigned n_partons = 0;
    for (unsigned i = 0; i < parts.size(); ++i) {
      unsigned id = abs(parts[i]->pdgid());
      std::vector<bool> flags=parts[i]->statusFlags();
      if (!(flags[GenStatusBits::IsHardProcess] && flags[GenStatusBits::FromHardProcess] &&  flags[GenStatusBits::IsFirstCopy]) ) continue;
      //if (count_jets) {
      if (id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 6 || id == 21) partons++;
      //}
      //if (id == bosonid) count_jets = true;
    }
    if (partons > 4) {
      ++n_partons;
      //std::cerr << "Error making soup, event has " << partons << " partons!" << std::endl;
      //std::cout << " -- Warning making soup, event has " << partons << " partons!" << std::endl;
      //[31/10/16, 14:04:29] Anne-Marie Magnan: it's fine, it's just that at ME level I guess I was expecting at some point to have no more than 4 but pythia history tracking has changed so maybe we have now sometimes partons from PS....
    //throw;
    }
    return partons;
  }

  double HinvWeights::Efficiency(double m, double m0, double sigma, double alpha, double n, double norm){
    const double sqrtPiOver2 = 1.2533141373;
    const double sqrt2 = 1.4142135624;
    double sig = fabs((double) sigma);
    double t = (m - m0)/sig;
    if(alpha < 0)
      t = -t;
    double absAlpha = fabs(alpha/sig);
    double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    double b = absAlpha - n/absAlpha;
    double ApproxErf;
    double arg = absAlpha / sqrt2;
    if (arg > 5.) ApproxErf = 1;
    else if (arg < -5.) ApproxErf = -1;
    else ApproxErf = TMath::Erf(arg);
    double leftArea = (1 + ApproxErf) * sqrtPiOver2;
    double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
    double area = leftArea + rightArea;
    if( t <= absAlpha ){
      arg = t / sqrt2;
      if(arg > 5.) ApproxErf = 1;
      else if (arg < -5.) ApproxErf = -1;
      else ApproxErf = TMath::Erf(arg);
      return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
    }
    else{
      return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) -
        1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area;
    }
  }

  unsigned HinvWeights::findPtEtaBin(const double & aPt, const double & aEta, const std::vector<double> & ptbinvec, const std::vector<double> & etabinvec){
    if (etabinvec.size()==0 || ptbinvec.size()==0) {
      std::cout << " -- Error, uninitialised arrays for finding bin number! " << __FILE__ << " line " << __LINE__ << std::endl;
      exit(1);
    }

    unsigned nEta = etabinvec.size()-1;
    unsigned nPt = ptbinvec.size()-1;
    unsigned etabin = 0;
    unsigned ptbin = nPt-1;
    for (unsigned ieta(0); ieta<nEta;++ieta){
      if (aEta>=etabinvec[ieta] &&  aEta<etabinvec[ieta+1]) {
        etabin = ieta;
        break;
      }
    }
    for (unsigned ipt(0); ipt<nPt;++ipt){
      if (aPt>=ptbinvec[ipt] &&  aPt<ptbinvec[ipt+1]) {
        ptbin = ipt;
        break;
      }
    }


    unsigned bin = nPt*etabin+ptbin;
    //std::cout << " check bin: eta = " << aEta << " pt=" << aPt << " bineta = " << etabin << " ptbin= " << ptbin << " final bin = " << bin << std::endl;

    return bin;
  }


  void HinvWeights::fillVector(const std::string & aFileName, 
			       const unsigned nPtBins,
			       const unsigned nEtaBins,
			       std::vector<double> aVector[3],
			       std::vector<double> & ptbin,
			       std::vector<double> & etabin,
			       bool protect){
    //std::cout<<aFileName<<":"<<std::endl;//!!
    aVector[0].clear();
    aVector[1].clear();
    aVector[2].clear();
    ptbin.clear();
    etabin.clear();
    std::ifstream lInput;
    lInput.open(aFileName);
    if(!lInput.is_open()) {
      std::cerr << "Unable to open file: " << aFileName << ". Setting vector content to 1." << std::endl;
      //max expected size for e and mu is 33...
      aVector[0].resize(nPtBins*nEtaBins,1);
      aVector[1].resize(nPtBins*nEtaBins,1);
      aVector[2].resize(nPtBins*nEtaBins,1);
      return;
    }

    unsigned counter = 0;

    while(1){
      double pTmin = 0;
      double pTmax = 0;
      double etaMin = 0;
      double etaMax = 0;
      double SF = 0;
      double SFerrPlus = 0;
      double SFerrMinus = 0;
      lInput>>pTmin>>pTmax>>etaMin>>etaMax>>SF>>SFerrMinus>>SFerrPlus;
      //std::cout<<"  "<<pTmin<<" "<<pTmax<<" "<<etaMin<<" "<<etaMax<<" "<<SF<<std::endl;

      //protect against blank line at the end of the file
      if (pTmin > 1) {
	if (protect && SF>1) SF=1;
	if (SF<0) SF=0;
        aVector[0].push_back(SF);
	if (protect && (SF+SFerrPlus) > 1) aVector[1].push_back(1);
        else aVector[1].push_back(SF+SFerrPlus);
        if ((SF-SFerrMinus) < 0) aVector[2].push_back(0);
	else aVector[2].push_back(SF-SFerrMinus);

        if (counter%nPtBins==0) etabin.push_back(etaMin);
        if (counter<nPtBins) ptbin.push_back(pTmin);
        counter++;
        if (counter == nPtBins*nEtaBins) {
          etabin.push_back(etaMax);
          ptbin.push_back(pTmax);
        }
      }
      if(lInput.eof()){
        break;
      }
    }

    //std::cout << " ---- CHECK! Size of pt and eta vectors for file " << aFileName << " : ";
    //std::cout << " " << nPtBins;
    //std::cout << " " << nEtaBins << std::endl;

    lInput.close();

  }

  double HinvWeights::nloReweighting(const double & aMjj, const double & aYstar){

    double weight = 1.0;

    //double y_par0 = 8.49667e-01;
    //double y_par1 = 1.49687e-01;
    TF1 *f0c = new TF1("f0c","[0]+[1]*x "); // |Ystar|
    f0c->SetParameters(8.49667e-01, 1.49687e-01 );  //NLOmcfm/MadgraphGenJ 80M

    TF1 *f1l = new TF1("f1l","[0]+[1]*log(x)+[2]*x "); // mjj (in GeV)
    f1l->SetParameters( 3.92568e-01, 1.20734e-01, -2.55622e-04  ); //NLOmcfm/MadgraphGenJ 80M 

    weight = f0c->Eval(aYstar)*f1l->Eval(aMjj);

    return weight;
  };


}//namespace

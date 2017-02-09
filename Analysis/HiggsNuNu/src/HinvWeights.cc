#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvWeights.h"
#include "UserCode/ICHiggsTauTau/interface/TriggerPath.hh"
#include "UserCode/ICHiggsTauTau/interface/TriggerObject.hh"
#include "UserCode/ICHiggsTauTau/interface/EventInfo.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/SimpleParamParser.h"
#include "UserCode/ICHiggsTauTau/Analysis/Modules/interface/HTFromLHEParticles.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFile.h"
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
    do_1dparkedtrg_weights_           = false; //Parked analysis weights treating mjj,j2pt,l1 and hlt met separately from histogram
    do_fitted1dparkedtrg_weights_     = false; //Parked analysis weights treating mjj,j2pt,l1 and hlt met separately from fitted function
    do_3dtrg_weights_                 = false; //old parked analysis from 3D histograms from chayanit
    do_binnedin2d1dfittedtrg_weights_ = false; //Bin in two of 3 variables and fit in the other
    std::vector<std::string> thisvarorder;
    thisvarorder.push_back("Jet");
    thisvarorder.push_back("Mjj");
    thisvarorder.push_back("MET");
    binnedin2d1dfitweightvarorder_=thisvarorder;
    std::vector<double> thisvar1binning;
    thisvar1binning.push_back(30);
    thisvar1binning.push_back(40);
    thisvar1binning.push_back(50);
    thisvar1binning.push_back(60);
    thisvar1binning.push_back(1000);
    binnedin2d1dfitweightvar1binning_=thisvar1binning;
    std::vector<double> thisvar2binning;
    thisvar2binning.push_back(0);
    thisvar2binning.push_back(600);
    thisvar2binning.push_back(800);
    thisvar2binning.push_back(900);
    thisvar2binning.push_back(1000);
    thisvar2binning.push_back(5000);
    binnedin2d1dfitweightvar2binning_ = thisvar2binning;
    do_run2_                = true;
    do_metmht_              = true;
    trg_applied_in_mc_      = false;
    do_idiso_tight_weights_ = false;
    do_idiso_veto_weights_  = false;
    do_w_soup_              = false;
    do_w_reweighting_       = false;
    do_dy_soup_             = false;
    do_dy_soup_htbinned_    = false;
    do_dy_reweighting_      = false;
    do_idiso_err_           = false;
    do_idiso_errupordown_   = true;
    do_idiso_errmuore_      = true;
    do_lumixs_weights_      = false;
    save_lumixs_weights_    = true;
    input_params_ = "filelists/Apr04/ParamsApr04.dat";
    sample_name_= "test";
    input_met_ = "metNoMuons";
    input_jet_ = "pfJetsPFlow";
    trg_weight_file_="input/scale_factors/DataMCWeight_53X_v1.root";
    Alumi_  = -1;
    BClumi_ = -1;
    Dlumi_  = -1;

    errLabel.push_back("");
    errLabel.push_back("_v0Up");
    errLabel.push_back("_v0Down");
    errLabel.push_back("_v1Up");
    errLabel.push_back("_v1Down");
    errLabel.push_back("_v2Up");
    errLabel.push_back("_v2Down");

    // For v_nlo_Reweighting (kfactor_VBF_zjets_v2.root and kfactor_VBF_wjets_v2.root files in input/scalefactors from MIT group)
    kfactors_file_="input/scale_factors/kfactors.root";
    kfactor_VBF_zjets_v2_file_="input/scale_factors/kfactor_VBF_zjets_v2.root";
    kfactor_VBF_wjets_v2_file_="input/scale_factors/kfactor_VBF_wjets_v2.root";
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
    std::cout << "-- Do ID & iso weight errors ?: \t\t" << do_idiso_err_ <<std::endl;
    std::cout << "-- Do ID & iso weight error for muore ?: \t\t" << do_idiso_errmuore_ <<std::endl;
    std::cout << "-- Do ID & iso weight errors up or down?: \t\t" << do_idiso_errupordown_ <<std::endl;

    std::cout << "-- Input MET for MET HLT:  \t\t" << input_met_ << std::endl;
    std::cout << "-- Note: Input MET for MET L1 is always metNoMuons." << std::endl;
    std::cout << "-- A lumi: "<<Alumi_<<", BC lumi: "<<BClumi_<<", D lumi: "<<Dlumi_<<std::endl;

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
      if (mc_ == mc::spring16_80X){
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
    if (do_w_reweighting_ || do_dy_reweighting_) { // For v_nlo_Reweighting (kfactors.root file in input/scalefactors from MIT group)

      kfactors_ = TFile::Open(kfactors_file_.c_str());
      kfactor_VBF_zjets_v2_ = TFile::Open(kfactor_VBF_zjets_v2_file_.c_str());
      kfactor_VBF_wjets_v2_ = TFile::Open(kfactor_VBF_wjets_v2_file_.c_str());

      if (do_w_reweighting_) {
        std::cout << " -- Applying reweighting of W events to NLO from MIT (Raffaele)." << std::endl;
        hist_kfactors_N_W = (TH1F*)kfactor_VBF_wjets_v2_->Get("bosonPt_NLO_vbf");
        hist_kfactors_D_W = (TH1F*)kfactor_VBF_wjets_v2_->Get("bosonPt_LO_vbf");
//         hist_kfactors_N_W->Sumw2();
//         hist_kfactors_N_W->Scale(1./hist_kfactors_N_W->Integral());
//         hist_kfactors_D_W->Sumw2();
//         hist_kfactors_D_W->Scale(1./hist_kfactors_D_W->Integral());
        hist_kfactors_EWKcorr_W      = (TH1F*)kfactors_->Get("EWKcorr/W");
        hist_kfactors_WJets_012j_NLO = (TH1F*)kfactors_->Get("WJets_012j_NLO/nominal");
      }
      if (do_dy_reweighting_) {
        std::cout << " -- Applying reweighting of DY events to NLO from MIT (Raffaele)." << std::endl;
        hist_kfactors_N_Z = (TH1F*)kfactor_VBF_zjets_v2_->Get("bosonPt_NLO_vbf");
        hist_kfactors_D_Z = (TH1F*)kfactor_VBF_zjets_v2_->Get("bosonPt_LO_vbf");
        hist_kfactors_EWKcorr_Z      = (TH1F*)kfactors_->Get("EWKcorr/Z");
        hist_kfactors_ZJets_012j_NLO = (TH1F*)kfactors_->Get("ZJets_012j_NLO/nominal");
      }
    }



    if (save_weights_ && do_trg_weights_){
      // do weights even if not applied, to fill histo with weight for comparison !

      //get trigger scale factor histograms from file
      triggerSF_ = new TFile(trg_weight_file_.c_str());
      if(do_3dtrg_weights_){
        //OPEN 3D HISTOGRAMS
        hist_trigSF_3D.push_back((TH3F*)gDirectory->Get("h3DHLT_Dijet40_A"));
        hist_trigSF_3D.push_back((TH3F*)gDirectory->Get("h3DHLT_Dijet35_BC"));
        hist_trigSF_3D.push_back((TH3F*)gDirectory->Get("h3DHLT_Dijet30_D"));
      }
      else if(do_1dparkedtrg_weights_){
        //OPEN 1D PARKED HISTOGRAMS
        if(!do_fitted1dparkedtrg_weights_){
          hist_trigSF_METL1vec.push_back((TH1F*)gDirectory->Get("hData_MET_L1_A"));
          hist_trigSF_METL1vec.push_back((TH1F*)gDirectory->Get("hData_MET_L1_BC"));
          hist_trigSF_METL1vec.push_back((TH1F*)gDirectory->Get("hData_MET_L1_D"));

          hist_trigSF_METHLTvec.push_back((TH1F*)gDirectory->Get("hData_MET_1D_A"));
          hist_trigSF_METHLTvec.push_back((TH1F*)gDirectory->Get("hData_MET_1D_BC"));
          hist_trigSF_METHLTvec.push_back((TH1F*)gDirectory->Get("hData_MET_1D_D"));

          hist_trigSF_MjjHLTvec.push_back((TH1F*)gDirectory->Get("hData_MJJ_1D_A"));
          hist_trigSF_MjjHLTvec.push_back((TH1F*)gDirectory->Get("hData_MJJ_1D_BC"));
          hist_trigSF_MjjHLTvec.push_back((TH1F*)gDirectory->Get("hData_MJJ_1D_D"));

          hist_trigSF_JetHLTvec.push_back((TH1F*)gDirectory->Get("hData_JET2_1D_A"));
          hist_trigSF_JetHLTvec.push_back((TH1F*)gDirectory->Get("hData_JET2_1D_BC"));
          hist_trigSF_JetHLTvec.push_back((TH1F*)gDirectory->Get("hData_JET2_1D_D"));
        }
        else{
          //SET UP FITTED 1D FUNCTIONS
          func_trigSF_METL1vec.push_back((TF1*)gDirectory->Get("fData_MET_L1_A"));
          func_trigSF_METHLTvec.push_back((TF1*)gDirectory->Get("fData_MET_1D_A"));
          func_trigSF_MjjHLTvec.push_back((TF1*)gDirectory->Get("fData_MJJ_1D_A"));
          func_trigSF_JetHLTvec.push_back((TF1*)gDirectory->Get("fData_JET2_1D_A"));

          func_trigSF_METL1vec.push_back((TF1*)gDirectory->Get("fData_MET_L1_BC"));
          func_trigSF_METHLTvec.push_back((TF1*)gDirectory->Get("fData_MET_1D_BC"));
          func_trigSF_MjjHLTvec.push_back((TF1*)gDirectory->Get("fData_MJJ_1D_BC"));
          func_trigSF_JetHLTvec.push_back((TF1*)gDirectory->Get("fData_JET2_1D_BC"));

          func_trigSF_METL1vec.push_back((TF1*)gDirectory->Get("fData_MET_L1_D"));
          func_trigSF_METHLTvec.push_back((TF1*)gDirectory->Get("fData_MET_1D_D"));
          func_trigSF_MjjHLTvec.push_back((TF1*)gDirectory->Get("fData_MJJ_1D_D"));
          func_trigSF_JetHLTvec.push_back((TF1*)gDirectory->Get("fData_JET2_1D_D"));
        }
      }
      else if(do_binnedin2d1dfittedtrg_weights_){
        std::cout<<"Getting trigger efficiency functions"<<std::endl;
        for(unsigned iVar1=0;iVar1<(binnedin2d1dfitweightvar1binning_.size()-1);iVar1++){
          std::vector<std::vector<TF1*> > thisfuncvectorvector[7];
          for(unsigned iVar2=0;iVar2<(binnedin2d1dfitweightvar2binning_.size()-1);iVar2++){
            std::vector<TF1*> thisfuncvector[7];
            //HF bins
            for(unsigned iVar3=0;iVar3<(do_run2_?2:1);iVar3++){
              std::ostringstream convert;
              if (!do_metmht_) convert<<iVar1+1<<iVar2+1;
              //if (do_run2_) convert<<iVar3+1;
	      if (do_run2_) convert<<iVar2+1;
              std::string histnumber=convert.str();
              if(!do_run2_){
                thisfuncvector[0].push_back((TF1*)gDirectory->Get(("fData_"+binnedin2d1dfitweightvarorder_[2]+"_1D_"+histnumber+"A").c_str()));
                thisfuncvector[0].push_back((TF1*)gDirectory->Get(("fData_"+binnedin2d1dfitweightvarorder_[2]+"_1D_"+histnumber+"BC").c_str()));
                thisfuncvector[0].push_back((TF1*)gDirectory->Get(("fData_"+binnedin2d1dfitweightvarorder_[2]+"_1D_"+histnumber+"D").c_str()));
              }
              else{
                for (unsigned iErr(0); iErr<7;++iErr){
                  if (!do_metmht_) thisfuncvector[iErr].push_back((TF1*)gDirectory->Get(("fdata_"+binnedin2d1dfitweightvarorder_[2]+"_1d_"+histnumber+"Deff"+errLabel[iErr]).c_str()));
                  //else thisfuncvector[iErr].push_back((TF1*)gDirectory->Get(("METMHT_BIN"+histnumber+errLabel[iErr]).c_str()));
                  else thisfuncvector[iErr].push_back((TF1*)gDirectory->Get(("f1_METMHT120_MJJBIN"+histnumber).c_str()));
                }
              }
            }
            if(!do_run2_) thisfuncvectorvector[0].push_back(thisfuncvector[0]);
            for (unsigned iErr(0); iErr<7;++iErr){
              thisfuncvectorvector[iErr].push_back(thisfuncvector[iErr]);
            }
          }
          if (!do_run2_) func_trigSF_binnedin2d[0].push_back(thisfuncvectorvector[0]);
          for (unsigned iErr(0); iErr<7;++iErr){
            func_trigSF_binnedin2d[iErr].push_back(thisfuncvectorvector[iErr]);
          }
        }
        std::cout<<"-- Done!"<<std::endl;
      }
      else{
        hist_trigSF_METL1 = (TH1F*)gDirectory->Get("METL1");
        hist_trigSF_METHLT = (TH1F*)gDirectory->Get("METHLT");
        hist_trigSF_MjjHLT = (TH1F*)gDirectory->Get("MjjHLT");
        hist_trigSF_JetHLT = (TH1F*)gDirectory->Get("JetHLT");

        std::cout << " -- Content of histogram METL1 : " << std::endl;
        for (int i(0); i< hist_trigSF_METL1->GetNbinsX()+2; ++i){
          std::cout << " -- bin " << i << " [" 
              << hist_trigSF_METL1->GetXaxis()->GetBinLowEdge(i) << "-"
              << hist_trigSF_METL1->GetXaxis()->GetBinUpEdge(i) << "] : "
              <<  hist_trigSF_METL1->GetBinContent(i)
              << std::endl;
        }

        std::cout << " -- Content of histogram METHLT : " << std::endl;
        for (int i(0); i< hist_trigSF_METHLT->GetNbinsX()+2; ++i){
          std::cout << " -- bin " << i << " ["
              << hist_trigSF_METHLT->GetXaxis()->GetBinLowEdge(i) << "-"
              << hist_trigSF_METHLT->GetXaxis()->GetBinUpEdge(i) << "] : "
              << hist_trigSF_METHLT->GetBinContent(i);

          //change first bins to "1": no data/MC scale factor applied...
          if (i>0 && hist_trigSF_METHLT->GetBinContent(i) == 0) {
            hist_trigSF_METHLT->SetBinContent(i,1);
            std::cout << " -> changed to " <<  hist_trigSF_METHLT->GetBinContent(i);
          }
          std::cout	<< std::endl;

        }

        std::cout << " -- Content of histogram MjjHLT : " << std::endl;
        for (int i(0); i< hist_trigSF_MjjHLT->GetNbinsX()+2; ++i){
          std::cout << " -- bin " << i << " [" 
              << hist_trigSF_MjjHLT->GetXaxis()->GetBinLowEdge(i) << "-"
              << hist_trigSF_MjjHLT->GetXaxis()->GetBinUpEdge(i) << "] : "
              <<  hist_trigSF_MjjHLT->GetBinContent(i);
          std::cout	<< std::endl;
        }

        std::cout << " -- Content of histogram JetHLT : " << std::endl;
        for (int i(0); i< hist_trigSF_JetHLT->GetNbinsX()+2; ++i){
          std::cout << " -- bin " << i << " [" 
              << hist_trigSF_JetHLT->GetXaxis()->GetBinLowEdge(i) << "-" 
              << hist_trigSF_JetHLT->GetXaxis()->GetBinUpEdge(i) << "] : "
              <<  hist_trigSF_JetHLT->GetBinContent(i);

          //change first bins to "1": no data/MC scale factor applied...
          if (i>0 && hist_trigSF_JetHLT->GetBinContent(i) == 0) {
            hist_trigSF_JetHLT->SetBinContent(i,1);
            std::cout << " -> changed to " <<  hist_trigSF_JetHLT->GetBinContent(i);
          }
          std::cout	<< std::endl;
        }
      }
    }
    if (save_weights_){
      std::vector<double> dummypt;
      std::vector<double> dummyeta;
      fillVector("input/scale_factors/Spring16_80X_ele_tight_id_SF.txt",5,10,eTight_idisoSF_,e_ptbin_,e_etabin_);
      fillVector("input/scale_factors/Spring16_80X_gsf_id_SF.txt",1,28,e_gsfidSF_,gsf_ptbin_,gsf_etabin_);
      fillVector("input/scale_factors/Spring16_80X_mu_tight_id_SF.txt",7,8,muTight_idSF_,mu_ptbin_,mu_etabin_);
      fillVector("input/scale_factors/Spring16_80X_mu_trackingSF.txt",1,10,mu_tkSF_,tk_ptbin_,tk_etabin_);
      if(!do_idiso_err_ || (do_idiso_err_ && do_idiso_errmuore_) ){//Central value electrons
        fillVector("input/scale_factors/Spring16_80X_ele_veto_id_data_eff.txt",5,10,eVeto_idisoDataEff_,dummypt,dummyeta);
        fillVector("input/scale_factors/Spring16_80X_ele_veto_id_mc_eff.txt",5,10,eVeto_idisoMCEff_,dummypt,dummyeta);
        fillVector("input/scale_factors/Spring16_80X_gsf_id_data_eff.txt",1,28,e_gsfidDataEff_,dummypt,dummyeta);
        fillVector("input/scale_factors/Spring16_80X_gsf_id_mc_eff.txt",1,28,e_gsfidMCEff_,dummypt,dummyeta);
        //fillVector("input/scale_factors/Fall15_76X_ele_tight_id_SF.txt",eTight_idisoSF_);
        //fillVector("input/scale_factors/Fall15_76X_ele_veto_id_data_eff.txt",eVeto_idisoDataEff_);
        //fillVector("input/scale_factors/Fall15_76X_ele_veto_id_mc_eff.txt",eVeto_idisoMCEff_);
      }
      if (!do_idiso_err_ || (do_idiso_err_ && !do_idiso_errmuore_) ){//Central value muons
        fillVector("input/scale_factors/Spring16_80X_mu_tight_iso_SF.txt",7,8,muTight_isoSF_,dummypt,dummyeta);
        fillVector("input/scale_factors/Spring16_80X_mu_loose_id_data_eff.txt",7,8,muVeto_idDataEff_,dummypt,dummyeta);
        fillVector("input/scale_factors/Spring16_80X_mu_loose_iso_data_eff.txt",7,8,muVeto_isoDataEff_,dummypt,dummyeta);
        fillVector("input/scale_factors/Spring16_80X_mu_loose_id_mc_eff.txt",7,8,muVeto_idMCEff_,dummypt,dummyeta);
        fillVector("input/scale_factors/Spring16_80X_mu_loose_iso_mc_eff.txt",7,8,muVeto_isoMCEff_,dummypt,dummyeta);
        //for tracking, MC is 1, SF=dataEff
        fillVector("input/scale_factors/Spring16_80X_mu_trackingSF.txt",1,10,mu_tkDataEff_,dummypt,dummyeta);
      }
      if(do_idiso_err_ && do_idiso_errmuore_){//Muon eff varied
        fillVectorError("input/scale_factors/Spring16_80X_mu_tight_id_SF.txt",muTight_idSF_,do_idiso_errupordown_);
        fillVectorError("input/scale_factors/Spring16_80X_mu_tight_iso_SF.txt",muTight_isoSF_,do_idiso_errupordown_);
        fillVectorError("input/scale_factors/Spring16_80X_mu_loose_id_data_eff.txt",muVeto_idDataEff_,do_idiso_errupordown_);
        fillVectorError("input/scale_factors/Spring16_80X_mu_loose_iso_data_eff.txt",muVeto_isoDataEff_,do_idiso_errupordown_);
        fillVectorError("input/scale_factors/Spring16_80X_mu_loose_id_mc_eff.txt",muVeto_idMCEff_,do_idiso_errupordown_);
        fillVectorError("input/scale_factors/Spring16_80X_mu_loose_iso_mc_eff.txt",muVeto_isoMCEff_,do_idiso_errupordown_);
        fillVectorError("input/scale_factors/Spring16_80X_mu_trackingSF.txt",mu_tkSF_,do_idiso_errupordown_);
        fillVectorError("input/scale_factors/Spring16_80X_mu_trackingSF.txt",mu_tkDataEff_,do_idiso_errupordown_);
      }
      else if (do_idiso_err_ && !do_idiso_errmuore_) {//Electron eff varied
        fillVectorError("input/scale_factors/Spring16_80X_ele_tight_id_SF.txt",eTight_idisoSF_,do_idiso_errupordown_);
        fillVectorError("input/scale_factors/Spring16_80X_gsf_id_SF.txt",e_gsfidSF_,do_idiso_errupordown_);
        fillVectorError("input/scale_factors/Spring16_80X_ele_veto_id_data_eff.txt",eVeto_idisoDataEff_,do_idiso_errupordown_);
        fillVectorError("input/scale_factors/Spring16_80X_ele_veto_id_mc_eff.txt",eVeto_idisoMCEff_,do_idiso_errupordown_);
        fillVectorError("input/scale_factors/Spring16_80X_gsf_id_data_eff.txt",e_gsfidDataEff_,do_idiso_errupordown_);
        fillVectorError("input/scale_factors/Spring16_80X_gsf_id_mc_eff.txt",e_gsfidMCEff_,do_idiso_errupordown_);
        //fillVectorError("input/scale_factors/Fall15_76X_ele_tight_id_SF.txt",eTight_idisoSF_,do_idiso_errupordown_);
        //fillVectorError("input/scale_factors/Fall15_76X_ele_veto_id_data_eff.txt",eVeto_idisoDataEff_,do_idiso_errupordown_);
        //fillVectorError("input/scale_factors/Fall15_76X_ele_veto_id_mc_eff.txt",eVeto_idisoMCEff_,do_idiso_errupordown_);
      }


      for (unsigned iBin(0); iBin<muTight_idSF_.size();++iBin){
        muTight_idisoSF_.push_back(muTight_idSF_[iBin]*muTight_isoSF_[iBin]);
        if(muVeto_idDataEff_[iBin]>=1)muVeto_idDataEff_[iBin]=0.99999;
        if(muVeto_isoDataEff_[iBin]>=1)muVeto_isoDataEff_[iBin]=0.99999;
        if(muVeto_idMCEff_[iBin]>=1)muVeto_idMCEff_[iBin]=0.99999;
        if(muVeto_isoMCEff_[iBin]>=1)muVeto_isoMCEff_[iBin]=0.99999;
        if(muVeto_idDataEff_[iBin]<0)muVeto_idDataEff_[iBin]=0;
        if(muVeto_isoDataEff_[iBin]<0)muVeto_isoDataEff_[iBin]=0;
        if(muVeto_idMCEff_[iBin]<0)muVeto_idMCEff_[iBin]=0;
        if(muVeto_isoMCEff_[iBin]<0)muVeto_isoMCEff_[iBin]=0;
        muVeto_idisoDataEff_.push_back(muVeto_idDataEff_[iBin]*muVeto_isoDataEff_[iBin]);
        muVeto_idisoMCEff_.push_back(muVeto_idMCEff_[iBin]*muVeto_isoMCEff_[iBin]);
        //std::cout<<muVeto_idisoMCEff_.back()<<" "<<muVeto_idisoDataEff_.back()<<" "<<muTight_idisoSF_.back()<<std::endl;//!!
      }

      for (unsigned iBin(0); iBin<e_gsfidDataEff_.size();++iBin){
        if(mu_tkDataEff_[iBin]>=1)mu_tkDataEff_[iBin]=0.99999;
        if(mu_tkDataEff_[iBin]<0)mu_tkDataEff_[iBin]=0;
      }

      for (unsigned iBin(0); iBin<eVeto_idisoDataEff_.size();++iBin){
        if(eVeto_idisoDataEff_[iBin]>=1)eVeto_idisoDataEff_[iBin]=0.99999;
        if(eVeto_idisoMCEff_[iBin]>=1)eVeto_idisoMCEff_[iBin]=0.99999;
        if(eVeto_idisoDataEff_[iBin]<0)eVeto_idisoDataEff_[iBin]=0;
        if(eVeto_idisoMCEff_[iBin]<0)eVeto_idisoMCEff_[iBin]=0;
      }

      for (unsigned iBin(0); iBin<e_gsfidDataEff_.size();++iBin){
        if(e_gsfidDataEff_[iBin]>=1)e_gsfidDataEff_[iBin]=0.99999;
        if(e_gsfidMCEff_[iBin]>=1)e_gsfidMCEff_[iBin]=0.99999;
        if(e_gsfidDataEff_[iBin]<0)e_gsfidDataEff_[iBin]=0;
        if(e_gsfidMCEff_[iBin]<0)e_gsfidMCEff_[iBin]=0;
      }


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

      //get METnoMuons:
      Met const* metHLT = event->GetPtr<Met>(input_met_);
      Met const* metL1 = event->GetPtr<Met>("metNoMuons");

    double l1met = metL1->pt();
    double hltmet = metHLT->pt();

    double non3dtrgweight=1;
    if (do_trg_weights_){//do_trg_weights_
      if(!do_3dtrg_weights_&&!do_1dparkedtrg_weights_&&!do_binnedin2d1dfittedtrg_weights_){
        double lMax = hist_trigSF_METL1->GetXaxis()->GetBinCenter(hist_trigSF_METL1->GetNbinsX());
        double lMin = hist_trigSF_METL1->GetXaxis()->GetBinCenter(1);
        if (l1met > lMax)  l1met = lMax;
        if (l1met < lMin)  l1met = lMin;
        int lBin = hist_trigSF_METL1->GetXaxis()->FindFixBin(l1met);
        double metl1weight = hist_trigSF_METL1->GetBinContent(lBin);
        non3dtrgweight=non3dtrgweight*metl1weight;
        if (do_trg_weights_) eventInfo->set_weight("trig_metL1",metl1weight);
        else eventInfo->set_weight("!trig_metL1",metl1weight);
        //std::cout << " -- MET L1 " << l1met << " " << metl1weight << std::endl;

        lMax = hist_trigSF_METHLT->GetXaxis()->GetBinCenter(hist_trigSF_METHLT->GetNbinsX());
        lMin = hist_trigSF_METHLT->GetXaxis()->GetBinCenter(1);
        if (hltmet > lMax)  hltmet = lMax;
        if (hltmet < lMin)  hltmet = lMin;
        lBin = hist_trigSF_METHLT->GetXaxis()->FindFixBin(hltmet);
        double methltweight = hist_trigSF_METHLT->GetBinContent(lBin);
        non3dtrgweight=non3dtrgweight*methltweight;
        if (do_trg_weights_) eventInfo->set_weight("trig_metHLT",methltweight);
        else eventInfo->set_weight("!trig_metHLT",methltweight);
        //std::cout << " -- MET HLT " << hltmet << " " << methltweight << std::endl;
      }

      double mjjhltweight = 1.0;
      double jet1hltweight = 1.0;
      double jet2hltweight = 1.0;

      double mjj=0.;
      double jet1pt=0.;
      double jet2pt=0.;
      bool hasJetsInHF = false;

      //get 2 leading jets
      std::vector<CompositeCandidate *> const& dijet_vec = event->GetPtrVec<CompositeCandidate>("jjLeadingCandidates");
      if (dijet_vec.size() > 0) {//if dijets

        CompositeCandidate const* dijet = dijet_vec.at(0);

        Candidate const* jet1 = dijet->GetCandidate("jet1");
        Candidate const* jet2 = dijet->GetCandidate("jet2");

        mjj = dijet->M();
        jet1pt = jet1->pt();
        jet2pt = jet2->pt();
        hasJetsInHF = fabs(jet1->eta())>=3 || fabs(jet2->eta())>=3 ;
        //std::cout<<"mjj "<<mjj<<" j2pt "<<jet2pt<<" metl1 "<<l1met<<" hltmet "<<hltmet<<std::endl;
        unsigned nruns;
        if(!do_run2_) nruns=3;
        else nruns=2;

        if(do_3dtrg_weights_){
          //GET THE 3 3D TRIG WEIGHTS
          double trgweights[nruns];
          unsigned naxes=3;
          unsigned bins[nruns][naxes];

          //std::cout<<" Event # "<<eventInfo->event()<<std::endl;
          //std::cout<<" jet2pt: "<<jet2pt<<", hltmet: "<<hltmet<<", mjj: "<<mjj<<std::endl;
          for(unsigned irun=0;irun<nruns;irun++){
            unsigned nbins[3];
            nbins[0]=hist_trigSF_3D[irun]->GetXaxis()->GetNbins();
            nbins[1]=hist_trigSF_3D[irun]->GetYaxis()->GetNbins();
            nbins[2]=hist_trigSF_3D[irun]->GetZaxis()->GetNbins();
            bins[irun][0]=hist_trigSF_3D[irun]->GetXaxis()->FindFixBin(jet2pt);
            bins[irun][1]=hist_trigSF_3D[irun]->GetYaxis()->FindFixBin(hltmet);
            bins[irun][2]=hist_trigSF_3D[irun]->GetZaxis()->FindFixBin(mjj);
            bool axesallinrange=true;
            for(unsigned iaxis=0;iaxis<naxes;iaxis++){
              if(bins[irun][iaxis]==0){
                std::cout<<"Warning: value of axis "<<iaxis<<" for run "<<irun<<" is below trig weight histogram range setting weight to 1."<<std::endl;
                trgweights[irun]=1;
                axesallinrange=false;
              }
              else if(bins[irun][iaxis]>nbins[iaxis]){
                std::cout<<"Warning: value of axis "<<iaxis<<" for run "<<irun<<" is above trig weight histogram range setting weight to 1."<<std::endl;
                trgweights[irun]=1;
                axesallinrange=false;
              }
            }
            if(axesallinrange){
              trgweights[irun]=hist_trigSF_3D[irun]->GetBinContent(bins[irun][0],bins[irun][1],bins[irun][2]);
              //std::cout<<"Weight for run: "<<irun<<" is "<<trgweights[irun]<<std::endl;
            }
            //std::cout << " Bin Jet2pt"<<irun<<" " << bins[irun][0] << " Bin metHLT"<<irun<<" " << bins[irun][1] << " BinMjj"<<irun<<" " << bins[irun][2] << std::endl;
          }
          //std::cout<< " Weight 0 " << trgweights[0] << " Weight 1 " << trgweights[1] << " Weight 2 " << trgweights[2] << std::endl;
          //LUMI WEIGHTED AVERAGE OVER RUNS
          double trgweight3d=(trgweights[0]*Alumi_+trgweights[1]*BClumi_+trgweights[2]*Dlumi_)/(Alumi_+BClumi_+Dlumi_);
          //std::cout<<" Total Weight "<<trgweight3d<<std::endl;
          //SET TRIGGER WEIGHT
          if (do_trg_weights_) eventInfo->set_weight("trig_3d",trgweight3d);
          else eventInfo->set_weight("!trig_3d",trgweight3d);

        }//3D weights
        else if(do_1dparkedtrg_weights_){
          double trgweights[nruns];
          unsigned nvars=4;
          if(!do_fitted1dparkedtrg_weights_){
            unsigned bins[nruns][nvars];
            for(unsigned irun=0;irun<nruns;irun++){
              unsigned nbins[nvars];
              nbins[0]=hist_trigSF_METL1vec[irun]->GetXaxis()->GetNbins();
              nbins[1]=hist_trigSF_METHLTvec[irun]->GetXaxis()->GetNbins();
              nbins[2]=hist_trigSF_MjjHLTvec[irun]->GetXaxis()->GetNbins();
              nbins[3]=hist_trigSF_JetHLTvec[irun]->GetXaxis()->GetNbins();
              bins[irun][0]=hist_trigSF_METL1vec[irun]->GetXaxis()->FindFixBin(l1met);
              bins[irun][1]=hist_trigSF_METHLTvec[irun]->GetXaxis()->FindFixBin(hltmet);
              bins[irun][2]=hist_trigSF_MjjHLTvec[irun]->GetXaxis()->FindFixBin(mjj);
              bins[irun][3]=hist_trigSF_JetHLTvec[irun]->GetXaxis()->FindFixBin(jet2pt);
              bool axisinrange[4]={true,true,true,true};
              double varweights[nvars];
              for(unsigned iaxis=0;iaxis<nvars;iaxis++){
                if(bins[irun][iaxis]==0){
                  //std::cout<<"Warning: value of axis "<<iaxis<<" for run "<<irun<<" is below trig weight histogram range setting weight to 0."<<std::endl;
                  varweights[iaxis]=0;
                  axisinrange[iaxis]=false;
                }
                else if(bins[irun][iaxis]>nbins[iaxis]){
                  //std::cout<<"Warning: value of axis "<<iaxis<<" for run "<<irun<<" is above trig weight histogram range setting weight to 1."<<std::endl;
                  varweights[iaxis]=1;
                  axisinrange[iaxis]=false;
                }
              }
              if(axisinrange[0]==true)varweights[0]=hist_trigSF_METL1vec[irun]->GetBinContent(bins[irun][0]);
              if(axisinrange[1]==true)varweights[1]=hist_trigSF_METHLTvec[irun]->GetBinContent(bins[irun][1]);
              if(axisinrange[2]==true)varweights[2]=hist_trigSF_MjjHLTvec[irun]->GetBinContent(bins[irun][2]);
              if(axisinrange[3]==true)varweights[3]=hist_trigSF_JetHLTvec[irun]->GetBinContent(bins[irun][3]);
              trgweights[irun]=varweights[0]*varweights[1]*varweights[2]*varweights[3];
              //std::cout<<"Weight for run: "<<irun<<" is "<<trgweights[irun]<<std::endl;

              //std::cout << " Bin Jet2pt"<<irun<<" " << bins[irun][0] << " Bin metHLT"<<irun<<" " << bins[irun][1] << " BinMjj"<<irun<<" " << bins[irun][2] << std::endl;
            }
          }//1D parked
          else{
            //!!GET FITTED 1D WEIGHTS AND PUT IN TRGWEIGHTS[nruns]
            for(unsigned irun=0;irun<nruns;irun++){
              double xmins[4],xmaxes[4];
              double varweights[4];
              double vars[4];
              TF1* funcs[4];
              vars[0]=l1met;
              vars[1]=hltmet;
              vars[2]=mjj;
              vars[3]=jet2pt;
              funcs[0]=func_trigSF_METL1vec[irun];
              funcs[1]=func_trigSF_METHLTvec[irun];
              funcs[2]=func_trigSF_MjjHLTvec[irun];
              funcs[3]=func_trigSF_JetHLTvec[irun];
              func_trigSF_METL1vec[irun]->GetRange(xmins[0],xmaxes[0]);
              func_trigSF_METHLTvec[irun]->GetRange(xmins[1],xmaxes[1]);
              func_trigSF_MjjHLTvec[irun]->GetRange(xmins[2],xmaxes[2]);
              func_trigSF_JetHLTvec[irun]->GetRange(xmins[3],xmaxes[3]);

              //Get weight
              for(unsigned ivar=0;ivar<4;ivar++){
                if(vars[ivar]<=xmaxes[ivar]){
                  if(vars[ivar]>=xmins[ivar]){
                    varweights[ivar]=funcs[ivar]->Eval(vars[ivar]);
                  }
                  else varweights[ivar]=0;
                }
                else varweights[ivar]=1;
              }
            trgweights[irun]=varweights[0]*varweights[1]*varweights[2]*varweights[3];
            }
          }

          //LUMI WEIGHTED AVERAGE OVER RUNS
          double trgweight=(trgweights[0]*Alumi_+trgweights[1]*BClumi_+trgweights[2]*Dlumi_)/(Alumi_+BClumi_+Dlumi_);
          //std::cout<<" Total Weight "<<trgweight<<std::endl;
          //SET TRIGGER WEIGHT
          if (do_trg_weights_) eventInfo->set_weight("trig_1d",trgweight);
          else eventInfo->set_weight("!trig_1d",trgweight);
        }
        else if(do_binnedin2d1dfittedtrg_weights_){//2D-1D
          //if(l1met!=hltmet){
          //std::cout<<"Error: you must use metnomuons for both l1met and hltmet"<<std::endl;
          //return 1;
          //}
          double vars[3];
          bool found[3]={false,false,false};
          //Get the 3 variables
          for(unsigned iVar=0;iVar<binnedin2d1dfitweightvarorder_.size();iVar++){
            if(binnedin2d1dfitweightvarorder_[iVar]=="MET"){
              vars[iVar]=hltmet;
              found[0]=true;
            }
            if(binnedin2d1dfitweightvarorder_[iVar]=="Mjj"){
              vars[iVar]=mjj;
              found[1]=true;
            }
            if(binnedin2d1dfitweightvarorder_[iVar]=="Jet"){
              vars[iVar]=jet2pt;
              found[2]=true;
            }
          }
          if(!((found[0]==true)&&(found[1]==true)&&(found[2]==true))){
            std::cout<<"You must specify MET,Mjj and Jet as the variables used for 2d binned 1d trigger weights"<<std::endl;
            return 1;
          }
          //FIND WHICH BIN YOU'RE IN
          int var1bin=-10;
          int var2bin=-10;
          if(vars[0]<binnedin2d1dfitweightvar1binning_[0])var1bin=-1;
          else{
            for(unsigned iBin=0;iBin<(binnedin2d1dfitweightvar1binning_.size()-1);iBin++){
              if(vars[0]<binnedin2d1dfitweightvar1binning_[iBin+1]){
                var1bin=iBin+1;
                break;
              }
            }
            if(var1bin==-10)var1bin=binnedin2d1dfitweightvar1binning_.size()-1;
          }
          if(vars[1]<binnedin2d1dfitweightvar2binning_[0])var2bin=-1;
          else{
            for(unsigned iBin=0;iBin<(binnedin2d1dfitweightvar2binning_.size()-1);iBin++){
              if(vars[1]<binnedin2d1dfitweightvar2binning_[iBin+1]){
                var2bin=iBin+1;
                break;
              }
            }
            if(var2bin==-10)var2bin=binnedin2d1dfitweightvar2binning_.size()-1;
          }

          for (unsigned iErr(0); iErr<7;++iErr){//Loop on errors
            double trgweights[3]={0,0,0};
            double xmin,xmax;
            if((var1bin!=-1)&&(var2bin!=-1)){
              //!!READ OUT WEIGHT FOR EACH RUN
              TF1* funcs[3]={0,0,0};
              for(unsigned iRun=0;iRun<nruns;iRun++){
                funcs[iRun]=func_trigSF_binnedin2d[iErr][var1bin-1][var2bin-1][iRun];
              }

              if (!hasJetsInHF) funcs[0]->GetRange(xmin,xmax);
              else funcs[1]->GetRange(xmin,xmax);

              //Get weight
              for(unsigned iRun=0;iRun<nruns;iRun++){

                if(vars[2]<=xmax){
                  if(vars[2]>=xmin){
                    trgweights[iRun]=funcs[iRun]->Eval(vars[2]);
                  }
                  else trgweights[iRun]=0;
                }
                else trgweights[iRun]=funcs[iRun]->Eval(xmax);
              }
            }
            else{
              for(unsigned iRun=0;iRun<nruns;iRun++){
                trgweights[iRun]=0;
              }
            }
            double trgweight;
            if(!do_run2_){
              //LUMI WEIGHTED AVERAGE OVER RUNS                                                                                                      
              trgweight=(trgweights[0]*Alumi_+trgweights[1]*BClumi_+trgweights[2]*Dlumi_)/(Alumi_+BClumi_+Dlumi_);
            }
            else{
              if (!hasJetsInHF) trgweight=trgweights[0];
              else trgweight=trgweights[1];
            }
            /*if (var1bin>0&&var2bin>0) 
              std::cout<<" Total Weight "<<trgweight
              <<" vars[0]=" << vars[0] << "(" << var1bin << ")"
              <<" vars[1]=" << vars[1] << "(" << var2bin << ")"
              <<" vars[2]=" << vars[2] << "(" << xmin << "," << xmax << ")"
              <<" hasJetsInHF=" << hasJetsInHF
              <<std::endl;   */
            //SET TRIGGER WEIGHT
            eventInfo->set_weight(("!trig_2dbinned1d"+errLabel[iErr]).c_str(),trgweight);

          }//endof Loop on errors
        }//2D-1D
        else{
          double lMax = hist_trigSF_MjjHLT->GetXaxis()->GetBinCenter(hist_trigSF_MjjHLT->GetNbinsX());
          double lMin = hist_trigSF_MjjHLT->GetXaxis()->GetBinCenter(1);
          if (mjj > lMax)  mjj = lMax;
          if (mjj < lMin)  mjj = lMin;
          int lBin = hist_trigSF_MjjHLT->GetXaxis()->FindFixBin(mjj);
          mjjhltweight = hist_trigSF_MjjHLT->GetBinContent(lBin);
          non3dtrgweight=non3dtrgweight*mjjhltweight;
          if (do_trg_weights_) eventInfo->set_weight("trig_mjjHLT",mjjhltweight);
          else eventInfo->set_weight("!trig_mjjHLT",mjjhltweight);
          //std::cout << " -- Mjj HLT " << mjj << " " << mjjhltweight << std::endl;

          lMax = hist_trigSF_JetHLT->GetXaxis()->GetBinCenter(hist_trigSF_JetHLT->GetNbinsX());
          lMin = hist_trigSF_JetHLT->GetXaxis()->GetBinCenter(1);
          if (jet1pt < lMin)  jet1pt = lMin;
          if (jet1pt > lMax)  jet1pt = lMax;
          lBin = hist_trigSF_JetHLT->GetXaxis()->FindFixBin(jet1pt);
          jet1hltweight = hist_trigSF_JetHLT->GetBinContent(lBin);
          non3dtrgweight=non3dtrgweight*jet1hltweight;
          if (do_trg_weights_) eventInfo->set_weight("trig_jet1HLT",jet1hltweight);
          else eventInfo->set_weight("!trig_jet1HLT",jet1hltweight);
          //std::cout << " -- Jet1 HLT " << jet1pt << " " << jet1hltweight << std::endl;

          if (jet2pt > lMax)  jet2pt = lMax;
          if (jet2pt < lMin)  jet2pt = lMin;
          lBin = hist_trigSF_JetHLT->GetXaxis()->FindFixBin(jet2pt);
          jet2hltweight = hist_trigSF_JetHLT->GetBinContent(lBin);
          non3dtrgweight=non3dtrgweight*jet2hltweight;
          if (do_trg_weights_) eventInfo->set_weight("trig_jet2HLT",jet2hltweight);
          else eventInfo->set_weight("!trig_jet2HLT",jet2hltweight);
          //std::cout << " -- Jet2 HLT " << jet2pt << " " << jet2hlt << std::endl;
          //std::cout<<"Weight for run: 0 is "<<non3dtrgweight<<std::endl;
          //weight *= (ele_trg * tau_trg);
          //event->Add("trigweight_1", ele_trg);
          //event->Add("trigweight_2", tau_trg);
        }
      }//endof if dijets //dijet pair
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
    double ele_weight = 1.0;
    for (unsigned iEle(0); iEle<elecs.size();++iEle){
      ele_weight *= eTight_idisoSF_[findPtEtaBin(elecs[iEle]->pt(),elecs[iEle]->eta(),e_ptbin_,e_etabin_)];
      ele_weight *= e_gsfidSF_[findPtEtaBin(elecs[iEle]->pt(),elecs[iEle]->eta(),gsf_ptbin_,gsf_etabin_)];
    }
    //add veto which are not tight
    std::vector<Electron*> const& loose = event->GetPtrVec<Electron>("vetoElectrons");
    for (unsigned iEle(0); iEle<loose.size();++iEle){
      //check overlap with tight
      if (isTightElectron(loose[iEle],elecs)) continue;
      unsigned lBin = findPtEtaBin(loose[iEle]->pt(),loose[iEle]->eta(),e_ptbin_,e_etabin_);
      ele_weight *= eVeto_idisoDataEff_[lBin]/eVeto_idisoMCEff_[lBin];
      ele_weight *= e_gsfidSF_[findPtEtaBin(loose[iEle]->pt(),loose[iEle]->eta(),gsf_ptbin_,gsf_etabin_)];
    }
    eventInfo->set_weight("!eleTight_idisoSF",ele_weight);
    tighteleweight->Fill(ele_weight);

    //std::cout << " ele OK" << std::endl;

    std::vector<Muon*> const& mus = event->GetPtrVec<Muon>("selMuons");
    double mu_weight = 1.0;
    for (unsigned iEle(0); iEle<mus.size();++iEle){
      unsigned lBin = findPtEtaBin(mus[iEle]->pt(),mus[iEle]->eta(),mu_ptbin_,mu_etabin_);
      unsigned mBin = findPtEtaBin(mus[iEle]->pt(),mus[iEle]->eta(),tk_ptbin_,tk_etabin_);
      mu_weight *= muTight_idisoSF_[lBin] * mu_tkSF_[mBin];
    }
    //add veto which are not tight
    std::vector<Muon*> const& loosemus = event->GetPtrVec<Muon>("vetoMuons");
    for (unsigned iEle(0); iEle<loosemus.size();++iEle){
      //check overlap with tight
      if (isTightMuon(loosemus[iEle],mus)) continue;
      unsigned lBin = findPtEtaBin(loosemus[iEle]->pt(),loosemus[iEle]->eta(),mu_ptbin_,mu_etabin_);
      unsigned mBin = findPtEtaBin(loosemus[iEle]->pt(),loosemus[iEle]->eta(),tk_ptbin_,tk_etabin_);
      mu_weight *= muVeto_idisoDataEff_[lBin]/muVeto_idisoMCEff_[lBin] * mu_tkSF_[mBin];
    }
    eventInfo->set_weight("!muTight_idisoSF",mu_weight);
    tightmuweight->Fill(mu_weight);

    //std::cout << " mu OK" << std::endl;

    if(do_idiso_tight_weights_){
      eventInfo->set_weight("idisoTight",ele_weight*mu_weight);
    }
    else{
      eventInfo->set_weight("!idisoTight",ele_weight*mu_weight);
    }

    //std::cout << " IDISO tight done." << std::endl;

    //TO DO: id+iso veto leptons
    //first try: take leptons from W in pT,eta acceptance
    std::vector<GenParticle*> const& genParts = event->GetPtrVec<GenParticle>("genParticles");
    double ele_veto_weight = 1.0;
    double mu_veto_weight = 1.0;

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
          ele_veto_weight *= (1-(eVeto_idisoDataEff_[lBin]*e_gsfidDataEff_[lBinGsf]))/(1-(eVeto_idisoMCEff_[lBin]*e_gsfidMCEff_[lBinGsf]));
          if (isTau) eventsWithGenElectronFromTauInAcc_++;
          else eventsWithGenElectronInAcc_++;
        }
      }

      //doMuons
      if (abs(genParts[iEle]->pdgid())==13){
        //std::cout << "Muon, status = " << genParts[iEle]->status() << std::endl;
        if (isTau) eventsWithGenMuonFromTau_++;
        else eventsWithGenMuon_++;
        if (genParts[iEle]->pt() > 10 && fabs(genParts[iEle]->eta()) < 2.1) {
          unsigned lBin = findPtEtaBin(genParts[iEle]->pt(),genParts[iEle]->eta(),mu_ptbin_,mu_etabin_);
          unsigned lBinTk = findPtEtaBin(genParts[iEle]->pt(),genParts[iEle]->eta(),tk_ptbin_,tk_etabin_);
          mu_veto_weight *= (1-(muVeto_idisoDataEff_[lBin]*mu_tkDataEff_[lBinTk]))/(1-muVeto_idisoMCEff_[lBin]);
          //if(mu_veto_weight<0)std::cout<<"Below zero weight:"<<(1-muVeto_idisoDataEff_[lBin])/(1-muVeto_idisoMCEff_[lBin])<<" "<<muVeto_idisoDataEff_[lBin]<<" "<<muVeto_idisoMCEff_[lBin]<<std::endl;//!!
          //if(mu_veto_weight>10000)std::cout<<"Very high weight:"<<(1-muVeto_idisoDataEff_[lBin])/(1-muVeto_idisoMCEff_[lBin])<<" "<<muVeto_idisoDataEff_[lBin]<<" "<<muVeto_idisoMCEff_[lBin]<<" "<<genParts[iEle]->pt()<<" "<<genParts[iEle]->eta()<<std::endl;//!!
          if (isTau) eventsWithGenMuonFromTauInAcc_++;
          else eventsWithGenMuonInAcc_++;
        }
      }
    }//endof Loop on genparticles

    vetoeleweight->Fill(ele_veto_weight);
    vetomuweight->Fill(mu_veto_weight);
    if (do_idiso_veto_weights_) eventInfo->set_weight("idisoVeto",ele_veto_weight*mu_veto_weight);
    else eventInfo->set_weight("!idisoVeto",ele_veto_weight*mu_veto_weight);

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
      else if (mc_ == mc::spring16_80X){
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

    if (do_w_reweighting_ || do_dy_reweighting_) { // For v_nlo_Reweighting (kfactors.root file in input/scalefactors from MIT group)
      double v_nlo_Reweight = 1.0;
      double v_pt = 0.0;
      double v_pt_oldBinning = 0.0;

      std::vector<GenParticle*> const& parts = event->GetPtrVec<GenParticle>("genParticles");

      for (size_t idxPart = 0; idxPart < parts.size(); ++idxPart) {// Loop over genParticles

        unsigned absPdgId = TMath::Abs(parts[idxPart]->pdgid());
        std::vector<bool> flags = parts[idxPart]->statusFlags();
        if ( !(flags[GenStatusBits::IsHardProcess] && 
                flags[GenStatusBits::FromHardProcess] && 
                flags[GenStatusBits::IsFirstCopy]) ) continue;
        v_pt = parts[idxPart]->pt();
        v_pt_oldBinning = parts[idxPart]->pt();

        if (absPdgId==24) {// W+-
          if (v_pt<150) {
            //std::cout << " -- Underflow! v_pt = "<< v_pt << " has been re-set to v_pt = 151.0" << std::endl;
            v_pt = 151.0;
          }
          if (v_pt>=1000) {
            //std::cout << " -- Overflow! v_pt = "<< v_pt << " has been re-set to v_pt = 1249.0" << std::endl;
            v_pt = 999.0;
          }
          if (v_pt_oldBinning<150) {
            v_pt_oldBinning = 151.0;
          }
          if (v_pt_oldBinning>=1250) {
            v_pt_oldBinning = 1249.0;
          }
          v_nlo_Reweight = ( ( hist_kfactors_N_W->GetBinContent(hist_kfactors_N_W->FindBin(v_pt)) )/( hist_kfactors_D_W->GetBinContent(hist_kfactors_D_W->FindBin(v_pt)) ) )*( ( hist_kfactors_EWKcorr_W->GetBinContent(hist_kfactors_EWKcorr_W->FindBin(v_pt_oldBinning)) )/( hist_kfactors_WJets_012j_NLO->GetBinContent(hist_kfactors_WJets_012j_NLO->FindBin(v_pt_oldBinning)) ) );
          //std::cout << " -- The NLO weight of W is v_nlo_Reweight = "<< v_nlo_Reweight << std::endl;
        } 
        if (absPdgId==23) {// Z
          if (v_pt<150) {
            //std::cout << " -- Underflow! v_pt = "<< v_pt << " has been re-set to v_pt = 151.0" << std::endl;
            v_pt = 151.0;
          }
          if (v_pt>=1000) {
            //std::cout << " -- Overflow! v_pt = "<< v_pt << " has been re-set to v_pt = 1249.0" << std::endl;
            v_pt = 999.0;
          }
          if (v_pt_oldBinning<150) {
            v_pt_oldBinning = 151.0;
          }
          if (v_pt_oldBinning>=1250) {
            v_pt_oldBinning = 1249.0;
          }
          v_nlo_Reweight = ( ( hist_kfactors_N_Z->GetBinContent(hist_kfactors_N_Z->FindBin(v_pt)) )/( hist_kfactors_D_Z->GetBinContent(hist_kfactors_D_Z->FindBin(v_pt)) ) )*( ( hist_kfactors_EWKcorr_Z->GetBinContent(hist_kfactors_EWKcorr_Z->FindBin(v_pt_oldBinning)) )/( hist_kfactors_ZJets_012j_NLO->GetBinContent(hist_kfactors_ZJets_012j_NLO->FindBin(v_pt_oldBinning)) ) );
          //std::cout << " -- The NLO weight of Z is v_nlo_Reweight = "<< v_nlo_Reweight << std::endl;
        }
      }//endof Loop over genParticles


      eventInfo->set_weight("!v_nlo_Reweighting", v_nlo_Reweight);

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
    unsigned ptbin = 0;
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
			       std::vector<double> & aVector,
			       std::vector<double> & ptbin,
			       std::vector<double> & etabin){
    //std::cout<<aFileName<<":"<<std::endl;//!!
    aVector.clear();
    ptbin.clear();
    etabin.clear();
    std::ifstream lInput;
    lInput.open(aFileName);
    if(!lInput.is_open()) {
      std::cerr << "Unable to open file: " << aFileName << ". Setting vector content to 1." << std::endl;
      //max expected size for e and mu is 33...
      aVector.resize(nPtBins*nEtaBins,1);
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
        aVector.push_back(SF);
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

//     std::cout << " ---- Size of vector for file " << aFileName << " = " << aVector.size() << std::endl;
//     std::cout << " ---- Size of pt bin vector: " << ptbin.size() << std::endl;
//     std::cout << " ---- Size of eta bin vector: " << etabin.size() << std::endl;

    lInput.close();

  }

  void HinvWeights::fillVectorError(const std::string & aFileName, std::vector<double> & aVector, bool upordown){
    //std::cout<<aFileName<<":"<<std::endl;//!!
    std::ifstream lInput;
    lInput.open(aFileName);
    if(!lInput.is_open()){
      std::cerr << "Unable to open file: " << aFileName << ". Setting vector content to 1." << std::endl;
      //max expected size for e and mu is 33...
      aVector.resize(33,1);
      return;
    }
    while(1){
      double pTmin = 0;
      double pTmax = 0;
      double etaMin = 0;
      double etaMax = 0;
      double SF = 0;
      double SFerrPlus = 0;
      double SFerrMinus = 0;
      lInput>>pTmin>>pTmax>>etaMin>>etaMax>>SF>>SFerrMinus>>SFerrPlus;
      //protect against blank line at the end of the file
      if(upordown){
        //std::cout<<"  "<<pTmin<<" "<<pTmax<<" "<<etaMin<<" "<<etaMax<<" "<<SF+SFerrPlus<<std::endl;//!!
        if (pTmin > 1) aVector.push_back(SF+SFerrPlus);
      }
      else{
        //std::cout<<"  "<<pTmin<<" "<<pTmax<<" "<<etaMin<<" "<<etaMax<<" "<<SF-SFerrMinus<<std::endl;//!!
        if (pTmin > 1) aVector.push_back(SF-SFerrMinus);
      }

      if(lInput.eof()){
        break;
      }
    }

    std::cout << " ---- Size of vector for file " << aFileName << " = " << aVector.size() << std::endl;
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

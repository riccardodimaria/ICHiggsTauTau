#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvSampleValidationPlots.h"

#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"

#include "TMVA/Reader.h"
#include "TVector3.h"


namespace ic {

  HinvValidationPlots::HinvValidationPlots(TFileDirectory const& dir) {
    TH1F::SetDefaultSumw2();

    JetPt1_pt  = dir.make<TH1D>("JetPt1_pt", "JetPt1_pt", 250,0,500);
    JetPt2_pt  = dir.make<TH1D>("JetPt2_pt", "JetPt2_pt", 250,0,500);
    JetPt1_eta = dir.make<TH1D>("JetPt1_eta","JetPt1_eta",100,-5.0,5.0);
    JetPt2_eta = dir.make<TH1D>("JetPt2_eta","JetPt2_eta",100,-5.0,5.0);
    JetPt1_phi = dir.make<TH1D>("JetPt1_phi","JetPt1_phi", 64,-3.2,3.2);    
    JetPt2_phi = dir.make<TH1D>("JetPt2_phi","JetPt2_phi", 64,-3.2,3.2);
    
  };


  HinvSampleValidationPlots::HinvSampleValidationPlots(std::string const& name): ModuleBase(name){
    fs_ = NULL;
    met_label_ = "pfMet";
    dijet_label_ = "jjCandidates";
    genparticles_label_ = "genParticlesTaus";
    sel_label_ = "JetPair";
    channel_ = "nunu";
    is_data_ = false;
    is_embedded_ = false;
  }

  HinvSampleValidationPlots::~HinvSampleValidationPlots(){;}


  int  HinvSampleValidationPlots::PreAnalysis(){
    std::cout << "-------------------------------------------- " << std::endl
	      << "** PreAnalysis Info for HinvSampleValidationPlots **" << std::endl
	      << "-------------------------------------------- " << std::endl;
    if (fs_) {
      std::cout << "MET Label: " << met_label_ << std::endl;
      std::cout << "dijet Label: " << dijet_label_ << std::endl;
      std::cout << "gen particles label: " << genparticles_label_ << std::endl;
      std::cout << "Selection Label: " << sel_label_ << std::endl;
      std::cout << "Channel :" << channel_ << std::endl;
      if (is_embedded_ && is_data_) std::cout << "Processing set for embedded data !" << std::endl;
      else if (is_data_) std::cout << "Processing set for data !" << std::endl;
      else if (is_embedded_) std::cout << "Processing set for embedded MC !" << std::endl;
      else  std::cout << "Processing set for MC !" << std::endl;
    }


    InitValidationPlots();

    yields_ = 0;

    return 0;
  }

  int  HinvSampleValidationPlots::Execute(TreeEvent *event){

    // Get the objects we need from the event
    EventInfo const* eventInfo    = event->GetPtr<EventInfo>("eventInfo");
    EventInfo * eventInfoNonConst = event->GetPtr<EventInfo>("eventInfo");

    wt_ = eventInfo->total_weight();
    
    std::vector<CompositeCandidate *> const& dijet_vec = event->GetPtrVec<CompositeCandidate>(dijet_label_);

    Met const* met           = event->GetPtr<Met>(met_label_);
    Met const* met_noMuons   = event->GetPtr<Met>("metNoMuons");
    std::vector<PFJet*> jets = event->GetPtrVec<PFJet>("pfJetsPFlow");
    std::sort(jets.begin(), jets.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
    
    std::vector<PFJet*> alljets = event->GetPtrVec<PFJet>("AllpfJetsPFlow");
    std::sort(alljets.begin(), alljets.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
    
    
    if (dijet_vec.size() != 0) {
      
      CompositeCandidate const* dijet = dijet_vec.at(0);

      Candidate const* jet1 = dijet->GetCandidate("jet1");
      Candidate const* jet2 = dijet->GetCandidate("jet2");

      double jpt_1_  = jet1->pt();
      double jeta_1_ = jet1->eta();
      double jpt_2_  = jet2->pt();
      double jeta_2_ = jet2->eta();
      
      m_valPlots->JetPt1_pt ->Fill(jet1->pt(), wt_);
      m_valPlots->JetPt2_pt ->Fill(jet2->pt(), wt_);
      m_valPlots->JetPt1_eta->Fill(jet1->eta(),wt_);
      m_valPlots->JetPt2_eta->Fill(jet2->eta(),wt_);
      m_valPlots->JetPt1_phi->Fill(jet1->phi(),wt_);
      m_valPlots->JetPt2_phi->Fill(jet2->phi(),wt_);
 
    }
    else if (jets.size() > 1) {
      
    }

    return 0;
  }

  int  HinvSampleValidationPlots::PostAnalysis(){
    return 0;
  }

  void  HinvSampleValidationPlots::PrintInfo(){
    ;
  }

  void HinvSampleValidationPlots::InitValidationPlots() {
    m_valPlots = new HinvValidationPlots(fs_->mkdir(sel_label_));
    std::cout << " m_valPlots initialised" << std::endl;
  }
      
  void HinvSampleValidationPlots::FillYields() {
    yields_ = yields_ + wt_;
  }
}//namespace



#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvQCDVBFStudyGenLevel.h"

#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "UserCode/ICHiggsTauTau/interface/GenParticle.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"

#include "TMVA/Reader.h"
#include "TVector3.h"

using namespace std;

namespace ic {

  HinvQCDVBFStudyGenLevelMETPlots::HinvQCDVBFStudyGenLevelMETPlots(TFileDirectory const& dir) {
    TH1D::SetDefaultSumw2();
    TH2D::SetDefaultSumw2();
    
    cout << "TDir Dump:" << endl;
    cout << dir.fullPath() << endl;
    
    gen_sumNu_Et  = dir.make<TH1D>("gen_sumNu_Et", "gen_sumNu_Et", 250,0,500);
    reco_met      = dir.make<TH1D>("reco_met",     "reco_met",     250,0,500);
    genVsReco_met = dir.make<TH2D>("genVsReco_met","genVsReco_met",250,0,500,250,0,500);
  };


  HinvQCDVBFStudyGenLevel::HinvQCDVBFStudyGenLevel(std::string const& name): ModuleBase(name){
    fs_ = NULL;
    met_label_ = "pfMet";
    dijet_label_ = "jjCandidates";
    genparticles_label_ = "genParticlesTaus";
    sel_label_ = "JetPair";
    channel_ = "nunu";
    is_data_ = false;
    is_embedded_ = false;
  }

  HinvQCDVBFStudyGenLevel::~HinvQCDVBFStudyGenLevel(){;}


  int  HinvQCDVBFStudyGenLevel::PreAnalysis(){
    std::cout << "-------------------------------------------- " << std::endl
	      << "** PreAnalysis Info for HinvQCDVBFStudyGenLevel **" << std::endl
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

  int  HinvQCDVBFStudyGenLevel::Execute(TreeEvent *event){
// cout << "________________" << endl;
    // Get the objects we need from the event
    EventInfo const             *eventInfo         = event->GetPtr<EventInfo>("eventInfo");
    EventInfo                   *eventInfoNonConst = event->GetPtr<EventInfo>("eventInfo");
    Met const                   *met               = event->GetPtr<Met>(met_label_);
    Met const                   *genMet            = event->GetPtr<Met>("pfMetGen");
    vector<GenParticle*> const&  genPart           = event->GetPtrVec<GenParticle>("genParticles");
     
//     ROOT::Math::PtEtaPhiEVector vecSumNu(0,0,0,0);
// 
//     cout << "GenParticles size = " << genPart.size() << endl;
//     for (unsigned i=0; i < genPart.size(); i++) {
//    
//       cout << "p" << i << " st=" << genPart[i]->status() << " id=" << genPart[i]->pdgid()<< " pt=" << genPart[i]->vector() << endl;      
//       if (genPart[i]->status() == 1){
// 
// 	
//         int id = genPart[i]->pdgid();
//         if(abs(id)==12 || abs(id)==14 || abs(id)==16){
// 	  cout << "Found a neutrino = " << id << " pt=" << genPart[i]->vector() << endl;
// 	  vecSumNu = vecSumNu+genPart[i]->vector();
// 	}
//       }
//     }

//     cout << "Vec Sum (nu): " << vecSumNu.pt() << endl;
//     cout << "Gen met     : " << genMet->pt() << endl;
//     cout << "Reco met    : " << met->pt() << endl;
  
    m_genLevelMETPlots->gen_sumNu_Et ->Fill(genMet->pt());
    m_genLevelMETPlots->reco_met     ->Fill(met->pt());
    m_genLevelMETPlots->genVsReco_met->Fill(genMet->pt(),met->pt());

    return 0;
  }

  int  HinvQCDVBFStudyGenLevel::PostAnalysis(){
    return 0;
  }

  void  HinvQCDVBFStudyGenLevel::PrintInfo(){
    ;
  }

  void HinvQCDVBFStudyGenLevel::InitValidationPlots() {
    m_genLevelMETPlots = new HinvQCDVBFStudyGenLevelMETPlots(fs_->mkdir(sel_label_));
    std::cout << " m_valPlots initialised" << std::endl;
  }
      
  void HinvQCDVBFStudyGenLevel::FillYields() {
    yields_ = yields_ + wt_;
  }
}//namespace



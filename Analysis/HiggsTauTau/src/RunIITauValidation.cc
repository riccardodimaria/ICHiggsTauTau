#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/RunIITauValidation.h"
#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "UserCode/ICHiggsTauTau/interface/Tau.hh"
//#include "DataFormats/TauReco/interface/PFTauFwd.h"
//#include "DataFormats/TauReco/interface/PFTau.h"
//#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"

#include "TVector3.h"
#include "boost/format.hpp"
#include "TMath.h"
#include<string>

namespace ic {


	RunIITauValidation::RunIITauValidation(std::string const& name) : ModuleBase(name) { 
		//add a few things configurable from the .cpp script
		write_plots_ = false;
		write_tree_ = true;
		fs_ = NULL;
	}


	RunIITauValidation::~RunIITauValidation() {
		;
	}

int RunIITauValidation::PreAnalysis() {
	std::cout << "-------------------------------------" << std::endl;
	std::cout << "RunIITauValidation" << std::endl;
	std::cout << "-------------------------------------" << std::endl;    
	if (fs_) {
		std::cout << boost::format(param_fmt()) % "write_plots"     % write_plots_;
		std::cout << boost::format(param_fmt()) % "write_tree"      % write_tree_;


		//Setup a flat tree for the outputs, if this is what you want to use
  if (write_tree_) {
   	outtree_ = fs_->make<TTree>("tauntuple","tauntuple");
        //Electrons:
         outtree_->Branch("taupt", &taupt_);
         outtree_->Branch("taueta", &taueta_);
         outtree_->Branch("tauphi", &tauphi_);
         outtree_->Branch("taudm", &taudm_);
         outtree_->Branch("taudmf",&taudmf_);
         outtree_->Branch("taudmfnew",&taudmfnew_);
         outtree_->Branch("taucombinediso", &taucombined_);
         				//etc etc for whatever is important
			}
			if(write_plots_) {
				//Instead could initiate plots here which can be saved to the same file_
	}

		}

		return 0;
	}

	int RunIITauValidation::Execute(TreeEvent *event) {

		// Get the objects we need from the event
		EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
                std::vector<Tau*> taus = event->GetPtrVec<Tau>("taus");
               // std::vector<Muon*> muons = event->GetPtrVec<Muon>("muons");
		//std::vector<Tau*> taus = event->GetPtrVec<Tau>("taus");
		//Met const*  met = event->GetPtr<Met>("pfMet");
	//	std::vector<GenParticle*> genparticles;
	//	genparticles = event->GetPtrVec<GenParticle>("genParticles");
	//	std::vector<GenJet*> genjets = event->GetPtrVec<GenJet>("genJets");
	//	std::vector<Track*> tracks = event->GetPtrVec<Track>("tracks");
	//	std::vector<Vertex*> vertex = event->GetPtrVec<Vertex>("vertices");



		if(write_tree_){

                  //taus:
                  for(UInt_t ielec=0;ielec<taus.size();++ielec){
                   taupt_ = taus.at(ielec)->pt();
                   taueta_ = taus.at(ielec)->eta();
                   tauphi_ = taus.at(ielec)->phi();
                   taudmf_ = taus.at(ielec)->GetTauID("decayModeFinding");
                   taudmfnew_ = taus.at(ielec)->GetTauID("decayModeFindingNewDMs");
                   if(taudmf_>0.5 || taudmfnew_>0.5){
                   taudm_ = taus.at(ielec)->decay_mode();
                   } else taudm_ = -999;
                   taucombined_ = taus.at(ielec)->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
           
                           outtree_->Fill();
                   }

		}
		return 0;
	}

	int RunIITauValidation::PostAnalysis() {
		std::cout << "-------------------------------------" << std::endl;
		std::cout << "RunIITauValidation" << std::endl;
		std::cout << "-------------------------------------" << std::endl;      
		return 0;
	}
	void RunIITauValidation::PrintInfo() {
		;
	}
}


#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/RunIIMuonValidation.h"
#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "UserCode/ICHiggsTauTau/interface/Electron.hh"
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


	RunIIMuonValidation::RunIIMuonValidation(std::string const& name) : ModuleBase(name) { 
		//add a few things configurable from the .cpp script
		write_plots_ = false;
		write_tree_ = true;
		fs_ = NULL;
	}


	RunIIMuonValidation::~RunIIMuonValidation() {
		;
	}

int RunIIMuonValidation::PreAnalysis() {
	std::cout << "-------------------------------------" << std::endl;
	std::cout << "RunIIMuonValidation" << std::endl;
	std::cout << "-------------------------------------" << std::endl;    
	if (fs_) {
		std::cout << boost::format(param_fmt()) % "write_plots"     % write_plots_;
		std::cout << boost::format(param_fmt()) % "write_tree"      % write_tree_;


		//Setup a flat tree for the outputs, if this is what you want to use
  if (write_tree_) {
   	outtree_ = fs_->make<TTree>("muonntuple","muonntuple");
        //Electrons:
         outtree_->Branch("muonpt", &muonpt_);
         outtree_->Branch("muoneta", &muoneta_);
         outtree_->Branch("muonphi", &muonphi_);
         outtree_->Branch("muonchargediso",&muonchargediso_);
         outtree_->Branch("muonchargedisoall",&muonchargedisoall_);
         outtree_->Branch("muonneutraliso",&muonneutraliso_);
         outtree_->Branch("muongammaiso",&muongammaiso_);
         outtree_->Branch("muonpuiso",&muonpuiso_);
				//etc etc for whatever is important
			}
			if(write_plots_) {
				//Instead could initiate plots here which can be saved to the same file_
	}

		}

		return 0;
	}

	int RunIIMuonValidation::Execute(TreeEvent *event) {

		// Get the objects we need from the event
		EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
             //   std::vector<Electron*> electrons = event->GetPtrVec<Electron>("electrons");
                std::vector<Muon*> muons = event->GetPtrVec<Muon>("muons");
		//std::vector<Tau*> taus = event->GetPtrVec<Tau>("taus");
		//Met const*  met = event->GetPtr<Met>("pfMet");
	//	std::vector<GenParticle*> genparticles;
	//	genparticles = event->GetPtrVec<GenParticle>("genParticles");
	//	std::vector<GenJet*> genjets = event->GetPtrVec<GenJet>("genJets");
	//	std::vector<Track*> tracks = event->GetPtrVec<Track>("tracks");
	//	std::vector<Vertex*> vertex = event->GetPtrVec<Vertex>("vertices");



		if(write_tree_){

                  //muons:
                  for(UInt_t ielec=0;ielec<muons.size();++ielec){
                   muonpt_ = muons.at(ielec)->pt();
                   muoneta_ = muons.at(ielec)->eta();
                   muonphi_ = muons.at(ielec)->phi();
                   muonchargediso_ = muons.at(ielec)->dr03_pfiso_charged();
                   muonchargedisoall_ = muons.at(ielec)->dr03_pfiso_charged_all();
                   muonneutraliso_ = muons.at(ielec)->dr03_pfiso_neutral();
                   muongammaiso_ = muons.at(ielec)->dr03_pfiso_gamma();
                   muonpuiso_ = muons.at(ielec)->dr03_pfiso_pu();
                   outtree_->Fill();
                   }

		}
		return 0;
	}

	int RunIIMuonValidation::PostAnalysis() {
		std::cout << "-------------------------------------" << std::endl;
		std::cout << "RunIIMuonValidation" << std::endl;
		std::cout << "-------------------------------------" << std::endl;      
		return 0;
	}
	void RunIIMuonValidation::PrintInfo() {
		;
	}
}


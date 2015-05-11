#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/RunIIElectronValidation.h"
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


	RunIIElectronValidation::RunIIElectronValidation(std::string const& name) : ModuleBase(name) { 
		//add a few things configurable from the .cpp script
		write_plots_ = false;
		write_tree_ = true;
		fs_ = NULL;
	}


	RunIIElectronValidation::~RunIIElectronValidation() {
		;
	}

int RunIIElectronValidation::PreAnalysis() {
	std::cout << "-------------------------------------" << std::endl;
	std::cout << "RunIIElectronValidation" << std::endl;
	std::cout << "-------------------------------------" << std::endl;    
	if (fs_) {
		std::cout << boost::format(param_fmt()) % "write_plots"     % write_plots_;
		std::cout << boost::format(param_fmt()) % "write_tree"      % write_tree_;


		//Setup a flat tree for the outputs, if this is what you want to use
  if (write_tree_) {
   	outtree_ = fs_->make<TTree>("electronntuple","electronntuple");
        //Electrons:
         outtree_->Branch("electronpt", &electronpt_);
         outtree_->Branch("electroneta", &electroneta_);
         outtree_->Branch("electronphi", &electronphi_);
         outtree_->Branch("electronidmva",&electronidmva_);
         outtree_->Branch("electronconversion",&electronconversion_);
         outtree_->Branch("electronchargediso",&electronchargediso_);
         outtree_->Branch("electronchargedisoall",&electronchargedisoall_);
         outtree_->Branch("electronneutraliso",&electronneutraliso_);
         outtree_->Branch("electrongammaiso",&electrongammaiso_);
         outtree_->Branch("electronpuiso",&electronpuiso_);
				//etc etc for whatever is important
			}
			if(write_plots_) {
				//Instead could initiate plots here which can be saved to the same file_
	}

		}

		return 0;
	}

	int RunIIElectronValidation::Execute(TreeEvent *event) {

		// Get the objects we need from the event
		EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
                std::vector<Electron*> electrons = event->GetPtrVec<Electron>("electrons");
               // std::vector<Muon*> muons = event->GetPtrVec<Muon>("muons");
		//std::vector<Tau*> taus = event->GetPtrVec<Tau>("taus");
		//Met const*  met = event->GetPtr<Met>("pfMet");
	//	std::vector<GenParticle*> genparticles;
	//	genparticles = event->GetPtrVec<GenParticle>("genParticles");
	//	std::vector<GenJet*> genjets = event->GetPtrVec<GenJet>("genJets");
	//	std::vector<Track*> tracks = event->GetPtrVec<Track>("tracks");
	//	std::vector<Vertex*> vertex = event->GetPtrVec<Vertex>("vertices");



		if(write_tree_){

                  //electrons:
                  for(UInt_t ielec=0;ielec<electrons.size();++ielec){
                   electronpt_ = electrons.at(ielec)->pt();
                   electroneta_ = electrons.at(ielec)->eta();
                   electronphi_ = electrons.at(ielec)->phi();
                   electronidmva_ = electrons.at(ielec)->GetIdIso("mvaNonTrigV025nsPHYS14");
                   electronconversion_ = (int)electrons.at(ielec)->has_matched_conversion();
                   electronchargediso_ = electrons.at(ielec)->dr03_pfiso_charged();
                   electronchargedisoall_ = electrons.at(ielec)->dr03_pfiso_charged_all();
                   electronneutraliso_ = electrons.at(ielec)->dr03_pfiso_neutral();
                   electrongammaiso_ = electrons.at(ielec)->dr03_pfiso_gamma();
                   electronpuiso_ = electrons.at(ielec)->dr03_pfiso_pu();
                   outtree_->Fill();
                   }

		}
		return 0;
	}

	int RunIIElectronValidation::PostAnalysis() {
		std::cout << "-------------------------------------" << std::endl;
		std::cout << "RunIIElectronValidation" << std::endl;
		std::cout << "-------------------------------------" << std::endl;      
		return 0;
	}
	void RunIIElectronValidation::PrintInfo() {
		;
	}
}


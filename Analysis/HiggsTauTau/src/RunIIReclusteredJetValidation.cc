#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/RunIIReclusteredJetValidation.h"
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


	RunIIReclusteredJetValidation::RunIIReclusteredJetValidation(std::string const& name) : ModuleBase(name) { 
		//add a few things configurable from the .cpp script
		write_plots_ = false;
		write_tree_ = true;
                from_miniaod_=false;
		fs_ = NULL;
	}


	RunIIReclusteredJetValidation::~RunIIReclusteredJetValidation() {
		;
	}

int RunIIReclusteredJetValidation::PreAnalysis() {
	std::cout << "-------------------------------------" << std::endl;
	std::cout << "RunIIReclusteredJetValidation" << std::endl;
	std::cout << "-------------------------------------" << std::endl;    
	if (fs_) {
		std::cout << boost::format(param_fmt()) % "write_plots"     % write_plots_;
		std::cout << boost::format(param_fmt()) % "write_tree"      % write_tree_;


		//Setup a flat tree for the outputs, if this is what you want to use
  if (write_tree_) {
   	outtree_ = fs_->make<TTree>("reclusteredjetntuple","reclusteredjetntuple");
        //Electrons:
         outtree_->Branch("jetpt", &jetpt_);
         outtree_->Branch("jeteta", &jeteta_);
         outtree_->Branch("jetphi", &jetphi_);
         outtree_->Branch("jetpdg", &jetpdg_);
         outtree_->Branch("pujetidmva",&pujetidmva_);
         outtree_->Branch("btagsimplehigh",&btagsimplehigh_);
         outtree_->Branch("btagsimplepur",&btagsimplepur_);
         outtree_->Branch("btagcombined", &btagcombined_);

                  				//etc etc for whatever is important
			}
			if(write_plots_) {
				//Instead could initiate plots here which can be saved to the same file_
	}

		}

		return 0;
	}

	int RunIIReclusteredJetValidation::Execute(TreeEvent *event) {

		// Get the objects we need from the event
		EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
		std::vector<PFJet*> jets = event->GetPtrVec<PFJet>("ak4PFJetsCHSReclustered");
//                std::vector<Tau*> taus = event->GetPtrVec<Electron>("taus");
               // std::vector<Muon*> muons = event->GetPtrVec<Muon>("muons");
		//std::vector<Tau*> taus = event->GetPtrVec<Tau>("taus");
		//Met const*  met = event->GetPtr<Met>("pfMet");
	//	std::vector<GenParticle*> genparticles;
	//	genparticles = event->GetPtrVec<GenParticle>("genParticles");
	//	std::vector<GenJet*> genjets = event->GetPtrVec<GenJet>("genJets");
	//	std::vector<Track*> tracks = event->GetPtrVec<Track>("tracks");
	//	std::vector<Vertex*> vertex = event->GetPtrVec<Vertex>("vertices");



		if(write_tree_){
   for(UInt_t ielec=0;ielec<jets.size();++ielec){
                   jetpt_ = jets.at(ielec)->pt();
                   jeteta_ = jets.at(ielec)->eta();
                   jetphi_ = jets.at(ielec)->phi();
                   pujetidmva_ = jets.at(ielec)->pu_id_mva_value();
                   jetpdg_ = jets.at(ielec)->parton_flavour(); 
                   btagsimplehigh_ = jets.at(ielec)->GetBDiscriminator("simpleSecondaryVertexHighEffBJetTags");
                   btagsimplepur_ = jets.at(ielec)->GetBDiscriminator("simpleSecondaryVertexHighPurBJetTags");
                   btagcombined_ = jets.at(ielec)->GetBDiscriminator("combinedSecondaryVertexBJetTags");
               outtree_->Fill();
             }
         }

		return 0;
	}

	int RunIIReclusteredJetValidation::PostAnalysis() {
		std::cout << "-------------------------------------" << std::endl;
		std::cout << "RunIIReclusteredJetValidation" << std::endl;
		std::cout << "-------------------------------------" << std::endl;      
		return 0;
	}
	void RunIIReclusteredJetValidation::PrintInfo() {
		;
	}
}


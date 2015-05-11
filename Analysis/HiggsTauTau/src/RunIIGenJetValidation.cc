#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/RunIIGenJetValidation.h"
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


	RunIIGenJetValidation::RunIIGenJetValidation(std::string const& name) : ModuleBase(name) { 
		//add a few things configurable from the .cpp script
		write_plots_ = false;
		write_tree_ = true;
                from_miniaod_=false;
		fs_ = NULL;
	}


	RunIIGenJetValidation::~RunIIGenJetValidation() {
		;
	}

int RunIIGenJetValidation::PreAnalysis() {
	std::cout << "-------------------------------------" << std::endl;
	std::cout << "RunIIGenJetValidation" << std::endl;
	std::cout << "-------------------------------------" << std::endl;    
	if (fs_) {
		std::cout << boost::format(param_fmt()) % "write_plots"     % write_plots_;
		std::cout << boost::format(param_fmt()) % "write_tree"      % write_tree_;


		//Setup a flat tree for the outputs, if this is what you want to use
  if (write_tree_) {
   	outtree_ = fs_->make<TTree>("genjetntuple","genjetntuple");
        //Electrons:
         outtree_->Branch("genjetpt", &genjetpt_);
         outtree_->Branch("genjeteta", &genjeteta_);
         outtree_->Branch("genjetphi", &genjetphi_);
         outtree_->Branch("genjetconstituents_", &genjetconstituents_);

         outtree_->Branch("slimmedgenjetpt", &slimmedgenjetpt_);
         outtree_->Branch("slimmedgenjeteta", &slimmedgenjeteta_);
         outtree_->Branch("slimmedgenjetphi", &slimmedgenjetphi_);
         outtree_->Branch("slimmedgenjetconstituents", &slimmedgenjetconstituents_);

         				//etc etc for whatever is important
			}
			if(write_plots_) {
				//Instead could initiate plots here which can be saved to the same file_
	}

		}

		return 0;
	}

	int RunIIGenJetValidation::Execute(TreeEvent *event) {

		// Get the objects we need from the event
		EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
		std::vector<GenJet*> jets = event->GetPtrVec<GenJet>("genJets");
                std::vector<GenJet*> reclusteredjets;
                if(from_miniaod_){
                  reclusteredjets = event->GetPtrVec<GenJet>("genJetsReclustered");
                }
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
             if(!from_miniaod_){
   for(UInt_t ielec=0;ielec<jets.size();++ielec){
                   genjetpt_ = jets.at(ielec)->pt();
                   genjeteta_ = jets.at(ielec)->eta();
                   genjetphi_ = jets.at(ielec)->phi();
                   genjetconstituents_ = jets.at(ielec)->n_constituents();
               slimmedgenjetpt_=-999;
              slimmedgenjeteta_=-999;
              slimmedgenjetphi_=-999;
              slimmedgenjetconstituents_=-999;
            outtree_->Fill();
             }
         }


 
                  if(from_miniaod_){
                    if(jets.size()==reclusteredjets.size()){

                  //taus:
                  for(UInt_t ielec=0;ielec<jets.size();++ielec){
                   genjetpt_ = jets.at(ielec)->pt();
                   genjeteta_ = jets.at(ielec)->eta();
                   genjetphi_ = jets.at(ielec)->phi();
                   genjetconstituents_=jets.at(ielec)->n_constituents();

                   slimmedgenjetpt_ = reclusteredjets.at(ielec)->pt();
                   slimmedgenjeteta_ = reclusteredjets.at(ielec)->eta();
                   slimmedgenjetphi_ = reclusteredjets.at(ielec)->phi();
                   slimmedgenjetconstituents_ = reclusteredjets.at(ielec)->n_constituents();
             outtree_->Fill();
             }
            } else if (jets.size()>reclusteredjets.size()){
                  for(UInt_t ielec=0;ielec<jets.size();++ielec){

                   genjetpt_ = jets.at(ielec)->pt();
                   genjeteta_ = jets.at(ielec)->eta();
                   genjetphi_ = jets.at(ielec)->phi();
                   genjetconstituents_ = jets.at(ielec)->n_constituents();
                   if(ielec<reclusteredjets.size()){

                   slimmedgenjetpt_ = reclusteredjets.at(ielec)->pt();
                   slimmedgenjeteta_ = reclusteredjets.at(ielec)->eta();
                   slimmedgenjetphi_ = reclusteredjets.at(ielec)->phi();
                   slimmedgenjetconstituents_ = reclusteredjets.at(ielec)->n_constituents();
                  } else {
                 slimmedgenjetpt_=-999;
                 slimmedgenjeteta_=-999;
                 slimmedgenjetphi_=-999;
                 slimmedgenjetconstituents_=-999;
                 }
             outtree_->Fill();
             }
             } else if (jets.size()<reclusteredjets.size()){
                  for(UInt_t ielec=0;ielec<reclusteredjets.size();++ielec){
                   slimmedgenjetpt_ = reclusteredjets.at(ielec)->pt();
                   slimmedgenjeteta_ = reclusteredjets.at(ielec)->eta();
                   slimmedgenjetphi_ = reclusteredjets.at(ielec)->phi();
                   slimmedgenjetconstituents_ = reclusteredjets.at(ielec)->n_constituents();
                   if(ielec<jets.size()){

                   genjetpt_ = jets.at(ielec)->pt();
                   genjeteta_ = jets.at(ielec)->eta();
                   genjetphi_ = jets.at(ielec)->phi();
                   genjetconstituents_ = jets.at(ielec)->n_constituents();
                  } else {
                 genjetpt_=-999;
                 genjeteta_=-999;
                 genjetphi_=-999;
                 genjetconstituents_=-999;
                 }
             outtree_->Fill();
           }
         }
        } 
               

		}
		return 0;
	}

	int RunIIGenJetValidation::PostAnalysis() {
		std::cout << "-------------------------------------" << std::endl;
		std::cout << "RunIIGenJetValidation" << std::endl;
		std::cout << "-------------------------------------" << std::endl;      
		return 0;
	}
	void RunIIGenJetValidation::PrintInfo() {
		;
	}
}


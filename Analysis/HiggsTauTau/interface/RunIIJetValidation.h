#ifndef ICHiggsTauTau_HiggsHTohh_RunIIJetValidation_h
#define ICHiggsTauTau_HiggsHTohh_RunIIJetValidation_h
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/interface/Jet.hh"
#include "TH2.h"

#include <string>

namespace ic {

	class RunIIJetValidation : public ModuleBase {
		private:

			CLASS_MEMBER(RunIIJetValidation, bool, write_tree)
			CLASS_MEMBER(RunIIJetValidation, bool, write_plots)
                        CLASS_MEMBER(RunIIJetValidation, bool, from_miniaod)
			CLASS_MEMBER(RunIIJetValidation, fwlite::TFileService*, fs)

			TTree *outtree_;
                        double jetpt_;
                        double jeteta_;
                        double jetphi_;
                        double pujetidmva_;
                        int jetpdg_;
                        double btagsimplehigh_;
                        double btagsimplepur_;
                        double btagcombined_;

                        double slimmedjetpt_;
                        double slimmedjeteta_;
                        double slimmedjetphi_;
                        double slimmedpujetidmva_;
                        int slimmedjetpdg_;
                        double slimmedbtagsimplehigh_;
                        double slimmedbtagsimplepur_;
                        double slimmedbtagcombined_;


		public:
			RunIIJetValidation(std::string const& name);
			virtual ~RunIIJetValidation();

			virtual int PreAnalysis();
			virtual int Execute(TreeEvent *event);
			virtual int PostAnalysis();
			virtual void PrintInfo();

	};

}


#endif

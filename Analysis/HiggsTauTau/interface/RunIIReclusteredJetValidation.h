#ifndef ICHiggsTauTau_HiggsHTohh_RunIIReclusteredJetValidation_h
#define ICHiggsTauTau_HiggsHTohh_RunIIReclusteredJetValidation_h
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/interface/Jet.hh"
#include "TH2.h"

#include <string>

namespace ic {

	class RunIIReclusteredJetValidation : public ModuleBase {
		private:

			CLASS_MEMBER(RunIIReclusteredJetValidation, bool, write_tree)
			CLASS_MEMBER(RunIIReclusteredJetValidation, bool, write_plots)
                        CLASS_MEMBER(RunIIReclusteredJetValidation, bool, from_miniaod)
			CLASS_MEMBER(RunIIReclusteredJetValidation, fwlite::TFileService*, fs)

			TTree *outtree_;
                        double jetpt_;
                        double jeteta_;
                        double jetphi_;
                        double pujetidmva_;
                        int jetpdg_;
                        double btagsimplehigh_;
                        double btagsimplepur_;
                        double btagcombined_;


		public:
			RunIIReclusteredJetValidation(std::string const& name);
			virtual ~RunIIReclusteredJetValidation();

			virtual int PreAnalysis();
			virtual int Execute(TreeEvent *event);
			virtual int PostAnalysis();
			virtual void PrintInfo();

	};

}


#endif

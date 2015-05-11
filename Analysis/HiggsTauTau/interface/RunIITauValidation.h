#ifndef ICHiggsTauTau_HiggsHTohh_RunIITauValidation_h
#define ICHiggsTauTau_HiggsHTohh_RunIITauValidation_h
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/interface/Tau.hh"
#include "TH2.h"

#include <string>

namespace ic {

	class RunIITauValidation : public ModuleBase {
		private:

			CLASS_MEMBER(RunIITauValidation, bool, write_tree)
			CLASS_MEMBER(RunIITauValidation, bool, write_plots)
			CLASS_MEMBER(RunIITauValidation, fwlite::TFileService*, fs)

			TTree *outtree_;
                        double taupt_;
                        double taueta_;
                        double tauphi_;
                        int taudm_;
                        double taudmf_;
                        double taudmfnew_;
                        double taucombined_;

		public:
			RunIITauValidation(std::string const& name);
			virtual ~RunIITauValidation();

			virtual int PreAnalysis();
			virtual int Execute(TreeEvent *event);
			virtual int PostAnalysis();
			virtual void PrintInfo();

	};

}


#endif

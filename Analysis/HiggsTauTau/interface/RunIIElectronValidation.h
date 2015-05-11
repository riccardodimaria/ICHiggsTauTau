#ifndef ICHiggsTauTau_HiggsHTohh_RunIIElectronValidation_h
#define ICHiggsTauTau_HiggsHTohh_RunIIElectronValidation_h
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/interface/Electron.hh"
#include "TH2.h"

#include <string>

namespace ic {

	class RunIIElectronValidation : public ModuleBase {
		private:

			CLASS_MEMBER(RunIIElectronValidation, bool, write_tree)
			CLASS_MEMBER(RunIIElectronValidation, bool, write_plots)
			CLASS_MEMBER(RunIIElectronValidation, fwlite::TFileService*, fs)

			TTree *outtree_;
                        double electronpt_;
                        double electroneta_;
                        double electronphi_;
                        double electronidmva_;
                        int electronconversion_;
                        double electronchargediso_;
                        double electronchargedisoall_;
                        double electronneutraliso_;
                        double electrongammaiso_;
                        double electronpuiso_;

		public:
			RunIIElectronValidation(std::string const& name);
			virtual ~RunIIElectronValidation();

			virtual int PreAnalysis();
			virtual int Execute(TreeEvent *event);
			virtual int PostAnalysis();
			virtual void PrintInfo();

	};

}


#endif

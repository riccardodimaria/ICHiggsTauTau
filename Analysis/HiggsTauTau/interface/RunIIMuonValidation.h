#ifndef ICHiggsTauTau_HiggsHTohh_RunIIMuonValidation_h
#define ICHiggsTauTau_HiggsHTohh_RunIIMuonValidation_h
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/interface/Muon.hh"
#include "TH2.h"

#include <string>

namespace ic {

	class RunIIMuonValidation : public ModuleBase {
		private:

			CLASS_MEMBER(RunIIMuonValidation, bool, write_tree)
			CLASS_MEMBER(RunIIMuonValidation, bool, write_plots)
			CLASS_MEMBER(RunIIMuonValidation, fwlite::TFileService*, fs)

			TTree *outtree_;
                        double muonpt_;
                        double muoneta_;
                        double muonphi_;
                        double muonchargediso_;
                        double muonchargedisoall_;
                        double muongammaiso_;
                        double muonneutraliso_;
                        double muonpuiso_;

		public:
			RunIIMuonValidation(std::string const& name);
			virtual ~RunIIMuonValidation();

			virtual int PreAnalysis();
			virtual int Execute(TreeEvent *event);
			virtual int PostAnalysis();
			virtual void PrintInfo();

	};

}


#endif

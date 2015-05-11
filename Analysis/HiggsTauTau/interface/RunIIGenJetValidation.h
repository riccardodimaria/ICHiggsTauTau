#ifndef ICHiggsTauTau_HiggsHTohh_RunIIGenJetValidation_h
#define ICHiggsTauTau_HiggsHTohh_RunIIGenJetValidation_h
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/interface/Electron.hh"
#include "TH2.h"

#include <string>

namespace ic {

	class RunIIGenJetValidation : public ModuleBase {
		private:

			CLASS_MEMBER(RunIIGenJetValidation, bool, write_tree)
			CLASS_MEMBER(RunIIGenJetValidation, bool, write_plots)
                        CLASS_MEMBER(RunIIGenJetValidation, bool, from_miniaod)
			CLASS_MEMBER(RunIIGenJetValidation, fwlite::TFileService*, fs)

			TTree *outtree_;
                        double genjetpt_;
                        double genjeteta_;
                        double genjetphi_;
                        int genjetconstituents_;
                         double slimmedgenjetpt_;
                        double slimmedgenjeteta_;
                        double slimmedgenjetphi_;
                        int slimmedgenjetconstituents_;


		public:
			RunIIGenJetValidation(std::string const& name);
			virtual ~RunIIGenJetValidation();

			virtual int PreAnalysis();
			virtual int Execute(TreeEvent *event);
			virtual int PostAnalysis();
			virtual void PrintInfo();

	};

}


#endif

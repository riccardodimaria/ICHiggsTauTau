#ifndef ICHiggsTauTau_HiggsHTohh_RunIIReclusteredGenJetValidation_h
#define ICHiggsTauTau_HiggsHTohh_RunIIReclusteredGenJetValidation_h
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/interface/Electron.hh"
#include "TH2.h"

#include <string>

namespace ic {

	class RunIIReclusteredGenJetValidation : public ModuleBase {
		private:

			CLASS_MEMBER(RunIIReclusteredGenJetValidation, bool, write_tree)
			CLASS_MEMBER(RunIIReclusteredGenJetValidation, bool, write_plots)
                        CLASS_MEMBER(RunIIReclusteredGenJetValidation, bool, from_miniaod)
			CLASS_MEMBER(RunIIReclusteredGenJetValidation, fwlite::TFileService*, fs)

			TTree *outtree_;
                        double genjetpt_;
                        double genjeteta_;
                        double genjetphi_;
                        int genjetconstituents_;


		public:
			RunIIReclusteredGenJetValidation(std::string const& name);
			virtual ~RunIIReclusteredGenJetValidation();

			virtual int PreAnalysis();
			virtual int Execute(TreeEvent *event);
			virtual int PostAnalysis();
			virtual void PrintInfo();

	};

}


#endif

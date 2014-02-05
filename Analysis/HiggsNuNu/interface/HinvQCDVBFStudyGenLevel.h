#ifndef ICHiggsTauTau_Module_HinvQCDVBFStudyGenLevel_h
#define ICHiggsTauTau_Module_HinvQCDVBFStudyGenLevel_h

#include <string>
#include <vector>

#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"

#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/HistoSet.h"

#include "UserCode/ICHiggsTauTau/interface/TH2DAsymErr.h"
#include "TH1D.h"
#include "TH2D.h"


namespace ic {

  struct HinvQCDVBFStudyGenLevelMETPlots {
    

    TH1D *gen_sumNu_Et;    
    TH1D *reco_met;

    TH2D *genVsReco_met;
    
    HinvQCDVBFStudyGenLevelMETPlots(TFileDirectory const& dir);
    
  };

  class HinvQCDVBFStudyGenLevel : public ModuleBase {

  private:

    CLASS_MEMBER(HinvQCDVBFStudyGenLevel,fwlite::TFileService*, fs);
    CLASS_MEMBER(HinvQCDVBFStudyGenLevel,std::string, met_label);
    CLASS_MEMBER(HinvQCDVBFStudyGenLevel,std::string, dijet_label);
    CLASS_MEMBER(HinvQCDVBFStudyGenLevel,std::string, genparticles_label);
    CLASS_MEMBER(HinvQCDVBFStudyGenLevel,std::string, sel_label);
    CLASS_MEMBER(HinvQCDVBFStudyGenLevel,std::string, channel);
    CLASS_MEMBER(HinvQCDVBFStudyGenLevel,bool, is_data);
    CLASS_MEMBER(HinvQCDVBFStudyGenLevel,bool, is_embedded);


    double yields_;
    HinvQCDVBFStudyGenLevelMETPlots *m_genLevelMETPlots;

  // Event Properties

    double wt_;

 
  public:
    HinvQCDVBFStudyGenLevel(std::string const& name);
    virtual ~HinvQCDVBFStudyGenLevel();
    
    virtual int PreAnalysis();
    virtual int Execute(TreeEvent *event);
    virtual int PostAnalysis();
    virtual void PrintInfo();
 
    void FillYields();
    void InitValidationPlots();
   
    void Reset();

  };

}


#endif

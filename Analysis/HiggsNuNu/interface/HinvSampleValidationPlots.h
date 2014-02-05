#ifndef ICHiggsTauTau_Module_HinvSampleValidationPlots_h
#define ICHiggsTauTau_Module_HinvSampleValidationPlots_h

#include <string>

#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"

#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/HistoSet.h"

#include "UserCode/ICHiggsTauTau/interface/TH2DAsymErr.h"
#include "TH1F.h"
#include "TH2F.h"


namespace ic {

  struct HinvValidationPlots {
    
    TH1D *JetPt1_pt;
    TH1D *JetPt2_pt;
    TH1D *JetPt1_eta;
    TH1D *JetPt2_eta;
    TH1D *JetPt1_phi;    
    TH1D *JetPt2_phi;
    
    HinvValidationPlots(TFileDirectory const& dir);
    
  };

  class HinvSampleValidationPlots : public ModuleBase {

  private:

    CLASS_MEMBER(HinvSampleValidationPlots,fwlite::TFileService*, fs);
    CLASS_MEMBER(HinvSampleValidationPlots,std::string, met_label);
    CLASS_MEMBER(HinvSampleValidationPlots,std::string, dijet_label);
    CLASS_MEMBER(HinvSampleValidationPlots,std::string, genparticles_label);
    CLASS_MEMBER(HinvSampleValidationPlots,std::string, sel_label);
    CLASS_MEMBER(HinvSampleValidationPlots,std::string, channel);
    CLASS_MEMBER(HinvSampleValidationPlots,bool, is_data);
    CLASS_MEMBER(HinvSampleValidationPlots,bool, is_embedded);


    double yields_;
    HinvValidationPlots *m_valPlots;

  // Event Properties

    double wt_;

 
  public:
    HinvSampleValidationPlots(std::string const& name);
    virtual ~HinvSampleValidationPlots();
    
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

#ifndef ICHiggsTauTau_Module_HTTRun2RecoilCorrector_h
#define ICHiggsTauTau_Module_HTTRun2RecoilCorrector_h

#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/RecoilCorrectorRun2.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include <string>

namespace ic {

class HTTRun2RecoilCorrector : public ModuleBase {
 private:
 // std::string jets_label_;
  bool disable;
  bool is_ztt;
  bool is_wjets;
  RecoilCorrectorRun2 * corrector_;
  CLASS_MEMBER(HTTRun2RecoilCorrector, std::string, dilepton_label)
  CLASS_MEMBER(HTTRun2RecoilCorrector, std::string, met_label)
  CLASS_MEMBER(HTTRun2RecoilCorrector, std::string, jets_label)
  CLASS_MEMBER(HTTRun2RecoilCorrector, std::string, sample)
  CLASS_MEMBER(HTTRun2RecoilCorrector, ic::channel, channel)
  CLASS_MEMBER(HTTRun2RecoilCorrector, ic::strategy, strategy)
  CLASS_MEMBER(HTTRun2RecoilCorrector, ic::mc, mc)
  CLASS_MEMBER(HTTRun2RecoilCorrector, ic::era, era)



 public:
  HTTRun2RecoilCorrector(std::string const& name);
  virtual ~HTTRun2RecoilCorrector();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();

};

}


#endif
#ifndef ICHiggsTauTau_Module_QCDStudy_h
#define ICHiggsTauTau_Module_QCDStudy_h

#include "Core/interface/TreeEvent.h"
#include "Core/interface/ModuleBase.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Utilities/interface/HistoSet.h"
#include "HiggsTauTau/interface/HTTConfig.h"
#include "HiggsTauTau/interface/HTTTriggerFilter.h"
#include "boost/range/algorithm_ext/erase.hpp"

#include <string>
#include <functional>


namespace ic {



template <class T, class U>
std::vector<CompositeCandidate*> & CompositeProducer(
    TreeEvent* event, std::vector<T*> const& vec_first,
    std::vector<U*> const& vec_second, std::string output) {
  std::vector<std::pair<T*, U*> > pairs = MakePairs(vec_first, vec_second);
  std::vector<CompositeCandidate> vec_out;
  std::vector<CompositeCandidate*> ptr_vec_out;
  for (unsigned i = 0; i < pairs.size(); ++i) {
    vec_out.push_back(CompositeCandidate());
    CompositeCandidate& cand_ref = vec_out.back();
    cand_ref.AddCandidate("leg0", pairs[i].first);
    cand_ref.AddCandidate("leg1", pairs[i].second);
  }
  event->Add(output + "Product", vec_out);
  std::vector<CompositeCandidate>& vec_in =
      event->Get<std::vector<CompositeCandidate> >(output + "Product");
  ptr_vec_out.resize(vec_in.size());
  for (unsigned i = 0; i < vec_in.size(); ++i) {
    ptr_vec_out[i] = &(vec_in[i]);
  }
  event->Add(output, ptr_vec_out);
  return event->GetPtrVec<CompositeCandidate>(output);
}

template <class Container, class Pred>
Container copy_keep_if(Container& target, Pred pred) {
  Container res = target;
  keep_if(res, pred);
  return res;
}

template <class Container, class Pred>
Container& keep_if(Container& target, Pred pred) {
  return boost::remove_erase_if(target, [&](typename Container::value_type const& x) {
    return !pred(x);
  });
}

class QCDStudy : public ModuleBase {
 private:
  CLASS_MEMBER(QCDStudy, ic::channel, channel)
  CLASS_MEMBER(QCDStudy, fwlite::TFileService*, fs)

  TTree * tree_;

  float pt_1_ = 0.;
  float pt_2_ = 0.;
  float eta_1_ = 0.;
  float eta_2_ = 0.;
  float phi_1_ = 0.;
  float phi_2_ = 0.;
  int q_1_ = 0;
  int q_2_ = 0;
  float iso_1_ = 0.;
  float iso_2_ = 0.;
  bool matched_jet_1_ = false;
  bool matched_jet_2_ = false;
  int flavour_1_ = 0;
  int flavour_2_ = 0;

  bool anti_mu_2_ = false;


 private:
  static bool ETLooseEle(Electron const* e);
  static bool ETLooseTau(Tau const* t);
  static bool ETTightEle(Electron const* e);
  static bool ETTightTau(Tau const* t);

  static bool MTLooseMuo(Muon const* m);
  static bool MTLooseTau(Tau const* t);
  static bool MTTightMuo(Muon const* m);
  static bool MTTightTau(Tau const* t);

  static bool EMLooseEle(Electron const* e);
  static bool EMLooseMuo(Muon const* m);
  static bool EMTightEle(Electron const* e);
  static bool EMTightMuo(Muon const* m);

  static bool TTLooseTau(Tau const* t);
  static bool TTTightTau(Tau const* t);


  static bool SortByOSThenSumPt(CompositeCandidate const* c1,
                         CompositeCandidate const* c2);
  static bool SuperSortET(CompositeCandidate const* c1, CompositeCandidate const* c2);
  static bool SuperSortMT(CompositeCandidate const* c1, CompositeCandidate const* c2);
  static bool SuperSortEM(CompositeCandidate const* c1, CompositeCandidate const* c2);
  static bool SuperSortTT(CompositeCandidate const* c1, CompositeCandidate const* c2);

 public:
  QCDStudy(std::string const& name);
  virtual ~QCDStudy();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();

};

}

#endif

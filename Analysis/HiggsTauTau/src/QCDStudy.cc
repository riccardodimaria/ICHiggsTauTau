#include "HiggsTauTau/interface/QCDStudy.h"
#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "Utilities/interface/FnPredicates.h"
#include "Utilities/interface/FnPairs.h"
#include "HiggsTauTau/interface/HTTGenEvent.h"
#include <boost/functional/hash.hpp>
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/format.hpp"
#include "boost/range/algorithm_ext/erase.hpp"

namespace ic {

QCDStudy::QCDStudy(std::string const& name)
    : ModuleBase(name), channel_(ic::channel::et), fs_(nullptr) {}

QCDStudy::~QCDStudy() { ; }

int QCDStudy::PreAnalysis() {
  if (fs_) {
    tree_ = fs_->make<TTree>("qcdstudy", "qcdstudy");
    tree_->Branch("pt_1", &pt_1_);
    tree_->Branch("pt_2", &pt_2_);
    tree_->Branch("eta_1", &eta_1_);
    tree_->Branch("eta_2", &eta_2_);
    tree_->Branch("phi_1", &phi_1_);
    tree_->Branch("phi_2", &phi_2_);
    tree_->Branch("q_1", &q_1_);
    tree_->Branch("q_2", &q_2_);
    tree_->Branch("iso_1", &iso_1_);
    tree_->Branch("iso_2", &iso_2_);
    tree_->Branch("matched_jet_1", &matched_jet_1_);
    tree_->Branch("matched_jet_2", &matched_jet_2_);
    tree_->Branch("flavour_1", &flavour_1_);
    tree_->Branch("flavour_2", &flavour_2_);
    tree_->Branch("anti_mu_2", &anti_mu_2_);
    tree_->Branch("nbtags", &nbtags_);
    tree_->Branch("flavour_btag", &nbtags_);
    // TFileDirectory subdir = fs_->mkdir(this->ModuleName());
  }
  return 0;
}

int QCDStudy::Execute(TreeEvent* event) {
  if (channel_ == channel::mt) {
    // std::vector<Muon const*> sel_muo;
    // std::vector<Tau const*> sel_tau;
    ic::CompositeCandidate* ditau = nullptr;
    auto muons = event->GetPtrVec<Muon>("muonsPFlow");
    auto taus = event->GetPtrVec<Tau>("taus");

    auto muo_loose = copy_keep_if(muons, MTLooseMuo);
    auto tau_loose = copy_keep_if(taus, MTLooseTau);
    auto& mt_pairs =
        CompositeProducer(event, muo_loose, tau_loose, "ditau_loose_mt");
    ic::keep_if(mt_pairs, [&](CompositeCandidate const* c) {
      return DR(c->At(0), c->At(1)) > 0.5;
    });
      std::sort(mt_pairs.begin(), mt_pairs.end(), SuperSortMT);
   if (mt_pairs.size() == 0) return 1;

   ditau = mt_pairs[0];
   auto sel_muo = static_cast<Muon const*>(ditau->At(0));
   auto sel_tau = static_cast<Tau const*>(ditau->At(1));

   pt_1_ = sel_muo->pt();
   pt_2_ = sel_tau->pt();
   eta_1_ = sel_muo->eta();
   eta_2_ = sel_tau->eta();
   phi_1_ = sel_muo->phi();
   phi_2_ = sel_tau->phi();
   q_1_ = sel_muo->charge();
   q_2_ = sel_tau->charge();
   iso_1_ = PF04IsolationVal(sel_muo, 0.5);
   iso_2_ = sel_tau->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
   anti_mu_2_ = sel_tau->GetTauID("againstMuonTight") > 0.5;

   auto jets = event->GetPtrVec<ic::PFJet>("pfJetsPFlow");

   std::vector<ic::Muon const*> sel_muo_vec = {sel_muo};
   std::vector<ic::Tau const*> sel_tau_vec = {sel_tau};
   std::vector<ic::Candidate const*> cands = {sel_muo, sel_tau};

   auto muo_jet_matches = ic::MatchByDR(sel_muo_vec, jets, 0.5, true, true);
   auto tau_jet_matches = ic::MatchByDR(sel_tau_vec, jets, 0.5, true, true);

   if (muo_jet_matches.size() > 0) {
    matched_jet_1_ = true;
    flavour_1_ = muo_jet_matches[0].second->parton_flavour();
   }
   if (tau_jet_matches.size() > 0) {
    matched_jet_2_ = true;
    flavour_2_ = tau_jet_matches[0].second->parton_flavour();
   }

   auto btag_jets = ic::copy_keep_if(jets, [&](ic::PFJet* jet) {
     return ic::MinPtMaxEta(jet, 20., 2.4) &&
            jet->GetBDiscriminator("combinedSecondaryVertexBJetTags") > 0.679 &&
            ic::MinDRToCollection(jet, cands, 0.5);
   });

   nbtags_ = btag_jets.size();

   tree_->Fill();
  }
  return 0;
}
int QCDStudy::PostAnalysis() {
  std::cout << "*** " << this->ModuleName() << "\n";
  return 0;
}

void QCDStudy::PrintInfo() {}

// ET selection
// ----------------------------------------------------------------
bool QCDStudy::ETLooseEle(Electron const* e) {
  return e->pt() > 24. && std::fabs(e->eta()) < 2.1 &&
         std::fabs(e->dxy_vertex()) < 0.045 &&
         std::fabs(e->dz_vertex()) < 0.2 && ElectronHTTId(e, false);
}

bool QCDStudy::ETTightEle(Electron const* e) {
  return ETLooseEle(e) && PF04IsolationVal(e, 0.5) < 0.1;
}

bool QCDStudy::ETLooseTau(Tau const* t) {
  return t->pt() > 20. && std::fabs(t->eta()) < 2.3 &&
         std::fabs(t->lead_dz_vertex()) < 0.2 &&
         t->GetTauID("decayModeFinding") > 0.5;
}

bool QCDStudy::ETTightTau(Tau const* t) {
  return ETLooseTau(t) &&
         t->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 1.5 &&
         passAntiEMVA(t, 0) && t->GetTauID("againstMuonLoose") > 0.5;
}

// MT selection
// ----------------------------------------------------------------
bool QCDStudy::MTLooseMuo(Muon const* m) {
  return m->pt() > 20. && std::fabs(m->eta()) < 2.1 &&
         std::fabs(m->dxy_vertex()) < 0.045 &&
         std::fabs(m->dz_vertex()) < 0.2 && MuonTight(m);
}

bool QCDStudy::MTTightMuo(Muon const* m) {
  return MTLooseMuo(m) && PF04IsolationVal(m, 0.5) < 0.1;
}

bool QCDStudy::MTLooseTau(Tau const* t) {
  return t->pt() > 20. && std::fabs(t->eta()) < 2.3 &&
         std::fabs(t->lead_dz_vertex()) < 0.2 &&
         t->GetTauID("decayModeFinding") > 0.5;
}

bool QCDStudy::MTTightTau(Tau const* t) {
  return MTLooseTau(t) &&
         t->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < 1.5 &&
         t->GetTauID("againstElectronLoose") && TauEoverP(t, 0.2) &&
         t->GetTauID("againstMuonTight") > 0.5;
}

// EM selection
// ----------------------------------------------------------------
bool QCDStudy::EMLooseEle(Electron const* e) {
  return e->pt() > 10. && std::fabs(e->eta()) < 2.3 &&
         std::fabs(e->dxy_vertex()) < 0.02 && std::fabs(e->dz_vertex()) < 0.1 &&
         ElectronHTTId(e, true);
}

bool QCDStudy::EMTightEle(Electron const* e) {
  return EMLooseEle(e) && PF04IsolationEBElec(e, 0.5, 0.15, 0.1);
}

bool QCDStudy::EMLooseMuo(Muon const* m) {
  return m->pt() > 10. && std::fabs(m->eta()) < 2.1 &&
         std::fabs(m->dxy_vertex()) < 0.02 && std::fabs(m->dz_vertex()) < 0.1 &&
         MuonTight(m);
}

bool QCDStudy::EMTightMuo(Muon const* m) {
  return EMLooseMuo(m) && PF04IsolationEB(m, 0.5, 0.15, 0.1);
}

// TT selection
// ----------------------------------------------------------------
bool QCDStudy::TTLooseTau(Tau const* t) {
  return t->pt() > 45. && std::fabs(t->eta()) < 2.1 &&
         std::fabs(t->lead_dz_vertex()) < 0.2 &&
         t->GetTauID("decayModeFinding") > 0.5;
}

bool QCDStudy::TTTightTau(Tau const* t) {
  return TTLooseTau(t) &&
         t->GetTauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") > 0.5 &&
         t->GetTauID("againstMuonLoose") > 0.5 &&
         t->GetTauID("againstElectronLoose") > 0.5;
}

// Sorting
// ----------------------------------------------------------------
bool QCDStudy::SortByOSThenSumPt(CompositeCandidate const* c1,
                       CompositeCandidate const* c2) {
  if (c1->charge() == 0 && c2->charge() != 0) return true;
  if (c1->charge() != 0 && c2->charge() == 0) return false;
  return ScalarPtSum(c1->AsVector()) > ScalarPtSum(c2->AsVector());
}

bool QCDStudy::SuperSortET(CompositeCandidate const* c1, CompositeCandidate const* c2) {
  // First we sort the electrons
  Electron const* e1 = static_cast<Electron const*>(c1->At(0));
  Electron const* e2 = static_cast<Electron const*>(c2->At(0));
  double e_iso1 = PF04IsolationVal(e1, 0.5);
  double e_iso2 = PF04IsolationVal(e2, 0.5);
  // If the iso is different we just use this
  if (e_iso1 != e_iso2) return e_iso1 < e_iso2;
  // If not try the pT
  if (e1->pt() != e2->pt()) return e1->pt() > e2->pt();
  // If both of these are the same then try the taus
  Tau const* t1 = static_cast<Tau const*>(c1->At(1));
  Tau const* t2 = static_cast<Tau const*>(c2->At(1));
  double t_iso1 = t1->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  double t_iso2 = t2->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  if (t_iso1 != t_iso2) return t_iso1 < t_iso2;
  return (t1->pt() > t2->pt());
}

bool QCDStudy::SuperSortMT(CompositeCandidate const* c1, CompositeCandidate const* c2) {
  // First we sort the electrons
  Muon const* m1 = static_cast<Muon const*>(c1->At(0));
  Muon const* m2 = static_cast<Muon const*>(c2->At(0));
  double m_iso1 = PF04IsolationVal(m1, 0.5);
  double m_iso2 = PF04IsolationVal(m2, 0.5);
  // If the iso is different we just use this
  if (m_iso1 != m_iso2) return m_iso1 < m_iso2;
  // If not try the pT
  if (m1->pt() != m2->pt()) return m1->pt() > m2->pt();
  // If both of these are the same then try the taus
  Tau const* t1 = static_cast<Tau const*>(c1->At(1));
  Tau const* t2 = static_cast<Tau const*>(c2->At(1));
  double t_iso1 = t1->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  double t_iso2 = t2->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  if (t_iso1 != t_iso2) return t_iso1 < t_iso2;
  return (t1->pt() > t2->pt());
}

bool QCDStudy::SuperSortEM(CompositeCandidate const* c1, CompositeCandidate const* c2) {
  // First we sort the muons
  Muon const* m1 = static_cast<Muon const*>(c1->At(1));
  Muon const* m2 = static_cast<Muon const*>(c2->At(1));
  double m_iso1 = PF04IsolationVal(m1, 0.5);
  double m_iso2 = PF04IsolationVal(m2, 0.5);
  // If the iso is different we just use this
  if (m_iso1 != m_iso2) return m_iso1 < m_iso2;
  // If not try the pT
  if (m1->pt() != m2->pt()) return m1->pt() > m2->pt();
  // If both of these are the same then try the electrons
  Electron const* e1 = static_cast<Electron const*>(c1->At(0));
  Electron const* e2 = static_cast<Electron const*>(c2->At(0));
  double e_iso1 = PF04IsolationVal(e1, 0.5);
  double e_iso2 = PF04IsolationVal(e2, 0.5);
  if (e_iso1 != e_iso2) return e_iso1 < e_iso2;
  return e1->pt() > e2->pt();
}

bool QCDStudy::SuperSortTT(CompositeCandidate const* c1, CompositeCandidate const* c2) {
  // First we sort the electrons
  Tau const* t1_1 = static_cast<Tau const*>(c1->At(0));
  Tau const* t1_2 = static_cast<Tau const*>(c2->At(0));
  double t1_iso1 = t1_1->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  double t1_iso2 = t1_2->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  // If the iso is different we just use this
  if (t1_iso1 != t1_iso2) return t1_iso1 < t1_iso2;
  // If not try the pT
  if (t1_1->pt() != t1_2->pt()) return t1_1->pt() > t1_2->pt();
  // If both of these are the same then try the taus
  Tau const* t2_1 = static_cast<Tau const*>(c1->At(1));
  Tau const* t2_2 = static_cast<Tau const*>(c2->At(1));
  double t2_iso1 = t2_1->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  double t2_iso2 = t2_2->GetTauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  // If the iso is different we just use this
  if (t2_iso1 != t2_iso2) return t2_iso1 < t2_iso2;
  // If not try the pT
  return t2_1->pt() > t2_2->pt();
}
}

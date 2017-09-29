#ifndef ICHiggsTauTau_Utilities_FnPredicates_h
#define ICHiggsTauTau_Utilities_FnPredicates_h
#include <algorithm>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <set>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include "boost/range/algorithm_ext/erase.hpp"
#include "Math/VectorUtil.h"
#include "TVector3.h"

#include "UserCode/ICHiggsTauTau/interface/Objects.hh"
#include "UserCode/ICHiggsTauTau/interface/SuperCluster.hh"
#include "UserCode/ICHiggsTauTau/interface/CompositeCandidate.hh"

namespace ic {

  template <class T, class U> void erase_if(T & t, U pred) {
    t.erase(std::remove_if(t.begin(), t.end(), pred), t.end());
  }

  template <class T, class U>
  void erase_if_not(T& t, U pred) {
    t.erase(std::remove_if(
                t.begin(), t.end(),
                [&](typename T::value_type const& x) { return !pred(x); }),
            t.end());
  }

  /**
   * Remove elements from a container that return false from a given predicate
   */
  template <class Container, class Pred>
  Container& keep_if(Container& target, Pred pred) {
    return boost::remove_erase_if(target, [&](typename Container::value_type const& x) {
      return !pred(x);
    });
  }

  /**
   * Copy a container, filter elements, then return this copy
   */
  template <class Container, class Pred>
  Container copy_keep_if(Container& target, Pred pred) {
    Container res = target;
    keep_if(res, pred);
    return res;
  }

  //Dummy
  template <class T> bool DummyFunction(T const* cand ){
    (void)cand;
    return true;
  }


  //----------------------------------------------------------
  // Candidate
  //----------------------------------------------------------
  bool MinPtMaxEta(Candidate const* cand, double const& minPt, double const& maxEta);
  bool MinPtMaxSCEta(Electron const* cand, double const& minPt, double const& maxEta);
  bool EtaOutsideRange(Candidate const* cand, double const& minEtaCut, double const& maxEtaCut);
  bool VertexDz(Tau const* cand, double const& vertexZ);
  //----------------------------------------------------------
  // Vertex
  //----------------------------------------------------------
  bool GoodVertex(Vertex const* vertex);

  //
  bool isTightMuon(Muon* veto,std::vector<Muon*> sel);
  bool isTightElectron(Electron* veto,std::vector<Electron*> sel);
  bool isTightTau(Tau* veto,std::vector<Tau*> sel);

  //----------------------------------------------------------
  // Electron
  //----------------------------------------------------------
  bool InEcalGap(Electron const* electron);
  bool MaxEtaSC(Electron const* electron, double const& maxEta);
  bool passAntiEMVA(Tau const * tau, int WP);
  bool TauEoverP(Tau const* tau, double cut);
  bool Electron2011WP85Iso(Electron const* elec, double const& rho);
  bool Electron2011WP85IsoNoRhoSub(Electron const* elec);
  bool Electron2011WP70ID(Electron const* elec);
  bool Electron2011WP85ID(Electron const* elec);
  bool Electron2011WP95ID(Electron const* elec);
  bool TightElectronIDSpring15(Electron const* elec);
  bool VetoElectronIDPhys14(Electron const* elec);
  bool VetoElectronIDSpring15(Electron const* elec);
  bool TightElectronFullIDSpring15(Electron const* elec, double const& rho);
  bool VetoElectronFullIDSpring15(Electron const* elec, double const& rho);
  bool VetoElectronID(Electron const* elec);
  bool VetoElectronIDSpring16(Electron const* elec);
  bool ElectronSimpleWP85Iso(Electron const* elec);
  bool ElectronSimpleWP85ID(Electron const* elec);
  bool ElectronHTTIdIso(Electron const* elec, unsigned const& mode);
  bool ElectronHTTVHID(Electron const* elec);
  bool ElectronHTTId(Electron const* elec, bool loose_wp);
  bool ElectronHTTIdPhys14(Electron const* elec, bool loose_wp);
  bool ElectronHTTIdSpring15(Electron const* elec, bool loose_wp);
  bool ElectronHTTIdTrigSpring15(Electron const* elec, bool loose_wp);
  bool ElectronHTTTrigNoIPId(Electron const* elec, bool loose_wp);
  bool ElectronHTTIdSpring16(Electron const* elec, bool loose_wp);

  bool TightElectronFullID16(Electron const* elec, double const& rho);
  bool VetoElectronFullID16(Electron const* elec, double const& rho);

  double getTotalEA(const double & eta);
  double getTotalEA2016(const double & eta);
  std::pair<double,double> getEA(const double & eta);
  double getEA2016(const double & eta, const unsigned type);
  bool TightPhotonIDSpring15(Photon const* photon,double const& rho);
  bool MediumPhotonIDSpring15(Photon const* photon,double const& rho);
  bool LoosePhotonIDSpring15(Photon const* photon,double const& rho);

  bool LoosePhotonID2016(Photon const* photon,double const& rho);
  bool MediumPhotonID2016(Photon const* photon,double const& rho);
  bool TightPhotonID2016(Photon const* photon,double const& rho);

  bool HttEMuFakeElectron(Electron const* elec);
  bool HttEMuFakeMuon(Muon const* muon);
  bool HttEMuFakeMuon2011(Muon const* muon);

  bool ElectronHTT2011IdIso(Electron const* elec);

  bool ElectronZbbID(Electron const * elec);
  bool MVAMETElectron(Electron const * elec);

  bool GenParticleInMassBand(GenParticle const* part, int pdgid, double m_low, double m_high);




  bool MinPtMaxEtaSC(SuperCluster const* cand, double const& minPt, double const& maxEta);

  bool InEcalGapSC(SuperCluster const* electron);



  bool CaloJetID(CaloJet const* jet);

  // Standard particle-flow jet id
  bool PFJetID(PFJet const* jet);
  // Standard particle-flow jet id for 2015
  bool PFJetID2015(PFJet const* jet);
  // Standard particle-flow jet id for 2016
  bool PFJetID2016(PFJet const* jet);
  bool PFJetID2016v2(PFJet const* jet);


  // Particle-flow jet id without the HF energy in the neutral energy cut
  bool PFJetIDNoHFCut(PFJet const* jet);
  bool PUJetID(PFJet const* jet, bool is_2012);
  // Pileup jet id cuts.  0 = 2011, 1 = 2012, 2 = 2012 (december training), 3 = 2015 training
  //loose working points
  bool PileupJetID(PFJet const* jet, unsigned training, bool doTight=false);

  bool IsReBTagged(Jet const* jet, std::map<std::size_t, bool> const& tag_map);

  bool CSVMediumWP(Jet const* cand, double const& minDiscri);

  double PZeta(CompositeCandidate const* cand, Candidate const* met, double const& alpha);
  double PZetaVis(CompositeCandidate const* cand);

  double MT(Candidate const* cand1, Candidate const* cand2);

  bool IsFilterMatched(Candidate const* cand, std::vector<TriggerObject *> const& objs, std::string const& filter, double const& max_dr);
  std::pair <bool,unsigned> IsFilterMatchedWithIndex(Candidate const* cand, std::vector<TriggerObject *> const& objs, std::string const& filter, double const& max_dr);
  std::pair <bool,std::vector<unsigned>> IsFilterMatchedWithMultipleIndexs(Candidate const* cand, std::vector<TriggerObject *> const& objs, std::string const& filter, double const& max_dr);

  template<class T>
  double PF04IsolationVal(T const* cand, double const& dbeta, bool allcharged) {
	  double charged_iso = allcharged ? cand->dr04_pfiso_charged_all() : cand->dr04_pfiso_charged();
    double iso =  charged_iso
                  + std::max(cand->dr04_pfiso_neutral() + cand->dr04_pfiso_gamma() - dbeta * cand->dr04_pfiso_pu(), 0.0);
    iso = iso / cand->pt();
    return iso;
  }

  template<class T>
  double PF04IsolationVal(T const* cand, double const& dbeta) {
    double iso =  cand->dr04_pfiso_charged_all()
                  + std::max(cand->dr04_pfiso_neutral() + cand->dr04_pfiso_gamma() - dbeta * cand->dr04_pfiso_pu(), 0.0);
    iso = iso / cand->pt();
    return iso;
  }

  
  template<class T>
  bool PF04Isolation(T const* cand, double const& dbeta, double const& cut) {
    double iso =  PF04IsolationVal(cand, dbeta);
    return (iso < cut);
  }
  bool PF04IsolationEBElec(Electron const* cand, double const& dbeta, double const& cut_barrel, double const& cut_endcap);

  template<class T>
  bool PF04IsolationEB(T const* cand, double const& dbeta, double const& cut_barrel, double const& cut_endcap) {
    double iso =  PF04IsolationVal(cand, dbeta);
    if (fabs(cand->eta()) < 1.479) {
      return (iso < cut_barrel);
    } else {
      return (iso < cut_endcap);
    }
  }

  double PUW03IsolationVal(Muon const* muon);
  double PUW04IsolationVal(Muon const* muon);
  
  template<class T>
  double PF03IsolationVal(T const* cand, double const& dbeta, bool allcharged=false) {
    double charged_iso = allcharged ? cand->dr03_pfiso_charged_all() : cand->dr03_pfiso_charged();
    double iso =  charged_iso 
                  + std::max(cand->dr03_pfiso_neutral() + cand->dr03_pfiso_gamma() - dbeta * cand->dr03_pfiso_pu(), 0.0);
    iso = iso / cand->pt();
    return iso;
  }
  template<class T>
  bool PF03Isolation(T const* cand, double const& dbeta, double const& cut, bool allcharged=false) {
    double iso =  PF03IsolationVal(cand, dbeta, allcharged);
    return (iso < cut);
  }

  template<class T> 
  double GetEffectiveArea(T const* cand){
    return 0;
  }
  
  inline double GetEffectiveArea(Electron const* cand){
    double cand_eta = cand->eta();
    using std::abs;
//From   https://github.com/ikrav/cmssw/blob/egm_id_747_v2/RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt
    if(abs(cand_eta)<1.) return 0.1752;
    else if(abs(cand_eta)<1.479) return 0.1862;
    else if(abs(cand_eta)<2.) return 0.1411;
    else if(abs(cand_eta)<2.2) return 0.1534;
    else if(abs(cand_eta)<2.3) return 0.1903;
    else if(abs(cand_eta)<2.4) return 0.2243;
    else if(abs(cand_eta)<5.) return 0.2687;
    return 0;

  }

  inline double GetEffectiveArea(Muon const* cand){
     double cand_eta = cand->eta();
     using std::abs;
    //From https://indico.cern.ch/event/451110/contribution/1/attachments/1165046/1679258/miniIsoEffectiveAreas_muonPOG_051015.pdf
    if(abs(cand_eta)<0.8) return 0.0735;
    else if(abs(cand_eta)<1.3) return 0.0619;
    else if(abs(cand_eta)<2.) return 0.0465;
    else if(abs(cand_eta)<2.2) return 0.0443;
    else if(abs(cand_eta)<2.5) return 0.0577;
    return 0;
 
  }

  template<class T>
  double PF03EAIsolationVal(T const* cand, EventInfo const* evt, bool jet_rho=true){
  double charged_iso = cand->dr03_pfiso_charged();
    double eff_area = GetEffectiveArea(cand);
    double rho=0;
    if(jet_rho){
      evt->jet_rho();
    } else {
      evt->lepton_rho();
    }
    double iso =  charged_iso 
                  + std::max(cand->dr03_pfiso_neutral() + cand->dr03_pfiso_gamma() - rho*eff_area, 0.0);
    iso = iso / cand->pt();
    return iso;
   }


  template<class T>
  bool PF03EAIsolation(T const* cand, EventInfo const* evt, double const& cut, bool jet_rho=true) {
    double iso =  PF03EAIsolationVal(cand, evt, jet_rho);
    return (iso < cut);
  }




  bool ElectronZbbIso(Electron const* elec, bool is_data, double const& rho, double const& cut);

  template <class T>
  bool DBeamspotLessThan(T const* cand, double const& cut){
    return (fabs(cand->dxy_beamspot()) < cut);
  }

  template <class T>
  double DeltaZ(T const* cand, double const& vz){
    return fabs(cand->vz() - vz);
  }

  bool PairPtSelection(CompositeCandidate const* cand, double const& jetpt1, double const& jetpt2);
  bool OrderedPairPtSelection(CompositeCandidate const* cand, double const& jetpt1, double const& jetpt2, bool const& aboveorbelow);
  bool PairMassInRange(CompositeCandidate const* cand, double const& mLow, double const& mHigh);
  bool PairEtaProdLessThan(CompositeCandidate const* cand, double const& max);
  bool PairDEtaLessThan(CompositeCandidate const* cand, double const& max);
  bool PairDPhiLessThan(CompositeCandidate const* cand, double const& max);
  bool PairAbsDPhiLessThan(CompositeCandidate const* cand, double const& max);

  bool PairOneWithPt(CompositeCandidate const* cand, double const& ptMin);

  bool PairOppCharge(CompositeCandidate const& cand);
  
  bool PairOppSign(CompositeCandidate const* cand);
  bool PairSameSign(CompositeCandidate const* cand);
  double PairMass(CompositeCandidate const* cand);
  double PairMassDiff(CompositeCandidate const* cand, double const &mass);

  bool MuonTight(Muon const* muon);
  bool MuonMedium(Muon const* muon);
  bool MuonMediumHIPsafe(Muon const* muon);
  bool MuonLoose(Muon const* muon);
  bool MuonIso(Muon const* muon);
  bool MuonTkIso(Muon const* muon);
  double MuonTkIsoVal(Muon const* muon);


  template<class T, class U>
  bool EtaProdLessThan(std::pair<T,U> const& p, double const& max) {
    if ((p.first->eta() * p.second->eta()) < max) return true;
    return false;
  }

  template<class T, class U>
  double EtaProd(std::pair<T,U> const& p) {
    return (p.first->eta() * p.second->eta()) ;
  }

  template<class T, class U>
  bool DEtaLessThan(std::pair<T,U> const& p, double const& max) {
    if (fabs(p.first->eta() - p.second->eta()) < max) return true;
    return false;
  }

  template<class T, class U>
  bool DPhiLessThan(std::pair<T,U> const& p, double const& max) {
    if (ROOT::Math::VectorUtil::DeltaPhi(p.first->vector(), p.second->vector()) < max) return true;
    return false;
  }

  template<class T, class U>
  bool MassLessThan(std::pair<T,U> const& p, double const& max) {
    if ((p.first->vector() + p.second->vector()).M() < max) return true;
    return false;
  }

  template<class T> 
  bool MinDRToCollection(Candidate const* cand, 
    std::vector<T> const& coll, double const& cut) {
    BOOST_FOREACH(T const& ele, coll) {
      if (ROOT::Math::VectorUtil::DeltaR(cand->vector(), ele->vector()) < cut) return false;
    }
    return true;
  }


  template<class T, class U> 
  bool FoundIdInCollection(T const* cand, 
    std::vector<U *> const& coll) {
    BOOST_FOREACH(U * ele, coll) {
      if (ele->id() == cand->id()) return true;
    }
    return false;
  }

  template<class T, class U>
  bool DRLessThan(std::pair<T,U> const& cPair, double const& dr) {
    if (DR(cPair.first,cPair.second) < dr) return true;
    return false;
  }

  template<class T, class U>
  double DR(T const& cand1, U const& cand2) {
    return ROOT::Math::VectorUtil::DeltaR(cand1->vector(),cand2->vector());
  }

  template<class T, class U>
  bool DRCompare(std::pair<T,U> const& p1, std::pair<T,U> const& p2) {
      double d1 = DR (p1.first,p1.second);
      double d2 = DR (p2.first,p2.second);
      return (d1 < d2);
  }

  template<class T>
  double ScalarPtSum(std::vector<T> const& vec) {
    double result = 0.0;
    for (unsigned i = 0; i < vec.size(); ++i) {
      result += vec[i]->pt();
    }
    return result;
  }

  template<class T>
  double ScalarEnergySum(std::vector<T> const& vec) {
    double result = 0.0;
    for (unsigned i = 0; i < vec.size(); ++i) {
      result += vec[i]->energy();
    }
    return result;
  }


  std::vector<Track *> GetTracksAtVertex(std::vector<Track *> const& trks, std::vector<Vertex *> const& vtxs, unsigned idx, double const& dz);
  

  bool MassDiffCompare(Candidate const* p1, Candidate const* p2, double const& mass);

  std::vector<GenParticle *> ExtractStableDaughters(GenParticle * part, std::vector<GenParticle *> const& input);

  std::vector<GenParticle *> ExtractMothers(GenParticle * part, std::vector<GenParticle *> const& input);

  std::vector<GenParticle *> ExtractDaughters(GenParticle * part, std::vector<GenParticle *> const& input);

  std::vector<GenParticle *> ExtractDaughtersRecursive(GenParticle * part, std::vector<GenParticle *> const& input);

  std::vector<GenJet> BuildTauJets(std::vector<GenParticle *> const& parts, bool include_leptonic, bool use_prompt);

  ROOT::Math::PtEtaPhiEVector reconstructWboson(Candidate const*  lepton, Candidate const* met);

  template <class T, class U>
    void getGenRecoMatches(const std::vector<T*> & recovec,
			   const std::vector<U*> & genvec, 
			   std::vector<std::pair<unsigned,bool> > & recotogenmatch,
			   const double dRmin=1000){
    recotogenmatch.resize(recovec.size(),std::pair<unsigned,bool>(1000,false));
    unsigned nReco = recovec.size();
    unsigned nGen = genvec.size();
    unsigned nRecNotMatched = 0;
    std::vector<bool> ignoreIndex;
    ignoreIndex.resize(nGen,false);
    for (unsigned i = 0; i < nReco; ++i) {//loop on recovec
      double mindR = 1000;
      unsigned genvecid=1000;
      for (unsigned ig = 0; ig < nGen; ++ig) {//loop on genvec
	if (ignoreIndex[ig]) continue;
	double dR = ROOT::Math::VectorUtil::DeltaR(genvec[ig]->vector(),recovec[i]->vector());
	if (dR<mindR){
	  mindR = dR;
	  genvecid = ig;
	}
      }//loop on genvec
      if (genvecid==1000){
	nRecNotMatched++;
	continue;
      }
      recotogenmatch[i]=std::pair<unsigned,bool>(genvecid,mindR<dRmin?true:false);
      ignoreIndex[genvecid] = true;
      //loop again to find minimum dpT/pT
      /*double ptgen = genvec[genvecid]->pt();
      double mindpT = fabs(recovec[i]->pt()-ptgen)/ptgen;
      for (unsigned j = 0; j < i; ++j) {//loop on recovec
	if (recotogenmatch[j].first!=genvecid) continue;
	double dR = ROOT::Math::VectorUtil::DeltaR(genvec[genvecid]->vector(),recovec[j]->vector());	
	double dptrel = fabs(recovec[j]->pt()-ptgen)/ptgen;
	if (dptrel<mindpT && dR<0.4) {
	  mindpT = dptrel;
	  recotogenmatch[i]=std::pair<unsigned,bool>(2000,false);
	} else if (dR<mindR){
	  recotogenmatch[i]=std::pair<unsigned,bool>(2000,false);
	} else {
	  recotogenmatch[j]=std::pair<unsigned,bool>(2000,false);
	}
	nRecNotMatched++;
	}*/
    }

    if (nReco<=nGen && nRecNotMatched>0) std::cout << " -- Warning: " << __FILE__ << " line " << __LINE__ << ", found " << nRecNotMatched << "/" << nReco << " recovec not matched, genvec size is: " << nGen << std::endl;
    
  }



  /*
  //---------------------------------------------------------------------------
  //Typedefs
  //---------------------------------------------------------------------------
  typedef std::pair<GenParticle*, GenParticle*> GenPGenP;
  typedef std::pair<GenJet*, GenJet*>           GenJGenJ;
  typedef std::pair<GenJet*, GenParticle*>      GenJGenP;
  typedef std::pair<GenParticle*, Electron*>    GenPElec;
  typedef std::pair<Electron*, Electron*>       ElecElec;
  typedef std::pair<GenParticle*, Muon*>        GenPMuon;
  typedef std::pair<Muon*, Muon*>               MuonMuon;
  typedef std::pair<Jet*, GenJet*>              RecJGenJ;

  //---------------------------------------------------------------------------
  //STL-like
  //---------------------------------------------------------------------------
  template <class T, class U> void erase_if(T & t, U pred) {
    t.erase(std::remove_if(t.begin(), t.end(), pred), t.end());
  }

  //---------------------------------------------------------------------------
  //Predicate (may require binding for STL)
  //---------------------------------------------------------------------------
  bool HasPdgId(GenParticle const* part, int pdgid, bool ignoreSign);
    
  
  bool MinPtMaxEtaConst(GenJet const* cand, double const& minPt, double const& maxEta);
  
  bool InEcalGap(Electron const* electron);
  
  bool ElectronIdByWP(Electron const* electron, std::string const& wp);
  
  bool ElectronIsoByWP(Electron const* electron,std::string const& wp);
  
  
  
  bool LeptonDBLessThan(Lepton const* lepton, double const& cut);
  
  bool BTaggedHE(Jet const* jet);
  
  bool BTaggedHP(Jet const* jet);
  
  bool PairMassInRange(Candidate const* cand1, Candidate const* cand2, 
                       double const& mLow, double const& mHigh);
  
  bool GenPairOppCharge(GenParticle const* part1, GenParticle const* part2);
  
  bool LeptonPairOppCharge(Lepton const* lepton1, Lepton const* lepton2);
  
  bool PartonJetFlavour(GenJet const* genJet, unsigned const& flavour);
  
  template<class T> 
    bool MinDRToCollection(Candidate const* cand, 
                      std::vector<T> const& coll, double const& cut) {
    BOOST_FOREACH(T const& ele, coll) {
      if (DR(cand, ele) < cut) return false;
    }
    return true;
  }
  



  //---------------------------------------------------------------------------
  //Comparator (may require binding for STL)
  //---------------------------------------------------------------------------
  inline bool PtComparator(Candidate const* c1, Candidate const* c2) {
    return (c1->Pt() < c2->Pt());
  }

  template<class T, class U>
  bool DRCompare(std::pair<T,U> const& p1, std::pair<T,U> const& p2) {
      double d1 = DR (p1.first,p1.second);
      double d2 = DR (p2.first,p2.second);
      return (d1 < d2);
  }
  

  //---------------------------------------------------------------------------
  //Mutator
  //---------------------------------------------------------------------------
  void JetEnergyUncertainty(Jet * jet, bool shift_up);
  
  void JetEnergyShift(Jet * jet, double unc);
  
  void SmearResolution(Jet * jet, TRandom2 * rand, double jerUnc);
  

  //---------------------------------------------------------------------------
  //Matching
  //---------------------------------------------------------------------------
  double DR(Candidate const* cand1, Candidate const* cand2);
  
  template<class T, class U>
    std::vector< std::pair<T,U> > MatchByDR(std::vector<T> const& c1,
                                              std::vector<U> const& c2,
                                              double const& maxDR,
                                              bool const& uniqueFirst,
                                              bool const& uniqueSecond) {
      std::vector< std::pair<T,U> > pairVec = MakePairs(c1,c2);
      erase_if(pairVec, !boost::bind(DRLessThan<T,U>, _1, maxDR));
      std::sort(pairVec.begin(), pairVec.end(), DRCompare<T,U>);
      if (uniqueFirst && uniqueSecond) {
        std::vector< std::pair<T,U> > uPairVec;
        std::vector<T> fVec;
        std::vector<U> sVec;
        std::pair<T,U> aPair;
        BOOST_FOREACH(aPair, pairVec) {
          bool inFVec = std::count(fVec.begin(),fVec.end(),aPair.first);
          bool inSVec = std::count(sVec.begin(),sVec.end(),aPair.second);
          if (!inFVec && !inSVec) {
            uPairVec.push_back(aPair);
            fVec.push_back(aPair.first);
            sVec.push_back(aPair.second);
          }
        }
        return uPairVec;
      } else if (uniqueFirst) {
        std::vector< std::pair<T,U> > uPairVec;
        std::vector<T> fVec;
        std::pair<T,U> aPair;
        BOOST_FOREACH(aPair, pairVec) {
          bool inFVec = std::count(fVec.begin(),fVec.end(),aPair.first);
          if (!inFVec) {
            uPairVec.push_back(aPair);
            fVec.push_back(aPair.first);
          }
        }
        return uPairVec;
      } else if (uniqueSecond) {
        std::vector< std::pair<T,U> > uPairVec;
        std::vector<U> sVec;
        std::pair<T,U> aPair;
        BOOST_FOREACH(aPair, pairVec) {
          bool inSVec = std::count(sVec.begin(),sVec.end(),aPair.second);
          if (!inSVec) {
            uPairVec.push_back(aPair);
            sVec.push_back(aPair.second);
          }
        }
        return uPairVec;
      }
      return pairVec;
    }


  //---------------------------------------------------------------------------
  //Collection Generators
  //---------------------------------------------------------------------------
  std::vector<GenParticle*> MakeFinalBQuarksCollection(
      std::vector<GenParticle*> const& partVec);
  
  std::vector<GenParticle*> MakeInitialBQuarksCollection(
      std::vector<GenParticle*> const& partVec);

  std::vector<GenParticle*> MakeFinalBHadronsCollection(
      std::vector<GenParticle*> const& partVec);
  
  //Searches a collection of GenParticles for a specific pdgId which, when
  //found, initiates a search over the subsequent particles for any daughters.
  //If found, these are added to the returned vector. 
  std::vector<GenParticle*> MakeDaughterCollection(
      std::vector<GenParticle*> const& partVec, unsigned const& pdgid);
  
  std::vector<GenParticle*> MakeDaughterCollectionMCATNLO(
      std::vector<GenParticle*> const& partVec, unsigned const& pdgid);

  std::vector<GenParticle*> GetStatus1Muons(
      std::vector<GenParticle*> const& partVec);
  
  //---------------------------------------------------------------------------
  //Helpers
  //---------------------------------------------------------------------------
  template <class T>
    std::vector<T*> MakePtrVector(std::vector<T> & vec) {
      std::vector<T*> ptrVec(vec.size(),0);
      for (unsigned i = 0; i < vec.size(); ++i) ptrVec[i] = &(vec[i]);
      return ptrVec;
    }

  template <class T, class U>
    std::vector<T> ExtractFirst(std::vector< std::pair<T,U> > const& pairVec) {
      std::vector<T> vec;
      vec.resize(pairVec.size());
      for (unsigned i = 0; i < vec.size(); ++i) vec[i] = (pairVec[i].first);
      return vec;
    }

  template <class T, class U>
    std::vector<U> ExtractSecond(std::vector< std::pair<T,U> > const& pairVec) {
      std::vector<U> vec;
      vec.resize(pairVec.size());
      for (unsigned i = 0; i < vec.size(); ++i) vec[i] = (pairVec[i].second);
      return vec;
    }

  template <typename T>
    std::vector< std::pair<T,T> > MakePairs(std::vector<T> const& collection) {
      unsigned n = collection.size();
      std::vector< std::pair<T,T> > pairVec;
      if (n == 0) return pairVec;
      pairVec.resize((n*(n-1))/2);
      unsigned vecIndex = 0;
      for (unsigned i = 0; i < (n-1); ++i) {
        for (unsigned j = (i+1); j < n; ++j, ++vecIndex) {
          pairVec[vecIndex] = (std::pair<T,T>(collection[i],collection[j]));
        }
      }
      return pairVec;
    }
  template <class T, class U>
    std::vector< std::pair<T,U> > MakePairs(std::vector<T> const& collection1,
                                            std::vector<U> const& collection2) {
      unsigned n = collection1.size();
      unsigned m = collection2.size();
      std::vector< std::pair<T,U> > pairVec;
      pairVec.resize(n*m);
      unsigned vecIndex = 0;
      for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < m; ++j, ++vecIndex) {
          pairVec[vecIndex] = (std::pair<T,U>(collection1[i],collection2[j]));
        }
      }
      return pairVec;
    }

  bool IsBflavPdgId(int pdgid);


  double ScaleReso(double pT, double eta, double jerUnc);

*/
  
  union ui64 {
    uint64_t one;
    int16_t four[4];
  };
  
  std::set<int16_t> GetTriggerTypes(TriggerObject* obj);
  
} // namepsace
#endif

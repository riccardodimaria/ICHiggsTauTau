#ifndef UserCode_ICHiggsTauTau_ICMuonIDProducer_h
#define UserCode_ICHiggsTauTau_ICMuonIDProducer_h
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
//
// class declaration
//

using namespace std;
using namespace reco;
class ICMuonIDProducer : public edm::EDFilter {
public:
    explicit ICMuonIDProducer(const edm::ParameterSet&);
    ~ICMuonIDProducer();
    
private:
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    
    // ----------member data ---------------------------
    bool verbose_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    
};
#endif

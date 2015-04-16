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
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//
// class declaration
//

using namespace std;
using namespace reco;
class RecoMuonIdProducer : public edm::EDFilter {
public:
    explicit RecoMuonIdProducer(const edm::ParameterSet&);
    ~RecoMuonIdProducer();
    
private:
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    
    // ----------member data ---------------------------
    bool verbose_;
    edm::EDGetTokenT<edm::View<reco::Muon>> muonToken_;
    
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RecoMuonIdProducer::RecoMuonIdProducer(const edm::ParameterSet& iConfig) {
	muonToken_ = consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonTag"));
    
    produces<edm::ValueMap<float> >("");
    
}


RecoMuonIdProducer::~RecoMuonIdProducer()
{
    
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool RecoMuonIdProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	using namespace edm;
    std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>() );
    
	Handle<edm::View<reco::Muon>> muCollection;
	iEvent.getByToken(muonToken_,muCollection);
    const edm::View<reco::Muon> muCandidates = (*muCollection.product());
    
    std::vector<float> values;
    values.reserve(muCollection->size());
    for ( edm::View<reco::Muon>::const_iterator muIter = muCandidates.begin(); muIter != muCandidates.end(); ++muIter) {
        std::cout<<muon::segmentCompatibility(*muIter)<<std::endl;
        
        double muIdglb = 0.;
        double muId = 0.;
       muIdglb = (*muIter).isGlobalMuon()&&(*muIter).globalTrack()->normalizedChi2()<3&&(*muIter).combinedQuality().chi2LocalPosition<12&&(*muIter).combinedQuality().trkKink<20;
       if(((*muIter).innerTrack()).isNonnull()){;
         muId= (*muIter).innerTrack()->validFraction()>=0.8&&muon::segmentCompatibility(*muIter)>=(muIdglb ? 0.303:0.451); 
       }
        values.push_back( muId );
      }
    
    edm::ValueMap<float>::Filler filler(*out);
    filler.insert(muCollection, values.begin(), values.end() );
	filler.fill();
    
	iEvent.put(out);
    
	return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoMuonIdProducer);

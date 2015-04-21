// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "UserCode/ICHiggsTauTau/plugins/ICElectronIDMVAProducerMiniAOD.hh"
ICElectronIDMVAProducerMiniAOD::ICElectronIDMVAProducerMiniAOD(const edm::ParameterSet& iConfig) {
    verbose_ = iConfig.getUntrackedParameter<bool>("verbose", false);
  //  vertexToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexTag"));
	electronToken_ = consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronTag"));
//    ecalRechitEBToken_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
 //   ecalRechitEEToken_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEERecHitCollection"));
	method_ = iConfig.getParameter<string>("method");
	std::vector<string> fpMvaWeightFiles = iConfig.getParameter<std::vector<std::string> >("mvaWeightFile");
	Trig_ = iConfig.getParameter<bool>("Trig");
    
    produces<edm::ValueMap<float> >("");
    
    mvaID_ = new EGammaMvaEleEstimatorCSA14();
    
    EGammaMvaEleEstimatorCSA14::MVAType type_;
    
    if(Trig_){type_ = EGammaMvaEleEstimatorCSA14::kTrig;}
    
    if(!Trig_){
        if (fpMvaWeightFiles.size() > 4) type_ = EGammaMvaEleEstimatorCSA14::kNonTrigPhys14;
        else type_ = EGammaMvaEleEstimatorCSA14::kNonTrig;
    }
    bool manualCat_ = true;
    
	string path_mvaWeightFileEleID;
	for(unsigned ifile=0 ; ifile < fpMvaWeightFiles.size() ; ++ifile) {
        path_mvaWeightFileEleID = edm::FileInPath ( fpMvaWeightFiles[ifile].c_str() ).fullPath();
        mvaWeightFiles_.push_back(path_mvaWeightFileEleID);
	}
    
    mvaID_->initialize(method_, type_, manualCat_, mvaWeightFiles_);
    
}


ICElectronIDMVAProducerMiniAOD::~ICElectronIDMVAProducerMiniAOD()
{
    
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool ICElectronIDMVAProducerMiniAOD::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	using namespace edm;
    std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>() );
    
//	Handle<reco::VertexCollection>  vertexCollection;
//	iEvent.getByToken(vertexToken_, vertexCollection);
    
    
	Handle<pat::ElectronCollection> egCollection;
	iEvent.getByToken(electronToken_,egCollection);
    const pat::ElectronCollection egCandidates = (*egCollection.product());
    
    
    
    std::vector<float> values;
    values.reserve(egCollection->size());
    for ( pat::ElectronCollection::const_iterator egIter = egCandidates.begin(); egIter != egCandidates.end(); ++egIter) {
        
        double mvaVal = -999999;
//        mvaVal = mvaID_->mvaValue( *egIter, *pv,thebuilder,lazyTools, verbose_);
        mvaVal = mvaID_->mvaValue(*egIter,verbose_);
        
        values.push_back( mvaVal );
	}
    
    edm::ValueMap<float>::Filler filler(*out);
    filler.insert(egCollection, values.begin(), values.end() );
	filler.fill();
    
	iEvent.put(out);
    
	return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ICElectronIDMVAProducerMiniAOD);

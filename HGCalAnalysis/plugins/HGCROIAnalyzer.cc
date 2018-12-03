#include "RecoNtuples/HGCalAnalysis/plugins/HGCROIAnalyzer.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "SimG4CMS/Calo/interface/CaloHitID.h"

#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"

#include "DataFormats/GeometryVector/interface/Basic3DVector.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Plane3D.h"
#include "CLHEP/Geometry/Vector3D.h"

#include "TVector2.h"

#include <iostream>

using namespace std;


//
HGCROIAnalyzer::HGCROIAnalyzer( const edm::ParameterSet &iConfig ) : 
  mc_( consumes<edm::HepMCProduct>(edm::InputTag("generatorSmeared")) ),
  recHitsEE_( consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCEERecHits")) ),
  genParticles_( consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles")) ),
  simTracks_( consumes<std::vector<SimTrack>>(edm::InputTag("g4SimHits")) ),
  simVertices_( consumes<std::vector<SimVertex>>(edm::InputTag("g4SimHits")) ),
  genBarCodes_(consumes<std::vector<int>>(edm::InputTag("genParticles")) ),
  slimmedHits_(new std::vector<SlimmedHit>),
  slimmedROIs_(new std::vector<SlimmedROI>),
  genVertex_(new TLorentzVector),
  dEdXWeights_(iConfig.getParameter<std::vector<double>>("dEdXWeights")),
  thicknessCorrection_(iConfig.getParameter<std::vector<double>>("thicknessCorrection")),
  byClosest_(iConfig.getParameter<bool>("byClosest"))
{  
  edm::Service<TFileService> fs;
  tree_=fs->make<TTree>("data","data");
  tree_->Branch("Hits",      "std::vector<SlimmedHit>", &slimmedHits_);
  tree_->Branch("ROIs",      "std::vector<SlimmedROI>", &slimmedROIs_);
  tree_->Branch("GenVertex", "TLorentzVector",          &genVertex_);
}

//
HGCROIAnalyzer::~HGCROIAnalyzer()
{
}

//
void HGCROIAnalyzer::endJob()
{
}

//
GlobalPoint HGCROIAnalyzer::projectHitPositionAt(float z,float eta,float phi){
  float theta=2*TMath::ATan(exp(-eta));
  float rho=z*TMath::Tan(theta);
  GlobalPoint xyz(rho*TMath::Cos(phi),rho*TMath::Sin(phi),z);
  return xyz;
}
  
//
void HGCROIAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{

  recHitTools_.getEventSetup(iSetup);
  edm::ESHandle<HGCalGeometry> geomHandle;
  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",geomHandle);

  //read the generated primary vertex
  edm::Handle<edm::HepMCProduct> mcHandle;
  iEvent.getByToken(mc_, mcHandle);
  HepMC::GenVertex *primaryVertex = *(mcHandle)->GetEvent()->vertices_begin();
  genVertex_->SetXYZT(primaryVertex->position().x() / 10.,
                      primaryVertex->position().y() / 10.,
                      primaryVertex->position().z() / 10.,
                      primaryVertex->position().t());


  //check MC truth on photons (from hard process, converted, unconverted)
  /*
  edm::Handle<std::vector<SimTrack>> simTracksHandle;
  iEvent.getByToken(simTracks_, simTracksHandle);
  edm::Handle<std::vector<SimVertex>> simVerticesHandle;
  iEvent.getByToken(simVertices_, simVerticesHandle);
  PhotonMCTruthFinder mcTruth;
  std::vector<PhotonMCTruth> photons=mcTruth.find(*(simTracksHandle.product()),*(simVerticesHandle.product()));
  for(auto a:photons){
    bool isConverted(a.isAConversion()>0);
  }
  */
  
  //read gen level particles
  slimmedROIs_->clear();
  edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle;
  iEvent.getByToken(genParticles_, genParticlesHandle);
  for(size_t i = 0; i < genParticlesHandle->size(); ++i )  {    
    const reco::GenParticle &p = (*genParticlesHandle)[i];
    if(p.pdgId()!=22) continue;
    if(!p.isPromptFinalState()) continue;    
    if(fabs(p.eta())<1.5 || fabs(p.eta())>2.9) continue;
    SlimmedROI photon(p.pt(),p.eta(),p.phi(),p.mass(),22);
    slimmedROIs_->push_back(photon);
  }
  if(slimmedROIs_->size()==0) return;

  //shift ROIs in phi to create control ROIs for noise/pileup control
  for(size_t ir=0; ir<2; ir++) {
    TLorentzVector p4=slimmedROIs_->at(ir).p4();
    float dphi=2.*TMath::Pi()/6.;
    for(int iphi=1; iphi<=5; iphi++) {
      float newPhi=TVector2::Phi_mpi_pi(p4.Phi()+iphi*dphi);
      SlimmedROI noise(0,p4.Eta(),newPhi,0,-(ir+1));
      slimmedROIs_->push_back(noise);
    }
  }
  
  //collect rec hits in regions of interest
  slimmedHits_->clear();
  edm::Handle<HGCRecHitCollection> recHitsHandleEE;
  iEvent.getByToken(recHitsEE_, recHitsHandleEE);
  int nUnmatched(0);
  for(size_t i=0; i< recHitsHandleEE->size(); i++){
    const HGCRecHit &h = (*recHitsHandleEE)[i];
    const DetId did = h.detid();
    GlobalPoint xyz=recHitTools_.getPosition(did);
    TVector2 xy(xyz.x(),xyz.y());
    float z = xyz.z();
 
    //check if it is within ROI assuming linear propagation from genVertex
    int matchedROIidx(-1), signalRegion(-1);
    for(size_t ir=0; ir<slimmedROIs_->size(); ir++){
      float eta=slimmedROIs_->at(ir).eta();
      float phi=slimmedROIs_->at(ir).phi();      
      GlobalPoint xyzExp=projectHitPositionAt(z-genVertex_->Z(),eta,phi);      
      TVector2 xyExp(xyzExp.x(),xyzExp.y());
      float d=(xyExp-xy).Mod();

      //get closest detid if required
      if(byClosest_) {
        DetId closest_did=geomHandle.product()->getClosestCell(xyzExp);
        GlobalPoint xyzClosest(recHitTools_.getPosition(closest_did));
        TVector2 xyClosest(xyzClosest.x(),xyzClosest.y());        
        d=(xyClosest-xy).Mod();
      }

      if(d>5.3) continue;
      matchedROIidx=ir; 
      signalRegion=3;
      if(d<=2.6) signalRegion=2;
      if(d<=1.3) signalRegion=1;
      break;
    }
    if(matchedROIidx<0) {
      nUnmatched++;
      continue;
    }

    //store information
    unsigned int layer = recHitTools_.getLayerWithOffset(did);
    double thickness   = recHitTools_.getSiThickness(did);
    float en=h.energy();
    double mip = dEdXWeights_[layer] * 0.001;  // convert in GeV
    size_t thickCorrIdx(0);
    if (thickness > 199. && thickness < 201) thickCorrIdx++;
    if (thickness > 299. && thickness < 301) thickCorrIdx++;
    mip /= thicknessCorrection_[thickCorrIdx];
    float en_mips=en/mip;
    if(en_mips<1) continue;

    float t=h.time();

    SlimmedHit sh(layer,matchedROIidx,signalRegion,en,en_mips,t,0);
    slimmedHits_->push_back(sh);
  }
  cout << "Rejected " << nUnmatched << "/" << recHitsHandleEE->size() << " EE hits" << endl;

  tree_->Fill();
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCROIAnalyzer);

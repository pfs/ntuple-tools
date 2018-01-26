#ifndef _HGCOccupancyAnalyzer_h_
#define _HGCOccupancyAnalyzer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimG4CMS/Calo/interface/HGCNumberingScheme.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"

#include "SimCalorimetry/HGCSimProducers/interface/HGCDigitizerBase.h"  

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"

#include <string>

/**
   @class SubdetectorOccupancyHisto
   @author P. Silva (CERN)
*/
class SubdetectorOccupancyHisto
{
public:
  
  /**
     @short CTOR
   */
 SubdetectorOccupancyHisto(int sdcode,edm::Service<TFileService> *fs) : sdcode_(sdcode),fs_(fs), nEvents_(0)
  { 
    thr_.push_back(2);
    thr_.push_back(4);
    thr_.push_back(8);
    thr_.push_back(20);
    thr_.push_back(40);
  }
  
  /**
     @short check if layer exists already
   */
  bool hasLayer(int layer) 
  { 
    return normHistos_.find(layer)!=normHistos_.end();
  }
  
  /**
     @short initiate montoring histos for a layer
   */
  void initLayer(int layer)
  {
    if(hasLayer(layer)) return;
    TString name("sd_"); name += sdcode_; name += "_layer"; name += layer;
    normHistos_[layer] = (*fs_)->make<TH1F>(name+"_nch",";Pseudo-rapidity;Number of channels",10,1.5,3.0);

    mipHistos_[layer] = (*fs_)->make<TH3F>(name+"_mip",";# MIPs;Pseudo-rapidity;Time sample [bx]",200,0,50,10,1.5,3.0,10,0,10);
    mipHistos_[layer]->Sumw2();

    totHistos_[layer] = (*fs_)->make<TH3F>(name+"_tot",";ToT [bx];Pseudo-rapidity;Time sample [bx]",9,1,10,10,1.5,3.0,10,0,10);
    totHistos_[layer]->Sumw2();

    busyHistos_[layer] = (*fs_)->make<TH3F>(name+"_busy",";Busy state;Pseudo-rapidity;Time sample [bx]",2,0,2,10,1.5,3.0,10,0,10);
    busyHistos_[layer]->Sumw2();

    dataCountHistos_[layer] = new TH1F(name+"_datacount",";Pseudo-rapidity;Data volume [bit]",10,1.5,3.0);
    dataCountHistos_[layer]->SetDirectory(0);

    dataVolHistos_[layer] = (*fs_)->make<TH2F>(name+"_datavol",";Data volume [bit];Pseudo-rapidity",60,0,12,10,1.5,3.0);
    dataVolHistos_[layer]->Sumw2();

    dataVolTypeHistos_[layer] = (*fs_)->make<TH2F>(name+"_datavoltype",";Data volume [bit];Pseudo-rapidity",3,0,3,10,1.5,3.0);
    dataVolTypeHistos_[layer]->GetXaxis()->SetBinLabel(1,"=1 bit");
    dataVolTypeHistos_[layer]->GetXaxis()->SetBinLabel(2,"=8 bits");
    dataVolTypeHistos_[layer]->GetXaxis()->SetBinLabel(3,"=12 bits");
    dataVolTypeHistos_[layer]->Sumw2();

    trigCountHistos_[layer] = new TH1F(name+"_trigcount",";Trigger volume [bit];Pseudo-rapidity",10,1.5,3.0);
    trigCountHistos_[layer]->SetDirectory(0);

    trigVolHistos_[layer] = (*fs_)->make<TH2F>(name+"_trigvol",";Trigger volume [bit];Pseudo-rapidity",60,0,6,10,1.5,3.0);
    trigVolHistos_[layer]->Sumw2();

    std::map<int,TH1F *> layerCountHistos,layerCountTHistos;   
    std::map<int,TH2F *> layerOccHistos,layerOccTHistos;

    /* Double_t occbins[100]; */
    /* for(size_t ibin=0; ibin<10; ibin++)   occbins[ibin]=ibin*1e-4; */
    /* for(size_t ibin=10; ibin<19; ibin++)  occbins[ibin]=(ibin-9.)*1e-3; */
    /* for(size_t ibin=19; ibin<29; ibin++)  occbins[ibin]=(ibin-18)*1e-2; */
    /* for(size_t ibin=29; ibin<100; ibin++) occbins[ibin]=0.1+(1-0.1)/(100-29)*(ibin-28); */
    for(size_t ithr=0; ithr<thr_.size(); ithr++)
      {
	TString thrName(name); thrName += "_thr"; thrName += thr_[ithr];
	layerCountHistos[ thr_[ithr] ] = new TH1F(thrName+"_count",";Counts;Pseudo-rapidity",10,1.5,3.0);
	layerCountHistos[ thr_[ithr] ]->SetDirectory(0);
	layerCountTHistos[ thr_[ithr] ] = (TH1F *) layerCountHistos[ thr_[ithr] ]->Clone(thrName+"_countT");
	layerCountTHistos[ thr_[ithr] ]->SetDirectory(0);
	
	layerOccHistos[ thr_[ithr] ] = (*fs_)->make<TH2F>(thrName+"_occ",";Occupancy;Pseudo-rapidity",200,0,1,10,1.5,3.0);
	layerOccHistos[ thr_[ithr] ]->Sumw2();
	layerOccTHistos[ thr_[ithr] ] = (*fs_)->make<TH2F>(thrName+"_occT",";Occupancy;Pseudo-rapidity",200,0,1,10,1.5,3.0);
	layerOccTHistos[ thr_[ithr] ]->Sumw2();
	//layerOccHistos[ thr_[ithr] ] = (*fs_)->make<TH2F>(thrName+"_occ",";Occupancy;Pseudo-rapidity",99,occbins,10,1.5,3.0);

      }
    
    countHistos_[layer]=layerCountHistos;
    countTHistos_[layer]=layerCountTHistos;
    occHistos_[layer]=layerOccHistos;
    occTHistos_[layer]=layerOccTHistos;
  }


  /**
     @short accumulate counts for a new hit
  */
  void count(int layer,float eta,int adc)
  {
    eta=fabs(eta);
    float adcT=adc/TMath::CosH(eta);
    
    int dataVol=computeDataVol(adc,eta);
    dataCountHistos_[layer]->Fill(eta,dataVol);
    if(dataVol==1)  dataVolTypeHistos_[layer]->Fill(0.,eta);
    if(dataVol==8)  dataVolTypeHistos_[layer]->Fill(1.,eta);
    if(dataVol==12) dataVolTypeHistos_[layer]->Fill(2.,eta);

    trigCountHistos_[layer]->Fill(eta,computeTriggerVol(adc,eta));
    for(size_t ithr=0; ithr<thr_.size(); ithr++)
      {
	if(adc>=thr_[ithr])
          countHistos_[layer][ thr_[ithr] ]->Fill(eta);
        if(adcT>=thr_[ithr])
          countTHistos_[layer][ thr_[ithr] ]->Fill(eta);
      }
  }

  /**
     @short computes number of bits to send
   */
  int computeDataVol(float adc,float eta)
  {
    float nmips=0.25*adc;
    int nbits(1);
    if(nmips>6.4)      nbits=12;
    else if(nmips>0.4) nbits=8;
    return nbits;
  }
  int computeTriggerVol(float adc,float eta)
  {
    float nmips=0.25*adc;
    int nbits(1);
    if(fabs(eta)<2)
      {
	if(nmips>96)      nbits=10;
        else if(nmips>10) nbits=6;
      }
    else if(fabs(eta)<2.5)
      {
	if(nmips>96)      nbits=10;
	else if(nmips>10) nbits=6;
      }
    else
      {
	if(nmips>192)     nbits=10;
	else if(nmips>25) nbits=6;
      }
    return nbits;
  }

  /**
     @short to be called at the end of an event
  */
  float endEvent(bool reset=true)
  {
    nEvents_++;

    //iterate over layers
    float totalEvtSize(0);
    for(std::map< int, std::map<int,TH2F *> >::iterator it=occHistos_.begin();
	it!=occHistos_.end();
	it++)
      {

	//occupancies
	TH1F *normH=normHistos_[it->first];

	//iterate over thresholds
	for(std::map<int,TH2F *>::iterator jt=it->second.begin();
	    jt!=it->second.end();
	    jt++)
	  {
	    TH2F *occH=jt->second;
	    TH1F *countH=countHistos_[it->first][jt->first];

            TH2F *occTH=occTHistos_[it->first][jt->first];
	    TH1F *countTH=countTHistos_[it->first][jt->first];

	    //fill occupancy histos
	    for(int xbin=1; xbin<=countH->GetXaxis()->GetNbins(); xbin++)
	      {
		float eta=countH->GetXaxis()->GetBinCenter(xbin);
		float cts=countH->GetBinContent(xbin);
		float norm=normH->GetBinContent(xbin)*2; //we have two endcaps
		float occ(norm>0 ? cts/norm : 0.);
		occH->Fill(occ,eta);
                
		cts=countTH->GetBinContent(xbin);
		occ=(norm>0 ? cts/norm : 0.);
		occTH->Fill(occ,eta);
	      }
	    
	    //reset counts
	    if(reset){
              countH->Reset("ICE");
              countTH->Reset("ICE");
            }
          }
	
	//data volumes
	TH1F *dataH=dataCountHistos_[it->first];
	TH1F *trigH=trigCountHistos_[it->first];
	for(int xbin=1; xbin<=dataH->GetXaxis()->GetNbins(); xbin++)
	  {
	    float ncells=normH->GetBinContent(xbin)*2; //we have two endcaps
	    float eta=dataH->GetXaxis()->GetBinCenter(xbin);

	    float data=dataH->GetBinContent(xbin);
	    totalEvtSize += data;
	    dataVolHistos_[it->first]->Fill(ncells>0 ? data/ncells : 1. ,eta);

	    float trig=trigH->GetBinContent(xbin);
	    trigVolHistos_[it->first]->Fill(ncells>0 ? trig/ncells : 1. ,eta);
	  }
	
	dataH->Reset("ICE");
	trigH->Reset("ICE");
      }

    return  totalEvtSize;
  }

  /**
     @short normalize according to the number of events analyzed
   */
  int finalize()
  {
    if(nEvents_==0) return 0;

    for(std::map< int, TH3F *>::iterator it=mipHistos_.begin();
	it!=mipHistos_.end();
	it++)
      {

	//divide by the number of cells analyzed
	TH1F *normH=normHistos_[it->first];
	for(int ybin=1; ybin<=it->second->GetYaxis()->GetNbins(); ybin++)
	  {
	    float ncells=normH->GetBinContent(ybin)*2; //we have two endcaps
	    if(ncells==0) continue;

	    //energy spectrum 
	    for(int xbin=0; xbin<=it->second->GetXaxis()->GetNbins()+1; xbin++)
	      {
		for(int zbin=0; zbin<=it->second->GetZaxis()->GetNbins()+1; zbin++)
		  {
		    float counts=it->second->GetBinContent(xbin,ybin,zbin);
		    float countsErr=it->second->GetBinError(xbin,ybin,zbin);
		    it->second->SetBinContent(xbin,ybin,counts/(nEvents_*ncells));
		    it->second->SetBinError(xbin,ybin,countsErr/(nEvents_*ncells));
		  }
	      }

	    //data type sent
	    for(int xbin=0; xbin<=dataVolTypeHistos_[it->first]->GetXaxis()->GetNbins(); xbin++)
	      {
		float counts=dataVolTypeHistos_[it->first]->GetBinContent(xbin,ybin);
		float countsErr=dataVolTypeHistos_[it->first]->GetBinError(xbin,ybin);
		dataVolTypeHistos_[it->first]->SetBinContent(xbin,ybin,counts/(nEvents_*ncells));
		dataVolTypeHistos_[it->first]->SetBinError(xbin,ybin,countsErr/(nEvents_*ncells));
	      }
	  }



	//normalized by the number of events (already averaged per cells)
	for(std::map<int,TH2F *>::iterator jt=occHistos_[it->first].begin();
	    jt!=occHistos_[it->first].end();
	    jt++)
          {
            jt->second->Scale(1./nEvents_);
            occTHistos_[it->first][jt->first]->Scale(1./nEvents_);
          }

	dataVolHistos_[it->first]->Scale(1./nEvents_);
	trigVolHistos_[it->first]->Scale(1./nEvents_);

	mipHistos_[it->first]->Scale(1./nEvents_);
	totHistos_[it->first]->Scale(1./nEvents_);
	busyHistos_[it->first]->Scale(1./nEvents_);    
      }
    return nEvents_;
  }

  /**
     @short DTOR
   */
  ~SubdetectorOccupancyHisto() {}

  //all histos are public and can be manipulated...
  std::map< int, TH1F *>                normHistos_;
  std::map< int, TH1F *>                dataCountHistos_, trigCountHistos_;
  std::map< int, TH2F *>                dataVolHistos_,   trigVolHistos_, dataVolTypeHistos_;
  std::map< int, std::map<int,TH1F *> > countHistos_,countTHistos_;
  std::map< int, std::map<int,TH2F *> > occHistos_,occTHistos_;
  std::map< int, TH3F *>                mipHistos_, totHistos_, busyHistos_;
 
 private:
  
  int sdcode_;
  edm::Service<TFileService> *fs_;
  std::vector<int> thr_;
  int nEvents_;
};


/**
   @class HGCOccupancyAnalyzer
   @author P. Silva (CERN)
*/

class HGCOccupancyAnalyzer : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCOccupancyAnalyzer( const edm::ParameterSet& );
  ~HGCOccupancyAnalyzer();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob();

 private:
  /**
   */
  void prepareAnalysis(std::map<int,const HGCalGeometry *> &hgcGeometries);
  
  /**
     @short digi analyzers
   */
  float analyzeHEDigis(size_t isd,edm::Handle<HGCHEDigiCollection> &heDigis,const HGCalGeometry *geom);
  float analyzeEEDigis(size_t isd,edm::Handle<HGCEEDigiCollection> &eeDigis,const HGCalGeometry *geom);
  
  //gen level
  bool isInit_;

  //hgcal
  std::vector<std::string> geometrySource_;
  std::vector<std::string> digiCollections_;
  std::vector<SubdetectorOccupancyHisto> occHistos_;

  //event size 
  TH2F *evtSizeH_;
};
 

#endif

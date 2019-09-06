// system include files
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"


#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "TH1D.h"
#include "TH2D.h"

class HGCalScintLayerClusterStudy : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources> {

public:

  explicit HGCalScintLayerClusterStudy(const edm::ParameterSet&);
  ~HGCalScintLayerClusterStudy() override {}
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual void beginJob() override {}
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;

  // --------------member data ---------------
  std::string           nameDetector_;
  edm::EDGetToken       recHitSource_;
  edm::EDGetToken       layerClusterSource_;
  int                   verbosity_;
  unsigned int          layers_;
  int                   firstLayer_;

  std::vector<TH2D*>    EtaPhi_Plus_;
  std::vector<TH2D*>    EtaPhi_Minus_;
  std::vector<TH2D*>    XY_Plus_;
  std::vector<TH2D*>    XY_Minus_;
};


HGCalScintLayerClusterStudy::HGCalScintLayerClusterStudy(const edm::ParameterSet& iConfig) :
  nameDetector_(iConfig.getParameter<std::string>("DetectorName")), 
  verbosity_(iConfig.getUntrackedParameter<int>("Verbosity",0)),
  layers_(0), firstLayer_(1) {

  usesResource(TFileService::kSharedResource);

  auto rechit        = iConfig.getParameter<edm::InputTag>("RecHitSource");
  auto layercluster  = iConfig.getParameter<edm::InputTag>("LayerClusterSource");
  
  if (nameDetector_ == "HGCalEESensitive" || 
      nameDetector_ == "HGCalHESiliconSensitive" ||
      nameDetector_ == "HGCalHEScintillatorSensitive" ) {
    recHitSource_    = consumes<HGCRecHitCollection>(rechit);
    layerClusterSource_ = consumes<std::vector<reco::CaloCluster> >(layercluster);
  } else {
    throw cms::Exception("BadSource")
      << "HGCal DetectorName given as " << nameDetector_ << " must be: "
      << "\"HGCalHESiliconSensitive\", \"HGCalHESiliconSensitive\", or \"HGCalHEScintillatorSensitive\"!"; 
  }
  edm::LogVerbatim("HGCalStudies") << "Initialize HGCalScintLayerClusterStudy for " << nameDetector_;
}

void HGCalScintLayerClusterStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("DetectorName","HGCalHEScintillatorSensitive");
  desc.add<edm::InputTag>("RecHitSource",edm::InputTag("HGCalRecHit","HGCHEBRecHits"));
  desc.add<edm::InputTag>("LayerClusterSource",edm::InputTag("hgcalLayerClusters",""));
  desc.addUntracked<int>("Verbosity",0);
  descriptions.add("hgcalLayerClusterStudyHESintillator",desc);
}

void HGCalScintLayerClusterStudy::beginRun(edm::Run const&,
                                edm::EventSetup const& iSetup) {

  edm::ESHandle<HGCalDDDConstants>  pHGDC;
  iSetup.get<IdealGeometryRecord>().get(nameDetector_, pHGDC);
  const HGCalDDDConstants & hgcons_ = (*pHGDC);
  layers_     = hgcons_.layers(true);
  firstLayer_ = hgcons_.firstLayer();

  if (false) {
    std::cout << "layers: " << layers_ << " | firstLayer: " << firstLayer_ << std::endl;
    std::cout << "--- " << hgcons_.cellSizeTrap(1, 1).first << " | " << hgcons_.cellSizeTrap(1, 1).second << std::endl;
    for (unsigned int l=0; l < layers_; l++) {
      std::cout <<  hgcons_.getREtaRange(l).first << " | " << hgcons_.getREtaRange(l).second << std::endl;
    }
  }
  
  edm::Service<TFileService> fs;
  char histoname[100];
  for (unsigned int il = 0; il < firstLayer_+layers_; il++) {
    sprintf (histoname, "XY_Plus_Layer_%d", il);
    XY_Plus_.push_back(fs->make<TH2D>(histoname, "Occupancy", 250, -250., 250., 250, -250., 250.));
    sprintf (histoname, "XY_Minus_Layer_%d", il);
    XY_Minus_.push_back(fs->make<TH2D>(histoname, "Occupancy", 250, -250., 250., 250, -250., 250.));
    sprintf (histoname, "EtaPhi_Plus_Layer_%d", il);
    EtaPhi_Plus_.push_back(fs->make<TH2D>(histoname, "Occupancy", 74., 1.42, 3.0, 288., -CLHEP::pi, CLHEP::pi));
    sprintf (histoname, "EtaPhi_Minus_Layer_%d", il);
    EtaPhi_Minus_.push_back(fs->make<TH2D>(histoname, "Occupancy", 74., -3.0, -1.42, 288., -CLHEP::pi, CLHEP::pi));
  }
}

void HGCalScintLayerClusterStudy::analyze(const edm::Event& iEvent, 
			       const edm::EventSetup& iSetup) {

  edm::Handle<HGCRecHitCollection> theRecHits;
  iEvent.getByToken(recHitSource_, theRecHits);

  edm::Handle<std::vector<reco::CaloCluster> > theLayerClusters;
  iEvent.getByToken(layerClusterSource_, theLayerClusters);

  edm::ESHandle<HGCalGeometry> geom;
  iSetup.get<IdealGeometryRecord>().get(nameDetector_, geom);
  if (!geom.isValid()) edm::LogWarning("HGCalStudies") << "Cannot get valid HGCalGeometry Object for " << nameDetector_;
  const HGCalGeometry* geom0 = geom.product();

  int ncluster(0);
  for (const auto & cluster : *(theLayerClusters.product())) {
    ncluster+=1;
    const std::vector< std::pair<DetId, float> > & hits = cluster.hitsAndFractions();
    //std::cout << hits.size() << std::endl;

    for (const auto & hit : hits) {

      if (hit.first.det()==DetId::HGCalHSc) {
	DetId detId = hit.first;
      
	GlobalPoint global = geom0->getPosition(detId);
	int layer   = ((detId.det() == DetId::Forward) ? HGCalDetId(detId).layer() :
		       ((detId.det() == DetId::HGCalHSc) ? 
			HGCScintillatorDetId(detId).layer() : 
			HGCSiliconDetId(detId).layer()));
      
	float globalx   = global.x();
	float globaly   = global.y();
	float globalz   = global.z();
	float globaleta = global.eta();
	float globalphi = global.phi();

	(globalz>0) ? XY_Plus_.at(layer)->Fill(globalx, globaly, ncluster) : XY_Minus_.at(layer)->Fill(globalx, globaly, ncluster);
	(globalz>0) ? EtaPhi_Plus_.at(layer)->Fill(globaleta, globalphi, ncluster) : EtaPhi_Minus_.at(layer)->Fill(globaleta, globalphi, ncluster);
      }
    }
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

#include "FWCore/Framework/interface/MakerMacros.h"

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalScintLayerClusterStudy);

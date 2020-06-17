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
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

class HGCalHitNtuple: public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {

public:
  explicit HGCalHitNtuple(const edm::ParameterSet&);
  ~HGCalHitNtuple() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  virtual void beginJob() override {}
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;

private:
  
  void fill_rec_tree_(int event, const HGCRecHitCollection& hits);

  template<class T>
  void fill_sim_tree_(int event, std::vector<PCaloHit> const& hits, const T* geom, std::string sensorType);

  void fill_cp_tree_(int event, std::vector<CaloParticle> const& cps, std::map<DetId, const HGCRecHit*> hitmap);

  void fill_lc_tree_(int event, std::vector<reco::CaloCluster> const& lcs, std::map<DetId, const HGCRecHit*> hitmap);

  //void fill_lc_ticl_tree_(int event, std::vector<reco::CaloCluster> const& lcs, std::vector<reco::HGCalMultiCluster> const& mcs);

  void computeThreshold();

  void computeRegionBoundary(std::vector<CaloParticle> const& cps);

  struct HitsInfo {
    HitsInfo() {
      event = 0;
      x = y = z = eta = phi = 0.;
      time = energy = 0;
      layer = additional1 = additional2 = 0;
    }
    int event;
    float x, y, z, eta, phi;
    float energy, time;
    int layer, additional1, additional2;
    // additional1 = thickness, size, and cpidx for rechits, lcs, and hitsInSimClusters respectively
    // additional2 = scidx for hitsInSimClusters
  };

  struct CPInfo {
    CPInfo() {
      event = idx = id = 0;
      energy = pt = eta = phi = energy_rec = 0.;
    }
    int event, idx, id;
    float energy, pt, eta, phi, energy_rec;
  };

//  struct LCInfo {
//    LCInfo() {
//      event = 0;
//      x = y = z = eta = phi = 0.;
//      time = energy = 0.;
//      size = layer = 0;
//    }
//    int event;
//    float x, y, z, eta, phi;
//    double time, energy;
//    int size, layer;
//  };

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::PCaloHitContainer>  simHitSourceEE_;
  edm::EDGetTokenT<edm::PCaloHitContainer>  simHitSourceHESi_;
  edm::EDGetTokenT<edm::PCaloHitContainer>  simHitSourceHEScint_;
  edm::EDGetToken recHitSourceEE_;
  edm::EDGetToken recHitSourceHESi_;
  edm::EDGetToken recHitSourceHEScint_;
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticleSource_;
  edm::EDGetTokenT<std::vector<reco::HGCalMultiCluster>> multiClusterSource_;
  edm::EDGetToken layerClusterSource_;
  
  // ----------threshold calculation-----------
  double ecut_;
  std::vector<double> dEdXweights_;
  std::vector<double> thicknessCorrection_;
  std::vector<double> fcPerMip_;
  double fcPerEle_;
  std::vector<double> nonAgedNoises_;
  double noiseMip_;
  int verbosity_;
  
  HitsInfo rechitsInfo, rechitsSignalRegionInfo, simhitsInfo, lcInfo, simhitsInSimClusterInfo;
  CPInfo cpInfo;
  //LCInfo lcInfo, lcNoEMInfo;

  TTree* RecHitTree;
  TTree* RecHitSignalRegionTree;
  TTree* SimHitTree;
  TTree* LCTree;
  TTree* SimHitsInSimClusterTree;
  TTree* CPTree;
  //TTree* LCNoEMTree;

  hgcal::RecHitTools rhtools_;

  unsigned int maxlayer_;

  std::vector<std::vector<double>> thresholds_;
  std::vector<std::vector<double>> v_sigmaNoise_;

  double regionEtaMin_, regionEtaMax_, regionPhiMin_, regionPhiMax_;
};

HGCalHitNtuple::HGCalHitNtuple(const edm::ParameterSet& iConfig) :
  ecut_(iConfig.getParameter<double>("ecut")),
  dEdXweights_(iConfig.getParameter<std::vector<double>>("dEdXweights")),
  thicknessCorrection_(iConfig.getParameter<std::vector<double>>("thicknessCorrection")),
  fcPerMip_(iConfig.getParameter<std::vector<double>>("fcPerMip")),
  fcPerEle_(iConfig.getParameter<double>("fcPerEle")),
  nonAgedNoises_(iConfig.getParameter<edm::ParameterSet>("noises").getParameter<std::vector<double>>("values")),
  noiseMip_(iConfig.getParameter<edm::ParameterSet>("noiseMip").getParameter<double>("noise_MIP")),
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)) {
  usesResource(TFileService::kSharedResource);
  
  simHitSourceEE_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("simSourceEE"));
  simHitSourceHESi_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("simSourceHESi"));
  simHitSourceHEScint_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("simSourceHEScint"));
  recHitSourceEE_ = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("sourceEE"));
  recHitSourceHESi_ = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("sourceHESi"));
  recHitSourceHEScint_ = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("sourceHEScint"));

  caloParticleSource_ = consumes<std::vector<CaloParticle>>(iConfig.getParameter<edm::InputTag>("sourceCaloParticle"));

  //layerClusterSource_ = consumes<std::vector<reco::CaloCluster>>(iConfig.getParameter<edm::InputTag>("sourceLayerCluster"));

  //multiClusterSource_ = consumes<std::vector<reco::HGCalMultiCluster>>(iConfig.getParameter<edm::InputTag>("sourceMultiCluster"));
}

void HGCalHitNtuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<double>("ecut", 3.0);
  desc.add<std::vector<double>>("dEdXweights", {});
  desc.add<std::vector<double>>("thicknessCorrection", {});
  desc.add<std::vector<double>>("fcPerMip", {});
  desc.add<double>("fcPerEle", 0.00016020506);
  edm::ParameterSetDescription descNestedNoises;
  descNestedNoises.add<std::vector<double>>("values", {});
  desc.add<edm::ParameterSetDescription>("noises", descNestedNoises);
  edm::ParameterSetDescription descNestedNoiseMIP;
  descNestedNoiseMIP.add<bool>("scaleByDose", false);
  descNestedNoiseMIP.add<unsigned int>("scaleByDoseAlgo", 0);
  desc.add<edm::ParameterSetDescription>("scaleByDose", descNestedNoiseMIP);
  descNestedNoiseMIP.add<std::string>("doseMap", "");
  desc.add<edm::ParameterSetDescription>("doseMap", descNestedNoiseMIP);
  descNestedNoiseMIP.add<double>("noise_MIP", 1. / 100.);
  desc.add<edm::ParameterSetDescription>("noiseMip", descNestedNoiseMIP);
  
  desc.add<edm::InputTag>("simSourceEE", edm::InputTag("g4SimHits", "HGCHitsEE"));
  desc.add<edm::InputTag>("simSourceHESi", edm::InputTag("g4SimHits", "HGCHitsHEfront"));
  desc.add<edm::InputTag>("simSourceHEScint", edm::InputTag("g4SimHits", "HGCHitsHEback"));
  desc.add<edm::InputTag>("sourceEE", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
  desc.add<edm::InputTag>("sourceHESi", edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
  desc.add<edm::InputTag>("sourceHEScint", edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
  desc.add<edm::InputTag>("sourceCaloParticle", edm::InputTag("mix", "MergedCaloTruth"));
  //desc.add<edm::InputTag>("sourceLayerCluster",edm::InputTag("hgcalLayerClusters",""));
  //desc.add<edm::InputTag>("sourceMultiCluster",edm::InputTag("multiClustersFromTrackstersEM", "MultiClustersFromTracksterByCA"));
  desc.addUntracked<int>("verbosity", 0);
  descriptions.add("hgcalHitNtuple", desc);
}

void HGCalHitNtuple::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {		
  edm::Service<TFileService> fs;
  RecHitTree = fs->make<TTree>("RecHitTree", "");
  RecHitSignalRegionTree = fs->make<TTree>("RecHitSignalRegionTree", "");
  SimHitTree = fs->make<TTree>("SimHitTree", "");
  CPTree = fs->make<TTree>("CPTree", "");
  SimHitsInSimClusterTree = fs->make<TTree>("SimHitsInSimClusterTree", "");
  LCTree = fs->make<TTree>("LCTree", "");
  //LCNoEMTree = fs->make<TTree>("LCNoEMTree", "");
  
  RecHitTree->Branch("RecHitsInfo", &rechitsInfo, "event/I:x/F:y/F:z/F:eta/F:phi/F:energy/F:time/F:layer/I:thickness/I:dummy/I");
  RecHitSignalRegionTree->Branch("RecHitsSignalRegionInfo", &rechitsSignalRegionInfo, "event/I:x/F:y/F:z/F:eta/F:phi/F:energy/F:time/F:layer/I:thickness/I:dummy/I");
  SimHitTree->Branch("SimHitsInfo", &simhitsInfo, "event/I:x/F:y/F:z/F:eta/F:phi/F:energy/F:time/F:layer/I:additional/I:dummy/I");
  LCTree->Branch("lcInfo", &lcInfo, "event/I:x/F:y/F:z/F:eta/F:phi/F:energy/F:time/F:layer/I:size/I:dummy/I");
  SimHitsInSimClusterTree->Branch("simhitsInSimClusterInfo", &simhitsInSimClusterInfo, "event/I:x/F:y/F:z/F:eta/F:phi/F:energy/F:time/F:layer/I:cpidx/I:scidx/I");
  CPTree->Branch("cpInfo", &cpInfo, "event/I:idx/I:id/I:energy/F:pt/F:eta/F:phi/F:energy_rec/F");
  //LCNoEMTree->Branch("lcNoEMInfo", &lcNoEMInfo, "event/I:x/F:y/F:z/F:eta/F:phi/F:energy/D:size/I:layer/I");
}

void HGCalHitNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  rhtools_.getEventSetup(iSetup);
  maxlayer_ = rhtools_.lastLayerBH();

  computeThreshold();

  int Event = iEvent.id().event();
  
  edm::Handle<edm::PCaloHitContainer> handleTheSimHitsEE;
  iEvent.getByToken(simHitSourceEE_, handleTheSimHitsEE);
  
  edm::Handle<edm::PCaloHitContainer> handleTheSimHitsHESi;
  iEvent.getByToken(simHitSourceHESi_, handleTheSimHitsHESi);
  
  edm::Handle<edm::PCaloHitContainer> handleTheSimHitsHEScint;
  iEvent.getByToken(simHitSourceHEScint_, handleTheSimHitsHEScint);
  
  edm::Handle<HGCRecHitCollection> handleTheRecHitsEE;
  iEvent.getByToken(recHitSourceEE_, handleTheRecHitsEE);

  edm::Handle<HGCRecHitCollection> handleTheRecHitsHESi;
  iEvent.getByToken(recHitSourceHESi_, handleTheRecHitsHESi);

  edm::Handle<HGCRecHitCollection> handleTheRecHitsHEScint;
  iEvent.getByToken(recHitSourceHEScint_, handleTheRecHitsHEScint);

  const auto& rechitsEE = *handleTheRecHitsEE;
  const auto& rechitsFH = *handleTheRecHitsHESi;
  const auto& rechitsBH = *handleTheRecHitsHEScint;
  std::map<DetId, const HGCRecHit*> hitmap;
  for (unsigned int i = 0; i < rechitsEE.size(); ++i) {
    hitmap[rechitsEE[i].detid()] = &rechitsEE[i];
  }
  for (unsigned int i = 0; i < rechitsFH.size(); ++i) {
    hitmap[rechitsFH[i].detid()] = &rechitsFH[i];
  }
  for (unsigned int i = 0; i < rechitsBH.size(); ++i) {
    hitmap[rechitsBH[i].detid()] = &rechitsBH[i];
  }
  
  edm::Handle<std::vector<CaloParticle>> handleTheCaloParticle;
  iEvent.getByToken(caloParticleSource_, handleTheCaloParticle);

  //edm::Handle<std::vector<reco::CaloCluster> > handleTheLayerClusters;
  //iEvent.getByToken(layerClusterSource_, handleTheLayerClusters);

  //edm::Handle<std::vector<reco::HGCalMultiCluster> > handleTheMultiClusters;
  //iEvent.getByToken(multiClusterSource_, handleTheMultiClusters);
  
  edm::ESHandle<HGCalGeometry> geomEE;
  edm::ESHandle<HGCalGeometry> geomHESi;
  edm::ESHandle<HGCalGeometry> geomHEScint;
  
  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive", geomEE);
  iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive", geomHESi);
  iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive", geomHEScint);
  
  const HGCalGeometry* geom0EE = geomEE.product();
  const HGCalGeometry* geom0HESi = geomHESi.product();
  const HGCalGeometry* geom0HEScint = geomHEScint.product();

  fill_sim_tree_(Event, *handleTheSimHitsEE, geom0EE, "Hexagon");
  fill_sim_tree_(Event, *handleTheSimHitsHESi, geom0HESi, "Hexagon");
  fill_sim_tree_(Event, *handleTheSimHitsHEScint, geom0HEScint, "Trapezoid");

  computeRegionBoundary(*handleTheCaloParticle);
  //std::cout << "Debug2: " << regionEtaMin_ << " | " << regionEtaMax_ << " | " << regionPhiMin_ << " | " << regionPhiMax_ << "\n";
  //fill_rec_tree_(Event, *handleTheRecHitsEE);
  //fill_rec_tree_(Event, *handleTheRecHitsHESi);
  //fill_rec_tree_(Event, *handleTheRecHitsHEScint);

  fill_cp_tree_(Event, *handleTheCaloParticle, hitmap);

  //fill_lc_tree_(Event, *handleTheLayerClusters, hitmap);

  //fill_lc_ticl_tree_(Event, *handleTheLayerClusters, *handleTheMultiClusters);
}

void HGCalHitNtuple::fill_rec_tree_(int event, const HGCRecHitCollection& hits) {
  
  for (unsigned int i = 0; i < hits.size(); ++i) {

    const HGCRecHit& hit = hits[i];
    DetId detId = hit.detid();
    
    int ithickness = rhtools_.getSiThickIndex(detId);
    int ilayer = rhtools_.getLayerWithOffset(detId)-1;
    
    if (ithickness == -1)
        ithickness = 3;
    double storedThreshold = thresholds_[ilayer][ithickness];
    if (hit.energy() < storedThreshold) continue;
    
    rechitsInfo.event = event;
    rechitsInfo.energy = hit.energy();
    rechitsInfo.time = hit.time();
    //GlobalPoint global = geom->getPosition(detId);
    GlobalPoint global = rhtools_.getPosition(detId);
    rechitsInfo.x = global.x();
    rechitsInfo.y = global.y();
    rechitsInfo.z = global.z();
    rechitsInfo.eta = global.eta();
    rechitsInfo.phi = global.phi();
    rechitsInfo.layer = ilayer+1;
    rechitsInfo.additional1 = ithickness;
    RecHitTree->Fill();
    if (regionEtaMin_ < global.eta() &&
	global.eta() < regionEtaMax_ &&
	regionPhiMin_ < global.phi() &&
	global.phi() < regionPhiMax_) {
      rechitsSignalRegionInfo.event = event;
      rechitsSignalRegionInfo.energy = hit.energy();
      rechitsSignalRegionInfo.time = hit.time();
      //GlobalPoint global = geom->getPosition(detId);
      //GlobalPoint global = rhtools_.getPosition(detId);
      rechitsSignalRegionInfo.x = global.x();
      rechitsSignalRegionInfo.y = global.y();
      rechitsSignalRegionInfo.z = global.z();
      rechitsSignalRegionInfo.eta = global.eta();
      rechitsSignalRegionInfo.phi = global.phi();
      rechitsSignalRegionInfo.layer = ilayer+1;
      rechitsSignalRegionInfo.additional1 = ithickness;
      RecHitSignalRegionTree->Fill();
    }
  }
}

template<class T>
void HGCalHitNtuple::fill_sim_tree_(int event, const std::vector<PCaloHit>& hits, const T* geom, std::string sensorType) {
  for (unsigned int i = 0; i < hits.size(); ++i) {
    simhitsInfo.event = event;
    const PCaloHit& hit = hits[i];
    simhitsInfo.energy = hit.energy();
    simhitsInfo.time = hit.time();
    uint32_t id_ = hits[i].id();
    if(sensorType=="Hexagon"){
      HGCSiliconDetId detId = HGCSiliconDetId(id_);
      GlobalPoint global = geom->getPosition(detId);
      simhitsInfo.x = global.x();
      simhitsInfo.y = global.y();
      simhitsInfo.z = global.z();
      simhitsInfo.eta = global.eta();
      simhitsInfo.phi = global.phi();
    } else {
      HGCScintillatorDetId detId = HGCScintillatorDetId(id_);
      GlobalPoint global = geom->getPosition(detId);
      simhitsInfo.x = global.x();
      simhitsInfo.y = global.y();
      simhitsInfo.z = global.z();
      simhitsInfo.eta = global.eta();
      simhitsInfo.phi = global.phi();
    }
    SimHitTree->Fill();
  }
}

void HGCalHitNtuple::fill_cp_tree_(int event, const std::vector<CaloParticle>& cps, std::map<DetId, const HGCRecHit*> hitmap) {
  for (unsigned int i = 0; i < cps.size(); ++i) {
    const CaloParticle& cp = cps[i];
    if (abs(cp.eta())>1.5 and abs(cp.eta())<3.0) {
      cpInfo.event = event;
      cpInfo.idx = i;
      cpInfo.id = cp.pdgId();
      cpInfo.energy = cp.energy();
      cpInfo.pt = cp.pt();
      cpInfo.eta = cp.eta();
      cpInfo.phi = cp.phi();
      float energy_rec = 0.;
      for (auto const sc : cp.simClusters()) {
	for (auto const h_and_f : sc->hits_and_fractions()) {
	  if (hitmap.count(h_and_f.first)) {
	    energy_rec += hitmap[h_and_f.first]->energy() * h_and_f.second;
	  }
	}
      }
      cpInfo.energy_rec = energy_rec;
      CPTree->Fill();
      if (cp.pt()>1) {
	int scidx(0);
	for (auto const sc : cp.simClusters()) {
	  for (auto const hit: sc->hits_and_fractions()) {
	    DetId detId = hit.first;
	    int ilayer = rhtools_.getLayerWithOffset(detId)-1;
	    simhitsInSimClusterInfo.event = event;
	    simhitsInSimClusterInfo.layer = ilayer+1;
	    simhitsInSimClusterInfo.additional1 = i;
	    simhitsInSimClusterInfo.additional2 = scidx;
	    if (hitmap.count(detId)) {
	      simhitsInSimClusterInfo.energy = hitmap[detId]->energy() * hit.second;
	      simhitsInSimClusterInfo.time = hitmap[detId]->time();
	      GlobalPoint global = rhtools_.getPosition(detId);
	      //GlobalPoint global = geom->getPosition(detId);
	      simhitsInSimClusterInfo.x = global.x();
	      simhitsInSimClusterInfo.y = global.y();
	      simhitsInSimClusterInfo.z = global.z();
	      simhitsInSimClusterInfo.eta = global.eta();
	      simhitsInSimClusterInfo.phi = global.phi();	    
	      SimHitsInSimClusterTree->Fill();
	    }
	  }
	  scidx+=1;
	}
      }
    }
  }
}

void HGCalHitNtuple::fill_lc_tree_(int event, std::vector<reco::CaloCluster> const& lcs, std::map<DetId, const HGCRecHit*> hitmap) {
  for (unsigned int i = 0; i < lcs.size(); ++i) {
    const reco::CaloCluster& lc = lcs[i];
    lcInfo.event = event;
    lcInfo.x = lc.x();
    lcInfo.y = lc.y();
    lcInfo.z = lc.z();
    lcInfo.eta = lc.eta();
    lcInfo.phi = lc.phi();
    lcInfo.energy = lc.energy();
    const std::vector< std::pair<DetId, float> > & hits = lc.hitsAndFractions();
    DetId detId = hits[0].first;
    lcInfo.layer = rhtools_.getLayerWithOffset(detId);
    lcInfo.additional1 = hits.size();
    float time_tot = 0;
    int hits_tot = 0;
    for (auto const hit : hits) {
      if(hitmap[hit.first]->time()>0){
	time_tot += hitmap[hit.first]->time()/hitmap[hit.first]->energy();
	//time_tot += hitmap[hit.first]->time();
	hits_tot += 1;
      }
    }
    if (hits_tot!=0)
      lcInfo.time = time_tot/hits_tot;
    //lcInfo.time = lc.time();
    LCTree->Fill();
  }
}

//void HGCalHitNtuple::fill_lc_ticl_tree_(int event, std::vector<reco::CaloCluster> const& lcs, std::vector<reco::HGCalMultiCluster> const& mcs) {
//  for (unsigned int i = 0; i < lcs.size(); ++i) {
//    const reco::CaloCluster& lc = lcs[i];
//    uint32_t seed = lc.seed().rawId();
//    bool mask = false;
//    for (const auto& mc : mcs) {
//      for (const auto& lc_in_mc : mc) {
//	if (seed == lc_in_mc->seed().rawId()) {
//	  mask = true;
//	  break;
//	}
//      }
//    }
//    if (mask == false) {
//      lcNoEMInfo.event = event;
//      lcNoEMInfo.x = lc.x();
//      lcNoEMInfo.y = lc.y();
//      lcNoEMInfo.z = lc.z();
//      lcNoEMInfo.eta = lc.eta();
//      lcNoEMInfo.phi = lc.phi();
//      lcNoEMInfo.energy = lc.energy();
//      const std::vector< std::pair<DetId, float> > & hits = lc.hitsAndFractions();
//      lcNoEMInfo.size = hits.size();
//      DetId detId = hits[0].first;
//      lcNoEMInfo.layer = rhtools_.getLayerWithOffset(detId);
//      LCNoEMTree->Fill();
//    }
//  }
//}

void HGCalHitNtuple::computeThreshold() {

  std::vector<double> thickness;
  const unsigned maxNumberOfThickIndices = 3;
  thickness.resize(maxNumberOfThickIndices + 1, 0);  // +1 to accomodate for the Scintillators
  thresholds_.resize(maxlayer_, thickness);
  v_sigmaNoise_.resize(maxlayer_, thickness);

  for (unsigned ilayer = 1; ilayer <= maxlayer_; ++ilayer) {
    for (unsigned ithick = 0; ithick < maxNumberOfThickIndices; ++ithick) {
      float sigmaNoise = 0.001f * fcPerEle_ * nonAgedNoises_[ithick] * dEdXweights_[ilayer] /
                         (fcPerMip_[ithick] * thicknessCorrection_[ithick]);
      thresholds_[ilayer - 1][ithick] = sigmaNoise * ecut_;
      v_sigmaNoise_[ilayer - 1][ithick] = sigmaNoise;
    }
    float scintillators_sigmaNoise = 0.001f * noiseMip_ * dEdXweights_[ilayer];
    thresholds_[ilayer - 1][maxNumberOfThickIndices] = ecut_ * scintillators_sigmaNoise;
    v_sigmaNoise_[ilayer - 1][maxNumberOfThickIndices] = scintillators_sigmaNoise;
  }
}

void HGCalHitNtuple::computeRegionBoundary(std::vector<CaloParticle> const& cps) {
  const CaloParticle& cp = cps[0];
  //std::cout << cp.eta() << " | " << cp.phi() << " | " << cp.energy() << "\n";
  regionEtaMin_ = cp.eta()-0.3;
  regionEtaMax_ = cp.eta()+0.3;
  regionPhiMin_ = cp.phi()-0.3;
  regionPhiMax_ = cp.phi()+0.3;
  //std::cout << "Debug1: " << regionEtaMin_ << " | " << regionEtaMax_ << " | " << regionPhiMin_ << " | " << regionPhiMax_ << "\n";
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(HGCalHitNtuple);

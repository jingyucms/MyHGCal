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

#include "DataFormats/Math/interface/deltaPhi.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

class HGCalClusterNtuple: public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {

public:
  explicit HGCalClusterNtuple(const edm::ParameterSet&);
  ~HGCalClusterNtuple() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  virtual void beginJob() override {}
  virtual void endJob() override {}
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override {}
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;

private:
  
  void fill_rechits_tree_(int event, const HGCRecHitCollection& hits);

  void fill_rechits_in_cluster_tree_(int event, std::vector<reco::CaloCluster> const& cluster, std::map<DetId, const HGCRecHit*> hitmap);

  void fill_cluster_tree_(int event, std::vector<reco::CaloCluster> const& cls);

  void fill_hists_(std::vector<CaloParticle> const& cps, std::vector<reco::CaloCluster> const& cls, std::map<DetId, const HGCRecHit*> hitmap);

  void computeThreshold();
  void computeBoundary(std::vector<CaloParticle> const& cps);

  struct HitsInfo {
    HitsInfo() {
      event = 0;
      x = y = z = eta = phi = 0.;
      time = energy = 0;
      layer = additional1 = 0;
      additional2 = -1;
    }
    int event;
    float x, y, z, eta, phi;
    float energy, time;
    int layer, additional1;
    int additional2;
  };

  // ----------member data ---------------------------
  edm::EDGetToken recHitSourceEE_;
  edm::EDGetToken recHitSourceHESi_;
  edm::EDGetToken recHitSourceHEScint_;
  edm::EDGetToken hdbClusterSource_;
  
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticleSource_;
  
  // ----------threshold calculation-----------
  double ecut_;
  std::vector<double> dEdXweights_;
  std::vector<double> thicknessCorrection_;
  std::vector<double> fcPerMip_;
  double fcPerEle_;
  std::vector<double> nonAgedNoises_;
  double noiseMip_;
  bool is200PU_;
  int verbosity_;
  
  HitsInfo rechitsInfo_, rechitsInClusterInfo_, hdbclusterInfo_;

  TTree* RecHitsTree;
  TTree* RecHitsInClusterTree;
  TTree* HdbClusterTree;

  TH1F* h_reconstructable;
  TH1F* h_reconstructed;
  TH1F* h_ratio;
  
  hgcal::RecHitTools rhtools_;

  unsigned int maxlayer_;

  std::vector<std::vector<double>> thresholds_;
  std::vector<std::vector<double>> v_sigmaNoise_;

  float minEta_;
  float maxEta_;
  float minPhi_;
  float maxPhi_;
};

HGCalClusterNtuple::HGCalClusterNtuple(const edm::ParameterSet& iConfig) :
  ecut_(iConfig.getParameter<double>("ecut")),
  dEdXweights_(iConfig.getParameter<std::vector<double>>("dEdXweights")),
  thicknessCorrection_(iConfig.getParameter<std::vector<double>>("thicknessCorrection")),
  fcPerMip_(iConfig.getParameter<std::vector<double>>("fcPerMip")),
  fcPerEle_(iConfig.getParameter<double>("fcPerEle")),
  nonAgedNoises_(iConfig.getParameter<edm::ParameterSet>("noises").getParameter<std::vector<double>>("values")),
  noiseMip_(iConfig.getParameter<edm::ParameterSet>("noiseMip").getParameter<double>("noise_MIP")),
  is200PU_(iConfig.getUntrackedParameter<bool>("is200PU", true)),
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)) {
  usesResource(TFileService::kSharedResource);
  
  recHitSourceEE_ = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("sourceEE"));
  recHitSourceHESi_ = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("sourceHESi"));
  recHitSourceHEScint_ = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("sourceHEScint"));

  hdbClusterSource_ = consumes<std::vector<reco::CaloCluster>>(iConfig.getParameter<edm::InputTag>("sourceHDBCluster"));

  caloParticleSource_ = consumes<std::vector<CaloParticle>>(iConfig.getParameter<edm::InputTag>("sourceCaloParticle"));
}

void HGCalClusterNtuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
  descNestedNoiseMIP.add<std::string>("doseMap", "");
  descNestedNoiseMIP.add<double>("noise_MIP", 1. / 100.);
  desc.add<edm::ParameterSetDescription>("noiseMip", descNestedNoiseMIP);
  desc.addUntracked<bool>("is200PU", true);
  
  desc.add<edm::InputTag>("sourceEE", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
  desc.add<edm::InputTag>("sourceHESi", edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
  desc.add<edm::InputTag>("sourceHEScint", edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
  
  desc.add<edm::InputTag>("sourceHDBCluster",edm::InputTag("hgcalHDBClusters","HDBClusters", "RECO"));
  desc.add<edm::InputTag>("sourceCaloParticle", edm::InputTag("mix", "MergedCaloTruth"));

  desc.addUntracked<int>("verbosity", 0);
  descriptions.add("hgcalClusterNtuple", desc);
}

void HGCalClusterNtuple::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {

  //std::cout << "beginRun\n";
  edm::Service<TFileService> fs;
  RecHitsTree = fs->make<TTree>("RecHitsTree", "");
  RecHitsInClusterTree = fs->make<TTree>("RecHitsInClusterTree", "");
  HdbClusterTree = fs->make<TTree>("HdbClusterTree", "");

  h_reconstructable = fs->make<TH1F>("h_reconstructable", "", 1000, 0, 1000);
  h_reconstructed = fs->make<TH1F>("h_reconstructed", "", 1000, 0, 1000);
  h_ratio = fs->make<TH1F>("h_ratio", "", 100, 0, 5);
  
  RecHitsTree->Branch("RecHitsInfo", &rechitsInfo_, "event/I:x/F:y/F:z/F:eta/F:phi/F:energy/F:time/F:layer/I:thickness/I:dummy/I");
  RecHitsInClusterTree->Branch("RecHitsInClusterInfo", &rechitsInClusterInfo_, "event/I:x/F:y/F:z/F:eta/F:phi/F:energy/F:time/F:layer/I:thickness/I:label/I");
  HdbClusterTree->Branch("HdbClusterInfo", &hdbclusterInfo_, "event/I:x/F:y/F:z/F:eta/F:phi/F:energy/F:time/F:dummy1/I:dummy2/I:dummy3/I");
}

void HGCalClusterNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //std::cout << "analyze\n";
  
  rhtools_.getEventSetup(iSetup);
  maxlayer_ = rhtools_.lastLayerBH();

  computeThreshold();

  int Event = iEvent.id().event();
  
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
  
  edm::Handle<std::vector<reco::CaloCluster> > handleTheHDBClusters;
  iEvent.getByToken(hdbClusterSource_, handleTheHDBClusters);

  edm::Handle<std::vector<CaloParticle>> handleTheCaloParticle;
  iEvent.getByToken(caloParticleSource_, handleTheCaloParticle);

  computeBoundary(*handleTheCaloParticle);

  fill_rechits_tree_(Event, *handleTheRecHitsEE);
  fill_rechits_tree_(Event, *handleTheRecHitsHESi);
  fill_rechits_tree_(Event, *handleTheRecHitsHEScint);

  fill_rechits_in_cluster_tree_(Event, *handleTheHDBClusters, hitmap);

  fill_cluster_tree_(Event, *handleTheHDBClusters);

  fill_hists_(*handleTheCaloParticle, *handleTheHDBClusters, hitmap);
}

void HGCalClusterNtuple::fill_rechits_tree_(int event, const HGCRecHitCollection& hits) {
  //std::cout << "fill_rechits_tree_\n";
  for (unsigned int i = 0; i < hits.size(); ++i) {

    const HGCRecHit& hit = hits[i];
    DetId detId = hit.detid();
    
    int ithickness = rhtools_.getSiThickIndex(detId);
    int ilayer = rhtools_.getLayerWithOffset(detId)-1;
    
    if (ithickness == -1)
        ithickness = 3;
    double storedThreshold = thresholds_[ilayer][ithickness];
    if (hit.energy() < storedThreshold) continue;
   
    rechitsInfo_.event = event;
    rechitsInfo_.energy = hit.energy();
    rechitsInfo_.time = hit.time();
    //GlobalPoint global = geom->getPosition(detId);
    GlobalPoint global = rhtools_.getPosition(detId);
    if (is200PU_ and (global.eta() > maxEta_ or global.eta() < minEta_ or global.phi() > maxPhi_ or global.phi() < minPhi_)) continue;
    rechitsInfo_.x = global.x();
    rechitsInfo_.y = global.y();
    rechitsInfo_.z = global.z();
    rechitsInfo_.eta = global.eta();
    rechitsInfo_.phi = global.phi();
    rechitsInfo_.layer = ilayer+1;
    rechitsInfo_.additional1 = ithickness;
    RecHitsTree->Fill();
  }
}

void HGCalClusterNtuple::fill_rechits_in_cluster_tree_(int event, std::vector<reco::CaloCluster> const& clusters, std::map<DetId, const HGCRecHit*> hitmap) {
  //std::cout << "fill_rechits_in_cluster_tree_\n";
  for (unsigned int i = 0; i < clusters.size(); ++i) {
    const reco::CaloCluster& cluster = clusters[i];
    for (auto const h_and_f : cluster.hitsAndFractions()) {
      DetId detId = h_and_f.first;
      
      int ithickness = rhtools_.getSiThickIndex(detId);
      int ilayer = rhtools_.getLayerWithOffset(detId)-1;
      
      rechitsInClusterInfo_.event = event;
      rechitsInClusterInfo_.layer = ilayer+1;
      rechitsInClusterInfo_.additional1 = ithickness;
      rechitsInClusterInfo_.additional2 = i;
      if (hitmap.count(detId)) {
	rechitsInClusterInfo_.energy = hitmap[detId]->energy() * h_and_f.second;
	rechitsInClusterInfo_.time = hitmap[detId]->time();
	GlobalPoint global = rhtools_.getPosition(detId);
	if (is200PU_ and (global.eta() > maxEta_ or global.eta() < minEta_ or global.phi() > maxPhi_ or global.phi() < minPhi_)) continue;
	rechitsInClusterInfo_.x = global.x();
	rechitsInClusterInfo_.y = global.y();
	rechitsInClusterInfo_.z = global.z();
	rechitsInClusterInfo_.eta = global.eta();
	rechitsInClusterInfo_.phi = global.phi();
	RecHitsInClusterTree->Fill();
      }
    }
  }
}

void HGCalClusterNtuple::fill_cluster_tree_(int event, std::vector<reco::CaloCluster> const& cls) {
  for (unsigned int i = 0; i < cls.size(); ++i) {
    const reco::CaloCluster& cluster = cls[i];
    hdbclusterInfo_.event = event;
    hdbclusterInfo_.energy = cluster.energy();
    hdbclusterInfo_.time = 0;
    hdbclusterInfo_.x = cluster.x();
    hdbclusterInfo_.y = cluster.y();
    hdbclusterInfo_.z = cluster.z();
    hdbclusterInfo_.eta = cluster.eta();
    hdbclusterInfo_.phi = cluster.phi();
    hdbclusterInfo_.layer = 0;
    hdbclusterInfo_.additional1 = 0;
    hdbclusterInfo_.additional2 = 0;
    HdbClusterTree->Fill();
  }
}

void HGCalClusterNtuple::fill_hists_(const std::vector<CaloParticle>& cps,
				     std::vector<reco::CaloCluster> const& cls,
				     std::map<DetId, const HGCRecHit*> hitmap){

  int size = cps.size();
  if (size > 2) size = 1;

  std::vector<int> matched;
  matched.resize(0);
  for (int i = 0; i < size; ++i) {
    const CaloParticle& cp = cps[i];
  
    float energy_rec = 0.;
    for (auto const sc : cp.simClusters()) {
      for (auto const h_and_f : sc->hits_and_fractions()) {
	if (hitmap.count(h_and_f.first)) {
	  energy_rec += hitmap[h_and_f.first]->energy() * h_and_f.second;
	}
      }
    }
    h_reconstructable->Fill(energy_rec);
    
    int size = cls.size();
    float energy, distance;
    distance = 0.1;
    energy = 0.;
    int imatched = -1;
    for (int j = 1; j < size; ++j) {
      const reco::CaloCluster& cluster = cls[j];
      float dist = 9999;
      for (auto const h_and_f : cluster.hitsAndFractions()) {
	DetId detId = h_and_f.first;
	GlobalPoint global = rhtools_.getPosition(detId);

	//int ilayer = rhtools_.getLayerWithOffset(detId)-1;
	//if (ilayer > 28) continue;
	
	float deta = cp.eta() - global.eta();
	float dphi = reco::deltaPhi(cp.phi(), global.phi());
	float d = sqrt(deta*deta + dphi*dphi);
	//std::cout << dist << "\n";
	if (d < dist) dist = d;
      }
      if (dist < distance and std::find(matched.begin(), matched.end(), j) == matched.end()) {
	imatched = j;
	energy = cluster.energy();
	distance = dist;
      }
    }
    matched.push_back(imatched);
    h_reconstructed->Fill(energy);
    std::cout << matched << " " << energy << "\n";
    
    float ratio = energy/energy_rec;
    
    h_ratio->Fill(ratio);
    
    if (ratio < 0.25) std::cout << "------\n"; 
  }
}

void HGCalClusterNtuple::computeThreshold() {
  //std::cout << "computeThreshold\n";
  std::vector<double> dummy;
  const unsigned maxNumberOfThickIndices = 3;
  dummy.resize(maxNumberOfThickIndices + 1, 0);  // +1 to accomodate for the Scintillators
  thresholds_.resize(maxlayer_, dummy);
  v_sigmaNoise_.resize(maxlayer_, dummy);

  //std::cout << "maxlayer_: " << maxlayer_ << "\n";
  //std::cout << "maxNumberOfThickIndices: " << maxNumberOfThickIndices << "\n";

  for (unsigned ilayer = 1; ilayer <= maxlayer_; ++ilayer) {
    for (unsigned ithick = 0; ithick < maxNumberOfThickIndices; ++ithick) {
      float sigmaNoise = 0.001f * fcPerEle_ * nonAgedNoises_[ithick] * dEdXweights_[ilayer] / (fcPerMip_[ithick] * thicknessCorrection_[ithick]);
      //std::cout << "ilayer: " << ilayer << "\n";
      //std::cout << "ithick: " << ithick << "\n";
      //std::cout << "sigmaNoise: " << sigmaNoise << "\n";
      thresholds_[ilayer - 1][ithick] = sigmaNoise * ecut_;
      v_sigmaNoise_[ilayer - 1][ithick] = sigmaNoise;
      //std::cout << "ilayer: " << ilayer << " nonAgedNoises: " << nonAgedNoises_[ithick]
      //	<< " fcPerEle: " << fcPerEle_ << " fcPerMip: " << fcPerMip_[ithick]
      //	<< " noiseMip: " << fcPerEle_ * nonAgedNoises_[ithick] / fcPerMip_[ithick]
      //	<< " sigmaNoise: " << sigmaNoise << "\n";
    }
    float scintillators_sigmaNoise = 0.001f * noiseMip_ * dEdXweights_[ilayer];
    thresholds_[ilayer - 1][maxNumberOfThickIndices] = ecut_ * scintillators_sigmaNoise;
    v_sigmaNoise_[ilayer - 1][maxNumberOfThickIndices] = scintillators_sigmaNoise;
    //std::cout << "ilayer: " << ilayer << " noiseMip: " << noiseMip_
    //<< " scintillators_sigmaNoise: " << scintillators_sigmaNoise << "\n";
  }
}

void HGCalClusterNtuple::computeBoundary(const std::vector<CaloParticle>& cps){
  const CaloParticle& cp = cps[0];
  float cpEta = cp.eta();
  float cpPhi = cp.phi();
  minEta_ = cpEta - 0.2;
  maxEta_ = cpEta + 0.2;
  minPhi_ = cpPhi - 0.2;
  maxPhi_ = cpPhi + 0.2;
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(HGCalClusterNtuple);

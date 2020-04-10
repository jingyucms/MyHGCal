#ifndef __MyHGCal_HDBClustering_HGCalHDBClusterProducer_H__
#define __MyHGCal_HDBClustering_HGCalHDBClusterProducer_H__

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "MyHGCal/HDBClustering/interface/HGCalHDBClusteringAlgoFactory.h"

using Density = hgcal_hdb_clustering::Density;

class HGCalHDBClusterProducer : public edm::stream::EDProducer<> {
public:
  HGCalHDBClusterProducer(const edm::ParameterSet&);
  ~HGCalHDBClusterProducer() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<HGCRecHitCollection> hits_ee_token;
  edm::EDGetTokenT<HGCRecHitCollection> hits_fh_token;
  edm::EDGetTokenT<HGCRecHitCollection> hits_bh_token;

  //reco::CaloCluster::AlgoId algoId;

  std::unique_ptr<HGCalHDBClusteringAlgoBase> algo;
};

HGCalHDBClusterProducer::HGCalHDBClusterProducer(const edm::ParameterSet& ps) {
  hits_ee_token = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCEEInput"));
  hits_fh_token = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCFHInput"));
  hits_bh_token = consumes<HGCRecHitCollection>(ps.getParameter<edm::InputTag>("HGCBHInput"));

  auto pluginPSet = ps.getParameter<edm::ParameterSet>("plugin");
  algo = HGCalHDBClusteringAlgoFactory::get()->create(pluginPSet.getParameter<std::string>("type"), pluginPSet);
  
  produces<std::vector<reco::BasicCluster>>("HDBClusters");
  //density
  produces<Density>();
  //time for layer clusters
}

void HGCalHDBClusterProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // hgcalHDBClusters
  edm::ParameterSetDescription desc;
  edm::ParameterSetDescription pluginDesc;
  pluginDesc.addNode(edm::PluginDescription<HGCalHDBClusteringAlgoFactory>("type", "HDBv1", true));

  desc.add<edm::ParameterSetDescription>("plugin", pluginDesc);
  desc.add<edm::InputTag>("HGCEEInput", edm::InputTag("HGCalRecHit", "HGCEERecHits"));
  desc.add<edm::InputTag>("HGCFHInput", edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
  desc.add<edm::InputTag>("HGCBHInput", edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
  descriptions.add("hgcalHDBClusters", desc);
}

void HGCalHDBClusterProducer::produce(edm::Event& evt, const edm::EventSetup& es) {

  edm::Handle<HGCRecHitCollection> ee_hits;
  edm::Handle<HGCRecHitCollection> fh_hits;
  edm::Handle<HGCRecHitCollection> bh_hits;

  std::unique_ptr<std::vector<reco::BasicCluster>> clusters(new std::vector<reco::BasicCluster>);
  auto density = std::make_unique<Density>();

  algo->reset();

  algo->getEventSetup(es);

  //make a map detid-rechit
  // NB for the moment just host EE and FH hits
  // timing in digi for BH not implemented for now

  evt.getByToken(hits_ee_token, ee_hits);
  algo->populate(*ee_hits);

  evt.getByToken(hits_fh_token, fh_hits);
  algo->populate(*fh_hits);

  evt.getByToken(hits_bh_token, bh_hits);
  algo->populate(*bh_hits);

  algo->makeClusters();
  *clusters = algo->getClusters();
  evt.put(std::move(clusters), "HDBClusters");

  //Keep the density
  *density = algo->getDensity();
  evt.put(std::move(density));
}

DEFINE_FWK_MODULE(HGCalHDBClusterProducer);

#endif

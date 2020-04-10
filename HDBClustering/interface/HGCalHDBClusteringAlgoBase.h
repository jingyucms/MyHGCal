#ifndef MyHGCal_HDBClustering_HGCalHDBClusteringAlgoBase_h
#define MyHGCal_HDBClustering_HGCalHDBClusteringAlgoBase_h

#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

// C/C++ headers
#include <vector>
#include <numeric>

namespace hgcal_hdb_clustering {
  template <typename T>
  std::vector<size_t> sorted_indices(const std::vector<T> &v) {
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    std::iota(std::begin(idx), std::end(idx), 0);

    // sort indices based on comparing values in v
    std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });

    return idx;
  }

  template <typename T>
  size_t max_index(const std::vector<T> &v) {
    // initialize original index locations
    std::vector<size_t> idx(v.size(), 0);
    std::iota(std::begin(idx), std::end(idx), 0);

    // take the max index based on comparing values in v
    auto maxidx = std::max_element(
        idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1].data.rho < v[i2].data.rho; });

    return (*maxidx);
  }

  //Density collection
  typedef std::map<DetId, float> Density;

};  // namespace hgcal_db_clustering

class HGCalHDBClusteringAlgoBase {
public:
  enum VerbosityLevel { pDEBUG = 0, pWARNING = 1, pINFO = 2, pERROR = 3 };

  HGCalHDBClusteringAlgoBase(VerbosityLevel v) : verbosity_(v) {};
  virtual ~HGCalHDBClusteringAlgoBase() {};

  virtual void populate(const HGCRecHitCollection &hits) = 0;
  virtual void makeClusters() = 0;
  virtual std::vector<reco::BasicCluster> getClusters() = 0;
  virtual void reset() = 0;
  virtual hgcal_hdb_clustering::Density getDensity() = 0;
  virtual void getEventSetupPerAlgorithm(const edm::EventSetup &es) {}

  inline void getEventSetup(const edm::EventSetup &es) {
    rhtools_.getEventSetup(es);
    maxlayer_ = rhtools_.lastLayer(false);
    lastLayerEE_ = rhtools_.lastLayerEE(false);
    lastLayerFH_ = rhtools_.lastLayerFH();
    firstLayerBH_ = rhtools_.firstLayerBH();
    scintMaxIphi_ = rhtools_.getScintMaxIphi();
    getEventSetupPerAlgorithm(es);
  }
  inline void setVerbosity(VerbosityLevel the_verbosity) { verbosity_ = the_verbosity; }

  //max number of layers
  unsigned int maxlayer_;
  // last layer per subdetector
  unsigned int lastLayerEE_;
  unsigned int lastLayerFH_;
  unsigned int firstLayerBH_;
  int scintMaxIphi_;
  bool isNose_;

protected:
  // The verbosity level
  VerbosityLevel verbosity_;

  // The vector of clusters
  std::vector<reco::BasicCluster> clusters_v_;

  hgcal::RecHitTools rhtools_;

};

#endif

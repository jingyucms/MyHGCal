#ifndef MyHGCal_HDBClustering_HGCalHDBAlgo_h
#define MyHGCal_HDBClustering_HGCalHDBAlgo_h

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/Math/interface/deltaPhi.h"

#include "MyHGCal/HDBClustering/interface/HGCalHDBClusteringAlgoBase.h"
#include "MyHGCal/HDBClustering/interface/HGCalAngularTiles.h"
#include "MyHGCal/HDBClustering/interface/HGCalDisjointSets.h"

using Density=hgcal_hdb_clustering::Density;

#include <algorithm>

template <typename T>
class HGCalHDBAlgoT : public HGCalHDBClusteringAlgoBase {
 public :
  
 HGCalHDBAlgoT(const edm::ParameterSet& ps)
   : HGCalHDBClusteringAlgoBase((HGCalHDBClusteringAlgoBase::VerbosityLevel)ps.getUntrackedParameter<unsigned int>("verbosity", 1)),
     NSamples_(ps.getParameter<int>("NSamples")),
     NCluster_(ps.getParameter<int>("NCluster")),
     Algo_(ps.getParameter<unsigned int>("Algo")),
     delta_(ps.getParameter<double>("delta")),
     ecut_(ps.getParameter<double>("ecut")),
     dEdXweights_(ps.getParameter<std::vector<double>>("dEdXweights")),
     thicknessCorrection_(ps.getParameter<std::vector<double>>("thicknessCorrection")),
     fcPerMip_(ps.getParameter<std::vector<double>>("fcPerMip")),
     fcPerEle_(ps.getParameter<double>("fcPerEle")),
     nonAgedNoises_(ps.getParameter<edm::ParameterSet>("noises").getParameter<std::vector<double>>("values")),
     noiseMip_(ps.getParameter<edm::ParameterSet>("noiseMip").getParameter<double>("noise_MIP")) {}

  ~HGCalHDBAlgoT() override {}

  static void fillPSetDescription(edm::ParameterSetDescription& iDesc){
    iDesc.add<int>("NSamples", 10);
    iDesc.add<int>("NCluster", 40);
    iDesc.add<unsigned int>("Algo", 1);
    iDesc.add<double>("delta", 0.04);
    iDesc.add<double>("ecut", 3.0);
    iDesc.addUntracked<unsigned int>("verbosity", 0);
    iDesc.add<std::vector<double>>("dEdXweights", {});
    iDesc.add<std::vector<double>>("thicknessCorrection", {});
    iDesc.add<std::vector<double>>("fcPerMip", {});
    iDesc.add<double>("fcPerEle", 0.0);
    edm::ParameterSetDescription descNestedNoises;
    descNestedNoises.add<std::vector<double>>("values", {});
    iDesc.add<edm::ParameterSetDescription>("noises", descNestedNoises);
    edm::ParameterSetDescription descNestedNoiseMIP;
    descNestedNoiseMIP.add<bool>("scaleByDose", false);
    descNestedNoiseMIP.add<unsigned int>("scaleByDoseAlgo", 0);
    descNestedNoiseMIP.add<std::string>("doseMap", "");
    descNestedNoiseMIP.add<double>("noise_MIP", 1. / 100.);
    iDesc.add<edm::ParameterSetDescription>("noiseMip", descNestedNoiseMIP);
  }

  void populate(const HGCRecHitCollection &hits) override;

  void makeClusters() override;

  std::vector<reco::BasicCluster> getClusters() override;

  void reset() override {
    clusters_v_.clear();
    cells_.clear();
    min_span_tree_.clear();
    single_linkage_tree_.clear();
    condensed_tree_.clear();
  }
  
  Density getDensity() override;

 private:

  // The two parameters used to identify clusters
  int NSamples_;
  int NCluster_;

  // Which algo to use
  // 0: leaf clustering
  // 1: eom clustering
  unsigned int Algo_;
  
  // Search Radius
  double delta_;

  // The hit energy cutoff
  double ecut_;

  // For keeping the density per hit
  Density density_;

  // various parameters used for calculating the noise levels for a given sensor (and whether to use
  // them)
  std::vector<double> dEdXweights_;
  std::vector<double> thicknessCorrection_;
  std::vector<double> fcPerMip_;
  double fcPerEle_;
  std::vector<double> nonAgedNoises_;
  double noiseMip_;

  std::vector<std::vector<double>> thresholds_;
  std::vector<std::vector<double>> v_sigmaNoise_;
  void computeThreshold();

  HGCalDisjointSets set_;

  int nCells_;

  float maxLambda_ = 500;
  float minLambda_ = 0;

  struct Cells {
    std::vector<DetId> detid;
    std::vector<int> layer;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;
    std::vector<float> energy;
    std::vector<float> rho;
    std::vector<float> sigmaNoise;
    std::vector<int> label;
    
    void clear() {
      detid.clear();
      layer.clear();
      eta.clear();
      phi.clear();
      x.clear();
      y.clear();
      z.clear();
      energy.clear();
      rho.clear();
      sigmaNoise.clear();
      label.clear();
    }
  };

  Cells cells_;

  struct Graph {

    int src;
    int dest;
    float weight;

    bool operator<(const Graph& g) const { return weight < g.weight; }
  };

  std::vector<Graph> safeEdges_;
  std::vector<Graph> min_span_tree_;

  struct Tree {
    std::vector<int> parent;
    std::vector<int> left;
    std::vector<int> right;
    std::vector<float> weight;
    std::vector<int> size;
    void clear() {
      parent.clear();
      left.clear();
      right.clear();
      weight.clear();
      size.clear();
    }
  };

  Tree single_linkage_tree_;
  
  struct Hierarchy {
    std::vector<int> parent;
    std::vector<int> child;
    std::vector<float> lambda;
    std::vector<int> size;
    std::vector<float> mass;
    void clear() {
      parent.clear();
      child.clear();
      lambda.clear();
      size.clear();
      mass.clear();
    }
  };
  
  Hierarchy condensed_tree_;

  std::vector<int> nodes_in_slt_;
  std::vector<int> nodes_in_hierarchy_;
  
  inline int findPreSplitInSLT(int node, int start_idx) {
    if (verbosity_ > 2)
      std::cout << "findPreSplitInSLT: " << node << "\n";
    if (std::find(nodes_in_slt_.rbegin(), nodes_in_slt_.rend(), node) != nodes_in_slt_.rend()) return node;
    for (int i = start_idx; i<nCells_; i++) {
      if (single_linkage_tree_.left[i] != node and single_linkage_tree_.right[i] != node)
	continue;
      int new_node = single_linkage_tree_.parent[i];
      if (std::find(nodes_in_slt_.rbegin(), nodes_in_slt_.rend(), new_node) != nodes_in_slt_.rend()) {
	single_linkage_tree_.parent[start_idx] = new_node; // path compression
	return new_node;
      } else {
	node = new_node;
	start_idx = i;
	break;
      }
    }
    return findPreSplitInSLT(node, start_idx);
  }

  inline int getNodeSizeInSLT(int node, int start) {
    int node_index(0);
    for (int i = start - 1; i >= 0; i--) 
      if (node == single_linkage_tree_.parent[i]) {
	node_index = i;
	break;
      }
    return single_linkage_tree_.size[node_index];
  }

  inline int findPreSplitInHierarchy(int node) {
    auto it = std::find(nodes_in_slt_.rbegin(), nodes_in_slt_.rend(), node);
    int node_index = std::distance(begin(nodes_in_slt_), it.base())-1;
    return nodes_in_hierarchy_[node_index];
  }

  int hierarchySize_;
  
  std::vector<int> leaf_clusters_;

  inline int findClusterIdxInLeafClusters(int node) {
    std::vector<int>::iterator it = std::find(leaf_clusters_.begin(), leaf_clusters_.end(), node);
    if (it == leaf_clusters_.end()) return -1;
    return std::distance(leaf_clusters_.begin(), it);
  }
  
  std::vector<std::vector<int>> eom_clusters_;
  
  std::vector<int> findChildren(int node) {
    int rootNode = nCells_;
    std::vector<int> children;
    std::queue<int> q;
    q.emplace(node);
    int protection = 0;
    while (!q.empty() and children.size() < 5) {
      protection ++ ; 
      int inode = q.front();
      for (int i = rootNode; i < hierarchySize_+1; i++) {
  	if (inode == condensed_tree_.parent[i] and i != hierarchySize_) {
  	  q.emplace(condensed_tree_.child[i]);
  	  q.emplace(condensed_tree_.child[i]+1);
  	  children.emplace_back(condensed_tree_.child[i]);
  	  children.emplace_back(condensed_tree_.child[i]+1);
  	  q.pop();
  	  break;
  	} else if (i == hierarchySize_) {
  	  q.pop();
  	  break;
  	}
      }
    }
    return children;
  }
  
  int calcNChildren(int node, std::vector<int> children) {
    int size = 0;
    for (int child: children) {
      int idx = child - 1;
      size += condensed_tree_.size[idx];
    }
    return size;
  }
  
  float calcMassChildren(int node, std::vector<int> children) {
    float childrenMass = 0;
    for (int child: children) {
      int idx = child - 1;
      childrenMass += condensed_tree_.mass[idx];
    }
    return childrenMass;
  }

  inline int findClusterIdxInEOMClusters(int node) {
    int nClusters = eom_clusters_.size();
    for (int i = 0; i < nClusters; i++) {
      auto cluster = eom_clusters_[i];
      for (auto cl : cluster) {
	if (cl == node) return i;
      }
    }
    return -1;
  }
  
  void prepareCellDataStructures();
  void calculateDensity(const T& t, int NSamples, float delta);
  void findSafeEdges(const T& t, float delta);
  void buildMST(const T& t, float delta);
  void prepareSLTDataStructures();
  void buildSLT();
  void prepareClusterHierarchyStructures();
  void buildClusterHierarchy(int NCluster);
  void findLeafClusters();
  void findEOMClusters();
  void assignClusters();
  math::XYZPoint calculatePosition(const std::vector<int>& v) const;

  inline float distance2(int cell1, int cell2) const {  // distance squared
    float d2 = 9999;
    if (abs(cells_.layer[cell1] - cells_.layer[cell2]) < 10 and
	((cells_.layer[cell1] < 50 and cells_.layer[cell2] < 50) or
	 (cells_.layer[cell1] >= 50 and cells_.layer[cell2] >= 50))) {
      const float dphi = reco::deltaPhi(cells_.phi[cell1], cells_.phi[cell2]);
      const float deta = cells_.eta[cell1] - cells_.eta[cell2];
      d2 = deta * deta + dphi * dphi;
    } 
    return d2;
  }

  inline float distance(int cell1, int cell2) const {  // 2-d distance on the layer (x-y)
    float d2 = distance2(cell1, cell2);
    float d = 9999;
    if (d2 < 9000) {
      d = std::sqrt(d2);
      if (cells_.layer[cell1] > 77 and cells_.layer[cell1] < 90) d = d/2;
      if (cells_.layer[cell1] > 89) d = d/4;
    }
    return d;
  }  
};

extern template class HGCalHDBAlgoT<HGCalAngularTiles>;

using HGCalHDBAlgo = HGCalHDBAlgoT<HGCalAngularTiles>;

#endif

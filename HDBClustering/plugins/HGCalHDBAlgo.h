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

template <typename T, typename S>
class HGCalHDBAlgoT : public HGCalHDBClusteringAlgoBase {
 public :
  
 HGCalHDBAlgoT(const edm::ParameterSet& ps)
   : HGCalHDBClusteringAlgoBase((HGCalHDBClusteringAlgoBase::VerbosityLevel)ps.getUntrackedParameter<unsigned int>("verbosity", 3)),
     NSamples_(ps.getParameter<unsigned int>("NSamples")),
     NCluster_(ps.getParameter<unsigned int>("NCluster")),
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
    iDesc.add<unsigned int>("NSamples", 20);
    iDesc.add<unsigned int>("NCluster", 40);
    iDesc.add<double>("delta", 0.04);
    iDesc.add<double>("ecut", 3.0);
    iDesc.addUntracked<unsigned int>("verbosity", 3);
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
  unsigned int NSamples_;
  unsigned int NCluster_;

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

  int nCells_;

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
  void findLeafClusters() {
    std::cout << "hierarchySize_: " << hierarchySize_ << std::endl;
    for (int i = nCells_; i < hierarchySize_; i++) {
      std::cout << " parent: " << condensed_tree_.parent[i]
		<< " child: "  << condensed_tree_.child[i]
		<< " lambda: " << condensed_tree_.lambda[i]
		<< " size: "   << condensed_tree_.size[i]
		<< " mass: "   << condensed_tree_.mass[i] << "\n";
      int cluster = condensed_tree_.child[i];
      bool find = std::find(condensed_tree_.parent.begin() + nCells_, condensed_tree_.parent.end(), cluster) == condensed_tree_.parent.end() ;
      if (!find) leaf_clusters_.emplace_back(cluster);
    }
  }

  inline int findClusterIdxInLeafClusters(int node) {
    std::vector<int>::iterator it = std::find(leaf_clusters_.begin(), leaf_clusters_.end(), node);
    if (it == leaf_clusters_.end()) return -1;
    return std::distance(leaf_clusters_.begin(), it);
  }
  
  std::vector<std::vector<int>> eom_clusters_;
  void findEOMClusters() {
    //int hierarchySize = condensed_tree_.parent.size();
    int rootNode = nCells_;

    for (int i = rootNode; i < hierarchySize_; i++) {
      int node = condensed_tree_.child[i];
      std::cout << "node: " << node << "\n";
      auto children = findChildren(node);
      float lambda_birth = condensed_tree_.lambda[i];
      float lambda_death = -1;
      auto ifind = std::find(condensed_tree_.parent.begin() + rootNode, condensed_tree_.parent.end(), node);
      if (ifind != condensed_tree_.parent.end()) {
	int idx = std::distance(condensed_tree_.parent.begin(), ifind);
	lambda_death = condensed_tree_.lambda[idx];
      }
      int nChildren = calcNChildren(i, children);
      std::cout << "lambda_death: " << lambda_death
		<< " lambda_birth: " << lambda_birth
		<< " nChildren: " << nChildren << "\n";
      if (lambda_death != -1)
	condensed_tree_.mass[i] += (lambda_death - lambda_birth) * nChildren; 
    }

    std::vector<int> selected;
    for (int i = rootNode; i < hierarchySize_; i++) {
      std::cout << "Selected: ";
      for (int s : selected) std::cout << s << " ";
      std::cout << "\n";
      int node = condensed_tree_.child[i];
      std::cout << "node: " << node << "\n";
      auto iselected = std::find(selected.begin(), selected.end(), node);
      if (iselected != selected.end()) continue;
      auto children = findChildren(node);
      int sizeChildren = children.size();
      float massChildren = calcMassChildren(i, children);
      float massNode = condensed_tree_.mass[i];
      std::cout << "massNode: " << massNode << " mass children: " << massChildren << "\n";
      //if (massNode > massChildren and sizeChildren < 5) {
      if (massNode > massChildren and sizeChildren < 5) {
	children.emplace_back(node);
	eom_clusters_.emplace_back(children);
	//selected.emplace_back(node);
	for (auto child: children) selected.emplace_back(child);
      }
    }
  }

  std::vector<int> findChildren(int node) {
    int rootNode = nCells_;
    std::vector<int> children;
    std::queue<int> q;
    q.emplace(node);
    int protection = 0;
    while (!q.empty() and protection < 100) {
      protection ++ ; 
      int inode = q.front();
      //std::cout << "inode: " << inode << "\n";
      for (int i = rootNode; i < hierarchySize_+1; i++) {
	if (inode == condensed_tree_.parent[i] and i != hierarchySize_) {
	  //std::cout << "Debug1: " << i << " " << condensed_tree_.parent[i] << "\n";
	  q.emplace(condensed_tree_.child[i]);
	  q.emplace(condensed_tree_.child[i]+1);
	  children.emplace_back(condensed_tree_.child[i]);
	  children.emplace_back(condensed_tree_.child[i]+1);
	  q.pop();
	  break;
	} else if (i == hierarchySize_) {
	  //std::cout << "Debug2\n";
	  //queue.erase(queue.begin());
	  q.pop();
	  break;
	}
      }
      //std::cout << "inode: " << inode
      //	<< " size of Queue: " << q.size()
      //	<< " size of children: " << children.size() << "\n";
      //for (int s: children) std::cout << s << " ";
      //std::cout << "\n";
    }
    return children;
  }

  int calcNChildren(int node, std::vector<int> children) {
    int rootNode = nCells_;
    int size = 0;
    //int startingNode = node;
    for (int child: children) {
      auto ifind = std::find(condensed_tree_.child.begin() + rootNode, condensed_tree_.child.end(), child);
      int idx = std::distance(condensed_tree_.child.begin(), ifind);
      size += condensed_tree_.size[idx];
      //startingNode = *idx;
    }
    return size;
  }

  float calcMassChildren(int node, std::vector<int> children) {
    int rootNode = nCells_;
    float childrenMass = 0;
    //int startingNode = node;
    for (int child: children) {
      auto ifind = std::find(condensed_tree_.child.begin() + rootNode, condensed_tree_.child.end(), child);
      int idx = std::distance(condensed_tree_.child.begin(), ifind);
      childrenMass += condensed_tree_.mass[idx];
      //startingNode = *idx;
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
  void findSafeEdges(const T& t, S& set, float delta);
  void buildMST(const T& t, S& s, float delta);
  void prepareSLTDataStructures();
  void buildSLT(S& s);
  void prepareClusterHierarchyStructures();
  void buildClusterHierarchy(int NCluster);
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
      if (cells_.layer[cell1] > 77 and cells_.layer[cell1] < 90 and cells_.layer[cell2] > 77 and  cells_.layer[cell2] < 90) d = d/2;
      if (cells_.layer[cell1] > 89 and cells_.layer[cell2] > 89) d = d/4;
    }
    return d;
  }  
};

extern template class HGCalHDBAlgoT<HGCalAngularTiles, HGCalDisjointSets>;

using HGCalHDBAlgo = HGCalHDBAlgoT<HGCalAngularTiles, HGCalDisjointSets>;

#endif

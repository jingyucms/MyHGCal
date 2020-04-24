#include "MyHGCal/HDBClustering/plugins/HGCalHDBAlgo.h"

using namespace hgcal_hdb_clustering;

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::populate(const HGCRecHitCollection& hits) {
  // loop over all hits and create the Hexel structure, skip energies below ecut

  computeThreshold();

  std::cout << "populate:" << hits.size() << "\n";

  for (unsigned int i = 0; i < hits.size(); ++i) {
    const HGCRecHit& hgrh = hits[i];
    DetId detid = hgrh.detid();
    unsigned int layerOnSide = (rhtools_.getLayerWithOffset(detid) - 1);

    // set sigmaNoise default value 1 to use kappa value directly in case of
    // sensor-independent thresholds
    float sigmaNoise = 1.f;
    int thickness_index = rhtools_.getSiThickIndex(detid);
    if (thickness_index == -1)
      thickness_index = 3;
    double storedThreshold = thresholds_[layerOnSide][thickness_index];
    sigmaNoise = v_sigmaNoise_[layerOnSide][thickness_index];

    if (hgrh.energy() < storedThreshold)
      continue;  // this sets the ZS threshold at ecut times the sigma noise
                 // for the sensor

    const GlobalPoint position(rhtools_.getPosition(detid));
    int offset = ((rhtools_.zside(detid) + 1) >> 1) * maxlayer_;
    int layer = layerOnSide + offset;

    //std::cout << "layer: " << layer << std::endl;
    if (layer < 50) continue;

    cells_.detid.emplace_back(detid);
    cells_.layer.emplace_back(layer);
    cells_.eta.emplace_back(position.eta());
    cells_.phi.emplace_back(position.phi());
    cells_.x.emplace_back(position.x());
    cells_.y.emplace_back(position.y());
    cells_.z.emplace_back(position.z());
    cells_.energy.emplace_back(hgrh.energy());
    cells_.sigmaNoise.emplace_back(sigmaNoise);
  }
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::makeClusters() {
  nCells_ = cells_.detid.size();
  std::cout << "makeClusters: " << nCells_ << "\n";
  prepareCellDataStructures();
  T tile;
  tile.clear();
  tile.fill(cells_.eta, cells_.phi);
  unsigned int NSamples = NSamples_;
  unsigned int NCluster = NCluster_;
  float delta = delta_;
  calculateDensity(tile, NSamples, delta);
  S set;
  buildMST(tile, set, delta);
  prepareSLTDataStructures();
  buildSLT(set);
  prepareClusterHierarchyStructures();
  buildClusterHierarchy(NCluster);
  assignClusters();
}

template <typename T, typename S>
std::vector<reco::BasicCluster> HGCalHDBAlgoT<T, S>::getClusters() {
  //int totalNumberOfClusters = leaf_clusters_.size() + 1;
  int totalNumberOfClusters = eom_clusters_.size()+1;
  std::cout << "totalNumberOfClusters: " << totalNumberOfClusters << "\n";
  clusters_v_.resize(totalNumberOfClusters+1);  

  std::vector<std::pair<DetId, float>> thisCluster;

  std::vector<std::vector<int>> cellsIdInCluster;
  cellsIdInCluster.resize(totalNumberOfClusters);

  for (int i = 0; i < nCells_; ++i) {
    auto clusterIndex = cells_.label[i] + 1;  // 0 is noise
    //if (clusterIndex != -1)
    //std::cout << clusterIndex << " debug1\n";
    cellsIdInCluster[clusterIndex].push_back(i);
  }

  for (int i = 0; i < totalNumberOfClusters; i++) {
    //std::cout << i << " debug1\n";
    auto cl = cellsIdInCluster[i];
    //std::cout << cl.size() << " debug2\n";
    float energy = 0.f;
    math::XYZPoint position = calculatePosition(cl); // undefined for now

    for (auto cellIdx : cl) {
      energy += cells_.energy[cellIdx];
      thisCluster.emplace_back(cells_.detid[cellIdx], 1.f);
    }
    //std::cout << i << " debug3\n";
    clusters_v_[i] =
      reco::BasicCluster(energy, position, reco::CaloID::DET_HGCAL_ENDCAP, thisCluster, reco::CaloCluster::AlgoId::undefined);
    thisCluster.clear();
  }
  return clusters_v_;
}

template <typename T, typename S>
math::XYZPoint HGCalHDBAlgoT<T, S>::calculatePosition(const std::vector<int>& cl) const {
  float z = 0.f;
  float x = 0.f;
  float y = 0.f;

  int minLayer = 100;
  for (auto idx : cl) {
    if (cells_.layer[idx] < minLayer) minLayer = cells_.layer[idx];
  }

  std::vector<int> cellsOnLayer;
  for (auto idx : cl) {
    if (cells_.layer[idx] == minLayer) cellsOnLayer.push_back(idx);
  }

  std::vector<float> X;
  std::vector<float> Y;
  for (auto idx : cellsOnLayer) {
    X.push_back(cells_.x[idx]);
    Y.push_back(cells_.y[idx]);
    z = cells_.z[idx];
  }

  x = std::accumulate( X.begin(), X.end(), 0.0) / X.size();
  y = std::accumulate( Y.begin(), Y.end(), 0.0) / Y.size();

  return math::XYZPoint(x, y, z);
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::prepareCellDataStructures() {
  cells_.rho.resize(nCells_, 9999);
  cells_.label.resize(nCells_, -1);
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::calculateDensity(const T& tile, int NSamples, float delta) {
  std::cout << "calculateDensity: " << nCells_ << "\n";
  for (int i = 0; i < nCells_; i++) {
    std::array<int, 4> search_box = tile.searchBox(cells_.eta[i] - delta,
						   cells_.eta[i] + delta,
						   cells_.phi[i] - delta,
						   cells_.phi[i] + delta);
    std::vector<float> distances;
    for (int etaBin = search_box[0]; etaBin < search_box[1] + 1; ++etaBin) {
      for (int phiBin = search_box[2]; phiBin < search_box[3] + 1; ++phiBin) {
	int binId = tile.getGlobalBinByBinEtaPhi(etaBin, phiBin);
	size_t binSize = tile[binId].size();
	
	for (unsigned int j = 0; j < binSize; j++) {
	  unsigned int otherId = tile[binId][j];

	  distances.emplace_back(distance(i, otherId));
	}
      }
    }
    std::sort (distances.begin(), distances.end());
    //std::cout << distances.size() << "\n";
    if (distances.size() >= 10)
      cells_.rho[i] = distances[NSamples];
    //if (distances.size() > 1) std::cout << "nearest: " << distances[0] << " " << distances[1] << "\n";
    //std::cout << "density: " << cells_.rho[i] << "\n";
  }
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::findSafeEdges(const T& t, S& set, float delta) {

  for (int i = 0; i < nCells_; i++) {
    Graph safeEdge;
    float d(9999);
    safeEdge.src = i;
    safeEdge.dest = nCells_-1;
    //safeEdge.dest = 0;
    safeEdge.weight = d; 
    std::array<int, 4> search_box = t.searchBox(cells_.eta[i] - delta,
						cells_.eta[i] + delta,
						cells_.phi[i] - delta,
						cells_.phi[i] + delta);
    //std::cout << "1- safeEdge.weight: " << safeEdge.weight << std::endl;
    for (int etaBin = search_box[0]; etaBin < search_box[1] + 1; ++etaBin) {
      for (int phiBin = search_box[2]; phiBin < search_box[3] + 1; ++phiBin) {
	int binId = t.getGlobalBinByBinEtaPhi(etaBin, phiBin);
	size_t binSize = t[binId].size();
	
	for (unsigned int j = 0; j < binSize; j++) {
	  unsigned int otherId = t[binId][j];
	  //if (otherId == i) continue;
	  if (set.find(otherId) == set.find(i)) continue;
	  float maxRho = std::max(cells_.rho[i], cells_.rho[otherId]);
	  float weight = std::max(maxRho, distance(i, otherId));
	  if (weight < safeEdge.weight) {
	    safeEdge.dest = otherId;
	    safeEdge.weight = weight;
	  }
	}
      }
    }
    //std::cout << "2- safeEdge.weight: " << safeEdge.weight << std::endl;
    safeEdges_.emplace_back(safeEdge);
  }
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::buildMST(const T& tile, S& set, float delta) {
  int nTrees = nCells_;
  set.init(nCells_);
  while (nTrees > 1) {
    // find save edges
    //std::cout << "buildMST: " << nTrees << "\n";
    //std::cout << safeEdges_.size() << " safeEdge-1\n";
    findSafeEdges(tile, set, delta);
    //std::cout << safeEdges_.size() << " safeEdge-2\n";

    // find cheapest edges
    std::vector<int> cheapest;
    cheapest.resize(nCells_, -1);
    for (int i = 0; i<nCells_; i++) {
      Graph safeEdge = safeEdges_[i];
      int src = safeEdge.src;
      //int dest = safeEdge.dest;
      //if (dest == -1) continue;
      float weight = safeEdge.weight;
      int set_src = set.findWithPathCompression(src);
      int cheapest_in_set = cheapest[set_src];
      if (cheapest_in_set == -1) {
	cheapest[set_src] = i;
      } else if (cheapest_in_set != -1 and weight < safeEdges_[cheapest_in_set].weight) {
	cheapest[set_src] = i;
      }
      //std::cout << "cell: " << i << " weight: " << weight << " set: " << set_src << " cheapest: " << cheapest[set_src] << " weight: " << safeEdges_[cheapest[set_src]].weight << "\n";
    }
    for (int i = 0; i<nCells_; i++) {
      //std::cout << "cheapest " << i << ": " << cheapest[i] << "\n";
      if (cheapest[i] == -1) continue;
      Graph safeEdge = safeEdges_[i];
      int src = safeEdge.src;
      int dest = safeEdge.dest;
      //std::cout << src << " " << dest << "\n";
      int set_src = set.findWithPathCompression(src);
      int set_dest = set.findWithPathCompression(dest);
      if (set_src != set_dest) {
	set.Union(set_src, set_dest);
	min_span_tree_.emplace_back(safeEdge);
	nTrees--;
      }
    }
    safeEdges_.clear();
  }
  std::sort(min_span_tree_.begin(), min_span_tree_.end());
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::prepareSLTDataStructures() {
  int treeSize = nCells_ - 1;
  single_linkage_tree_.parent.resize(treeSize, -1);
  single_linkage_tree_.left.resize(treeSize, -1);
  single_linkage_tree_.right.resize(treeSize, -1);
  single_linkage_tree_.weight.resize(treeSize, -1);
  single_linkage_tree_.size.resize(treeSize, 1);
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::buildSLT(S& set) {
  int nEdges = nCells_-1;
  set.init(nCells_);
  int parent = nEdges;
  std::cout << "buildSLT: " << min_span_tree_.size() << "\n";
  //std::cout << "parent: " << nEdges << "\n";
  for (int i = 0; i < nEdges; i++) {
    std::cout << "MST:      " << min_span_tree_[i].src
	      << " " << min_span_tree_[i].dest
	      << " " << min_span_tree_[i].weight << "\n";
    parent++;
    single_linkage_tree_.weight[i] = min_span_tree_[i].weight;
    single_linkage_tree_.parent[i] = parent;
    int src = min_span_tree_[i].src;
    int dest = min_span_tree_[i].dest;
    int src_rank = set.getRootRank(src);
    int dest_rank = set.getRootRank(dest);
    if (src_rank == 1 and dest_rank == 1) {
      single_linkage_tree_.left[i] = src;
      single_linkage_tree_.right[i] = dest;
    }
    if (src_rank == 1 and dest_rank > 1) {
      single_linkage_tree_.left[i] = src;
      single_linkage_tree_.right[i] = set.getLinkageNode(dest);
    }
    if (src_rank > 1 and dest_rank == 1) {
      single_linkage_tree_.left[i] = set.getLinkageNode(src);
      single_linkage_tree_.right[i] = dest;
    }
    if (src_rank > 1 and dest_rank > 1) {
      single_linkage_tree_.left[i] = set.getLinkageNode(src);
      single_linkage_tree_.right[i] = set.getLinkageNode(dest);
    }
    
    set.Union(src, dest);
    //std::cout << "src dest: " << src << " " << dest << "\n";
    set.setLinkageNode(src, parent);
    single_linkage_tree_.size[i] = set.getRootRank(src);
    //std::cout << "size " << set.getRootRank(src) << "\n";
    std::cout << "SLT: " << single_linkage_tree_.parent[i]
	      << " " << single_linkage_tree_.left[i]
	      << " " << single_linkage_tree_.right[i]
	      << " " << single_linkage_tree_.weight[i]
	      << " " << single_linkage_tree_.size[i] << "\n";
  }
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::prepareClusterHierarchyStructures() {
  condensed_tree_.parent.resize(nCells_, -1);
  condensed_tree_.child.resize(nCells_, -1);
  condensed_tree_.lambda.resize(nCells_, -1);
  condensed_tree_.size.resize(nCells_, 1);
  condensed_tree_.mass.resize(nCells_, 0);
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::buildClusterHierarchy(int NCluster) {
  //auto nSteps = single_linkage_tree_.parent.size();
  int nSteps = nCells_ - 1;
  std::cout << "buildClusterHierarchy: " << nSteps << "\n";
  // initialize a queue
  nodes_in_slt_.clear();
  nodes_in_hierarchy_.clear();
  int root_in_slt = single_linkage_tree_.parent[nSteps-1];
  nodes_in_slt_.emplace_back(root_in_slt);
  int root_in_hierarchy = nCells_;
  nodes_in_hierarchy_.emplace_back(root_in_hierarchy);
  
  int j = 0;
  int cluster_hierarchy_count = root_in_hierarchy;
  bool isRoot;
  for (int i = nSteps-1; i>=0; i--) {
    if (i == nSteps - 1 ) isRoot = true;
    else isRoot = false;
    int left = single_linkage_tree_.left[i];
    int right = single_linkage_tree_.right[i];
    int parent_in_slt = single_linkage_tree_.parent[i];
    float lambda = 1/single_linkage_tree_.weight[i];
    bool split = (left >= nCells_ and right >= nCells_ and
		  getNodeSizeInSLT(left, i) > NCluster and
		  getNodeSizeInSLT(right, i) > NCluster and
		  lambda > 1.);
    //std::cout << i << " split: " << split
    //	      << " parent: " << parent_in_slt
    //	      << " left: " << left
    //	      << " right: " << right << "\n";
    int previous_split_in_slt, previous_split_in_hierarchy;
    if (isRoot) {
      previous_split_in_slt = root_in_slt;
      previous_split_in_hierarchy = root_in_hierarchy;
    } else {
      previous_split_in_slt = findPreSplitInSLT(parent_in_slt, i);
      previous_split_in_hierarchy = findPreSplitInHierarchy(previous_split_in_slt);
    }
    //std::cout << "previous split in slt/den: " << previous_split_in_slt << " " << previous_split_in_hierarchy << "\n";
    if (split) {
      std::cout << " ------ Split!" << "\n";
      for (int k = 1; k <=2; k++) {
	condensed_tree_.parent.emplace_back(previous_split_in_hierarchy);
	condensed_tree_.lambda.emplace_back(lambda);
	condensed_tree_.size.emplace_back(0);
	cluster_hierarchy_count++;
	std::cout << "left: " << left
		  << " right: " << right
		  << " parent: " << parent_in_slt
		  << " lambda: " << lambda
		  << " previous_split_in_hierarchy: " << previous_split_in_hierarchy
		  << " cluster_hierarchy_count: " << cluster_hierarchy_count << "\n";
	condensed_tree_.child.emplace_back(cluster_hierarchy_count);
	if (k == 1) nodes_in_slt_.emplace_back(left);
	else nodes_in_slt_.emplace_back(right);
	nodes_in_hierarchy_.emplace_back(cluster_hierarchy_count);
      }
    } else if (left < nCells_ or right < nCells_ ) {
      condensed_tree_.parent[j] = previous_split_in_hierarchy;
      condensed_tree_.lambda[j] = lambda;
      condensed_tree_.size[j] = 1;
      int parent_index_in_hierarchy = previous_split_in_hierarchy - 1;
      if (left < nCells_) {
	condensed_tree_.child[j] = left;
	if (previous_split_in_hierarchy != root_in_hierarchy) {
	  condensed_tree_.size[parent_index_in_hierarchy] ++;
	  condensed_tree_.mass[parent_index_in_hierarchy] += lambda - condensed_tree_.lambda[parent_index_in_hierarchy];
	  //std::cout << " parent_index_in_hierarchy: " << parent_index_in_hierarchy
	  //	    << " " << condensed_tree_.child[parent_index_in_hierarchy] << "\n";
	}
	j++;
      }
      if (right < nCells_) {
	condensed_tree_.child[j] = right;
	
	if (previous_split_in_hierarchy != root_in_hierarchy) {
	  condensed_tree_.size[parent_index_in_hierarchy] ++;
	  condensed_tree_.mass[parent_index_in_hierarchy] += lambda - condensed_tree_.lambda[parent_index_in_hierarchy];
	  //std::cout << " parent_index_in_hierarchy: " << parent_index_in_hierarchy
	  //	    << " " << condensed_tree_.child[parent_index_in_hierarchy] << "\n";
	}
	j++;
      }
    }
  }
  hierarchySize_ = condensed_tree_.parent.size();
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::assignClusters() {
  std::cout << "assignClusters: " << condensed_tree_.parent.size() << "\n";
  //leaf_clusters_.clear();
  findLeafClusters();
  //std::cout << "leafClusters: " << leaf_clusters_.size() << "\n";
  eom_clusters_.clear();
  findEOMClusters();
  std::cout << "eomClusters: " << eom_clusters_.size() << "\n";
  for (int i = 0; i < nCells_; i++) {
    int cell = condensed_tree_.child[i];
    int parent = condensed_tree_.parent[i];
    //int label = findClusterIdxInLeafClusters(parent);
    int label = findClusterIdxInEOMClusters(parent);
    cells_.label[cell] = label;
  }
  //std::cout << "Number of Cluster: " << *max_element(cells_.label.begin(), cells_.label.end()) << "\n"; 
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::computeThreshold() {

  std::vector<double> dummy;
  const unsigned maxNumberOfThickIndices = 3;
  dummy.resize(maxNumberOfThickIndices + 1, 0);  // +1 to accomodate for the Scintillators
  thresholds_.resize(maxlayer_, dummy);
  v_sigmaNoise_.resize(maxlayer_, dummy);

  for (unsigned ilayer = 1; ilayer <= maxlayer_; ++ilayer) {
    for (unsigned ithick = 0; ithick < maxNumberOfThickIndices; ++ithick) {
      float sigmaNoise = 0.001f * fcPerEle_ * nonAgedNoises_[ithick] * dEdXweights_[ilayer] / (fcPerMip_[ithick] * thicknessCorrection_[ithick]);
      thresholds_[ilayer - 1][ithick] = sigmaNoise * ecut_;
      v_sigmaNoise_[ilayer - 1][ithick] = sigmaNoise;
      LogDebug("HGCalHDBAlgo") << "ilayer: " << ilayer << " nonAgedNoises: " << nonAgedNoises_[ithick]
				<< " fcPerEle: " << fcPerEle_ << " fcPerMip: " << fcPerMip_[ithick]
				<< " noiseMip: " << fcPerEle_ * nonAgedNoises_[ithick] / fcPerMip_[ithick]
				<< " sigmaNoise: " << sigmaNoise << "\n";
    }
    float scintillators_sigmaNoise = 0.001f * noiseMip_ * dEdXweights_[ilayer];
    thresholds_[ilayer - 1][maxNumberOfThickIndices] = ecut_ * scintillators_sigmaNoise;
    v_sigmaNoise_[ilayer - 1][maxNumberOfThickIndices] = scintillators_sigmaNoise;
    LogDebug("HGCalHDBAlgo") << "ilayer: " << ilayer << " noiseMip: " << noiseMip_
			     << " scintillators_sigmaNoise: " << scintillators_sigmaNoise << "\n";
  }
}

template <typename T, typename S>
Density HGCalHDBAlgoT<T, S>::getDensity() {
  return density_;
}

template class HGCalHDBAlgoT<HGCalAngularTiles, HGCalDisjointSets>;

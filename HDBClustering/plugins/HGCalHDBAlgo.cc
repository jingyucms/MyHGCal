#include "MyHGCal/HDBClustering/plugins/HGCalHDBAlgo.h"

#include <chrono>

using namespace hgcal_hdb_clustering;

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::populate(const HGCRecHitCollection& hits) {
  // loop over all hits and create the Hexel structure, skip energies below ecut

  computeThreshold();

  if (verbosity_ > 0)
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
  if (verbosity_ > 0)
    std::cout << "makeClusters: " << nCells_ << "\n";
  prepareCellDataStructures();
  T tile;
  tile.clear();
  tile.fill(cells_.eta, cells_.phi);
  int NSamples = NSamples_;
  int NCluster = NCluster_;
  float delta = delta_;

  //std::chrono::system_clock::time_point now = std::chrono::system_clock::now();

  std::chrono::steady_clock::time_point start1 = std::chrono::steady_clock::now();
  calculateDensity(tile, NSamples, delta);
  std::chrono::steady_clock::time_point start2 = std::chrono::steady_clock::now();
  if (verbosity_ == 0)
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(start2 - start1).count() << "\n";
  S set;
  
  buildMST(tile, set, delta);
  std::chrono::steady_clock::time_point start3 = std::chrono::steady_clock::now();
  if (verbosity_ == 0)
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(start3 - start2).count() << "\n";

  
  prepareSLTDataStructures();
  buildSLT(set);
  std::chrono::steady_clock::time_point start4 = std::chrono::steady_clock::now();
  if (verbosity_ == 0)
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(start4 - start3).count() << "\n";
  
  prepareClusterHierarchyStructures();
  buildClusterHierarchy(NCluster);

  std::chrono::steady_clock::time_point start5 = std::chrono::steady_clock::now();
  if (verbosity_ == 0)
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(start5 - start4).count() << "\n";
  
  assignClusters();
  std::chrono::steady_clock::time_point start6 = std::chrono::steady_clock::now();
  if (verbosity_ == 0)
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(start6 - start5).count() << "\n";

  if (verbosity_ == 0)
    std::cout << "Total Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(start6 - start1).count() << "\n";
}

template <typename T, typename S>
std::vector<reco::BasicCluster> HGCalHDBAlgoT<T, S>::getClusters() {
  int totalNumberOfClusters(0);
  if (Algo_ == 0)
    totalNumberOfClusters = leaf_clusters_.size()+1;
  else if (Algo_ == 1)
    totalNumberOfClusters = eom_clusters_.size()+1;

  if (verbosity_ > 0)
    std::cout << "totalNumberOfClusters: " << totalNumberOfClusters << "\n";
  clusters_v_.resize(totalNumberOfClusters+1);  

  std::vector<std::pair<DetId, float>> thisCluster;

  std::vector<std::vector<int>> cellsIdInCluster;
  cellsIdInCluster.resize(totalNumberOfClusters);

  for (int i = 0; i < nCells_; ++i) {
    auto clusterIndex = cells_.label[i] + 1;  // 0 is noise
    cellsIdInCluster[clusterIndex].push_back(i);
  }

  for (int i = 0; i < totalNumberOfClusters; i++) {
    auto cl = cellsIdInCluster[i];
    float energy = 0.f;
    math::XYZPoint position = calculatePosition(cl); // undefined for now

    for (auto cellIdx : cl) {
      energy += cells_.energy[cellIdx];
      thisCluster.emplace_back(cells_.detid[cellIdx], 1.f);
    }

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

  if (verbosity_ > 0)
    std::cout << "calculatePosition: " << x << " " << y << " " << z << "\n";
  return math::XYZPoint(x, y, z);
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::prepareCellDataStructures() {
  cells_.rho.resize(nCells_, 9999);
  cells_.label.resize(nCells_, -1);
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::calculateDensity(const T& tile, int NSamples, float delta) {
  if (verbosity_ > 0)
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
	  //if (i == 327 || i == 393) 
	  //std::cout << "Debug Density: " << i << " "
	  //	      << otherId << " "
	  //	      << distance(i, otherId) << "\n";
	  distances.emplace_back(distance(i, otherId));
	}
      }
    }
    std::sort (distances.begin(), distances.end());
    if (distances.size() >= 10)
      cells_.rho[i] = distances[NSamples];
    else cells_.rho[i] = distances[distances.size()-1]*10/distances.size();
    //if (i == 327 || i == 666 || i == 448) 
      //std::cout << "DEBUG--" << i << ": " 
      //	<< cells_.eta[i] << " "
      //	<< cells_.phi[i] << " "
      //	<< cells_.layer[i] << " "
      //	<< cells_.rho[i] << " "
      //	<< distances.size() << "\n";
  }
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::findSafeEdges(const T& t, S& set, float delta) {

  for (int i = 0; i < nCells_; i++) {
    Graph safeEdge;
    float d(9999);
    safeEdge.src = i;
    safeEdge.dest = nCells_-1;
    safeEdge.weight = d; 
    std::array<int, 4> search_box = t.searchBox(cells_.eta[i] - delta,
						cells_.eta[i] + delta,
						cells_.phi[i] - delta,
						cells_.phi[i] + delta);
    for (int etaBin = search_box[0]; etaBin < search_box[1] + 1; ++etaBin) {
      for (int phiBin = search_box[2]; phiBin < search_box[3] + 1; ++phiBin) {
	int binId = t.getGlobalBinByBinEtaPhi(etaBin, phiBin);
	size_t binSize = t[binId].size();
	
	for (unsigned int j = 0; j < binSize; j++) {
	  unsigned int otherId = t[binId][j];
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
    //if (i == 327 || i == 666 || i == 448) 
    //  std::cout << "DEBUG!!!: " << safeEdge.src << " " << safeEdge.dest << " " << safeEdge.weight << "\n";
    
    safeEdges_.emplace_back(safeEdge);
  }
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::buildMST(const T& tile, S& set, float delta) {
  int nTrees = nCells_;
  set.init(nCells_);
  while (nTrees > 1) {
    if (verbosity_ > 0)
      std::cout << "buildMST: " << nTrees << "\n";
    findSafeEdges(tile, set, delta);

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
      if (verbosity_ > 3)
	std::cout << "cell: " << i << " weight: " << weight << " set: " << set_src << " cheapest: " << cheapest[set_src] << " weight: " << safeEdges_[cheapest[set_src]].weight << "\n";
    }
    for (int i = 0; i<nCells_; i++) {
      if (verbosity_ > 3)
	std::cout << "cheapest " << i << ": " << cheapest[i] << "\n";
      if (cheapest[i] == -1) continue;
      Graph safeEdge = safeEdges_[cheapest[i]];
      int src = safeEdge.src;
      int dest = safeEdge.dest;
      if (verbosity_ > 3)
	std::cout << src << " " << dest << "\n";
      int set_src = set.findWithPathCompression(src);
      int set_dest = set.findWithPathCompression(dest);
      if (set_src != set_dest) {
	set.Union(src, dest);
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
  if (verbosity_ > 0)
    std::cout << "buildSLT: " << min_span_tree_.size() << "\n";
  for (int i = 0; i < nEdges; i++) {
    if (verbosity_ > 3)
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
    set.setLinkageNode(src, parent);
    single_linkage_tree_.size[i] = set.getRootRank(src);
    //std::cout << "size " << set.getRootRank(src) << "\n";
    if (verbosity_ > 3)
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
  if (verbosity_ > 0)
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
    int cell1 = min_span_tree_[i].src;
    int cell2 = min_span_tree_[i].dest;
    float d = distance(cell1, cell2);
    float rho1 = cells_.rho[cell1];
    float rho2 = cells_.rho[cell2];
    if (verbosity_ > 3) {
      std::cout << "Debugging Condensed Tree for edge: " << i << "\n";
      std::cout << cell1 << " "
		<< cell2 << " "
		<< d << " "
		<< rho1 << " "
		<< rho2 << "\n";
    }

    float minimum(10);
    bool split = ((lambda > minimum and left >= nCells_ and right >= nCells_ and
		  getNodeSizeInSLT(left, i) > NCluster and
		  getNodeSizeInSLT(right, i) > NCluster) or
		  (lambda <= minimum and left >= nCells_ and right >= nCells_ and
		  getNodeSizeInSLT(left, i) > 0 and
		  getNodeSizeInSLT(right, i) > 0));
    

    int previous_split_in_slt, previous_split_in_hierarchy;
    if (isRoot) {
      previous_split_in_slt = root_in_slt;
      previous_split_in_hierarchy = root_in_hierarchy;
    } else {
      if (verbosity_ > 2) std::cout << "------ Building Hierarchy!" << "\n";
      previous_split_in_slt = findPreSplitInSLT(parent_in_slt, i);
      previous_split_in_hierarchy = findPreSplitInHierarchy(previous_split_in_slt);
      if (verbosity_ > 2) {
	std::cout << previous_split_in_slt << " " << previous_split_in_hierarchy << "\n";
	std::cout << "nodes_in_slt_: ";
	for (int i : nodes_in_slt_)
	  std::cout << " " << i;
	std::cout << "\n";
	std::cout << "nodes_in_hie_: ";
	for (int i : nodes_in_hierarchy_)
	  std::cout << " " << i;
	std::cout << "\n";
      }
    }
    if (split) {
      if (verbosity_ > 2)
	std::cout << " ------ Split!" << "\n";
      for (int k = 1; k <=2; k++) {
	condensed_tree_.parent.emplace_back(previous_split_in_hierarchy);
	condensed_tree_.lambda.emplace_back(lambda);
	condensed_tree_.size.emplace_back(0);
	cluster_hierarchy_count++;
	if (verbosity_ > 2)
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
	}
	//if (parent_index_in_hierarchy == 3285) 
	//std::cout << condensed_tree_.parent[parent_index_in_hierarchy] << " "
	//	    << lambda << " "
	//	    << condensed_tree_.lambda[parent_index_in_hierarchy] << "\n";
	j++;
      }
      if (right < nCells_) {
	condensed_tree_.child[j] = right;
	
	if (previous_split_in_hierarchy != root_in_hierarchy) {
	  condensed_tree_.size[parent_index_in_hierarchy] ++;
	  condensed_tree_.mass[parent_index_in_hierarchy] += lambda - condensed_tree_.lambda[parent_index_in_hierarchy];
	}
	//if (parent_index_in_hierarchy == 3285) 
	//std::cout << condensed_tree_.parent[parent_index_in_hierarchy] << " "
	//	    << lambda << " "
	//	    << condensed_tree_.lambda[parent_index_in_hierarchy] << "\n";
	j++;
      }
    }
  }
  hierarchySize_ = condensed_tree_.parent.size();
}

template <typename T, typename S>
void HGCalHDBAlgoT<T, S>::assignClusters() {
  if (verbosity_ > 0)
    std::cout << "assignClusters: " << condensed_tree_.parent.size() << "\n";
  leaf_clusters_.clear();
  findLeafClusters();
  if (verbosity_ > 0)
    std::cout << "leafClusters: " << leaf_clusters_.size() << "\n";
  eom_clusters_.clear();
  findEOMClusters();
  if (verbosity_ > 0)
    std::cout << "eomClusters: " << eom_clusters_.size() << "\n";
  for (int i = 0; i < nCells_; i++) {
    int cell = condensed_tree_.child[i];
    int parent = condensed_tree_.parent[i];
    int label(0);
    if (Algo_ == 0)
      label = findClusterIdxInLeafClusters(parent);
    else if (Algo_ == 1)
      label = findClusterIdxInEOMClusters(parent);
    //if (cells_.rho[cell] <  0.1 or
	//(cells_.rho[cell] > 0.1 and cells_.energy[cell] > 9 * cells_.sigmaNoise[cell]) or
	//(cells_.rho[cell] < 1 and cells_.layer[cell] > 90))
    if (cells_.rho[cell] < 1)
      cells_.label[cell] = label;
  }
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

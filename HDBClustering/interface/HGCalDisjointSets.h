#ifndef MyHGCal_HDBClustering_HGCalDisjointSets_h
#define MyHGCal_HDBClustering_HGCalDisjointSets_h

class HGCalDisjointSets {
 public:
  
  void init(int Nnode) {
    subsets_.clear();
    subsets_.resize(Nnode);
    for (int v = 0; v < Nnode; ++v) {
      subsets_[v].parent = v;
      subsets_[v].rank = 1;
      subsets_[v].linkageNode = -1;
    }
  }

  void clear() {
    subsets_.clear();
  }

  struct Subset { 
    int parent; 
    int rank;
    int linkageNode;
  };

  int find(int i) {
    if (subsets_[i].parent == i)  
      return i;  
    return find(subsets_[i].parent);  
  }
  
  int findWithPathCompression(int i) {
    if (subsets_[i].parent != i) 
      subsets_[i].parent = findWithPathCompression(subsets_[i].parent); 
  
    return subsets_[i].parent; 
  }
  
  void Union(int x, int y){
    int xroot = findWithPathCompression(x); 
    int yroot = findWithPathCompression(y); 
  
    if (subsets_[xroot].rank < subsets_[yroot].rank) {
      subsets_[xroot].parent = yroot;
      subsets_[yroot].rank += subsets_[xroot].rank;
    }
    else if (subsets_[xroot].rank > subsets_[yroot].rank) {
      subsets_[yroot].parent = xroot;
      subsets_[xroot].rank +=subsets_[yroot].rank;
    }
    else { 
      subsets_[yroot].parent = xroot; 
      subsets_[xroot].rank+=subsets_[yroot].rank; 
    } 
  }

  std::vector<struct Subset> getSet(){
    return subsets_;
  }

  int getRootRank(int i) {
    int root = findWithPathCompression(i);
    return subsets_[root].rank;
  }

  void setLinkageNode(int i, int node) {
    int root = findWithPathCompression(i);
    subsets_[root].linkageNode = node;
  }

  int getLinkageNode(int i) {
    int root = findWithPathCompression(i);
    return subsets_[root].linkageNode;
  }

 private:
  std::vector<struct Subset> subsets_;
  
};

#endif

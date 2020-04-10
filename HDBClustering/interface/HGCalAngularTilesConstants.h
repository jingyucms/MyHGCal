#ifndef MyHGCal_HDBClustering_HGCalTilesConstants_h
#define MyHGCal_HDBClustering_HGCalTilesConstants_h

#include "DataFormats/Math/interface/constexpr_cmath.h"
#include <cstdint>
#include <array>

struct HGCalAngularTilesConstants {
  
  static constexpr float tileSize = 0.1f;
  static constexpr float minEta = -3.2f;
  static constexpr float maxEta = 3.2f;
  //To properly construct search box for cells in phi=[-3.2,-3.1] and [3.1,3.2],
  // cells in phi=[3.1 ,3.2] are copied to the first bin and cells in phi=[-3.2,-3.1] are copied to the last bin
  static constexpr float minPhi = -3.3f;
  static constexpr float maxPhi = 3.3f;
  static constexpr int nColumns = reco::ceil((maxEta - minEta) / tileSize);
  static constexpr int nRows = reco::ceil((maxPhi - minPhi) / tileSize);
  static constexpr int nTiles = nColumns * nRows;
};

#endif

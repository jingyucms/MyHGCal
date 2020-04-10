#ifndef MyHGCal_HDBClustering_HGCalAngularTiles_h
#define MyHGCal_HDBClustering_HGCalAngularTiles_h

#include "MyHGCal/HDBClustering/interface/HGCalAngularTilesConstants.h"

//#include <vector>
//#include <array>
//#include <cmath>
//#include <algorithm>
//#include <cassert>

template <typename T>
class HGCalAngularTilesT {
public:
  void fill(const std::vector<float>& eta,
            const std::vector<float>& phi) {
    
    auto cellsSize = eta.size();
    for (unsigned int i = 0; i < cellsSize; ++i) {
      tiles_[getGlobalBinEtaPhi(eta[i], phi[i])].push_back(i);
      // Copy cells in phi=[-3.2,-3.1] to the last bin
      if (getPhiBin(phi[i]) == mPiPhiBin) {
	tiles_[getGlobalBinEtaPhi(eta[i], phi[i] + 2 * M_PI)].push_back(i);
      }
      // Copy cells in phi=[3.1,3.2] to the first bin
      if (getPhiBin(phi[i]) == pPiPhiBin) {
	tiles_[getGlobalBinEtaPhi(eta[i], phi[i] - 2 * M_PI)].push_back(i);
      }
    }
  }

  int getEtaBin(float eta) const {
    constexpr float etaRange = T::maxEta - T::minEta;
    static_assert(etaRange >= 0.);
    constexpr float r = T::nColumns / etaRange;
    int etaBin = (eta - T::minEta) * r;
    etaBin = std::clamp(etaBin, 0, T::nColumns - 1);
    return etaBin;
  }

  int getPhiBin(float phi) const {
    constexpr float phiRange = T::maxPhi - T::minPhi;
    static_assert(phiRange >= 0.);
    constexpr float r = T::nRows / phiRange;
    int phiBin = (phi - T::minPhi) * r;
    phiBin = std::clamp(phiBin, 0, T::nRows - 1);
    return phiBin;
  }

  int mPiPhiBin = getPhiBin(-M_PI);
  int pPiPhiBin = getPhiBin(M_PI);

  
  int getGlobalBinEtaPhi(float eta, float phi) const {
    return getEtaBin(eta) + getPhiBin(phi) * T::nColumns;
  }

  int getGlobalBinByBinEtaPhi(int etaBin, int phiBin) const {
    return etaBin + phiBin * T::nColumns;
  }

  std::array<int, 4> searchBox(float etaMin, float etaMax, float phiMin, float phiMax) const {
    int etaBinMin = getEtaBin(etaMin);
    int etaBinMax = getEtaBin(etaMax);
    int phiBinMin = getPhiBin(phiMin);
    int phiBinMax = getPhiBin(phiMax);
    return std::array<int, 4>({{etaBinMin, etaBinMax, phiBinMin, phiBinMax}});
  }

  void clear() {
    for (auto& t : tiles_)
      t.clear();
  }

  const std::vector<int>& operator[](int globalBinId) const { return tiles_[globalBinId]; }

private:
  std::array<std::vector<int>, T::nTiles> tiles_;
};

using HGCalAngularTiles = HGCalAngularTilesT<HGCalAngularTilesConstants>;

#endif

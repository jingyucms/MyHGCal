#ifndef MyHGCal_HDBClustering_HGCalHDBClusteringAlgoFactory_H
#define MyHGCal_HDBClustering_HGCalHDBClusteringAlgoFactory_H

#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MyHGCal/HDBClustering/interface/HGCalHDBClusteringAlgoBase.h"

typedef edmplugin::PluginFactory<HGCalHDBClusteringAlgoBase * (const edm::ParameterSet&)> HGCalHDBClusteringAlgoFactory;

#endif

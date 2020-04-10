#include "MyHGCal/HDBClustering/interface/HGCalHDBClusteringAlgoFactory.h"
#include "MyHGCal/HDBClustering/interface/HGCalHDBClusteringAlgoBase.h"
#include "MyHGCal/HDBClustering/plugins/HGCalHDBAlgo.h"

#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"

EDM_REGISTER_VALIDATED_PLUGINFACTORY(HGCalHDBClusteringAlgoFactory, "HGCalHDBClusteringAlgoFactory");
DEFINE_EDM_VALIDATED_PLUGIN(HGCalHDBClusteringAlgoFactory, HGCalHDBAlgo, "HDBv1");

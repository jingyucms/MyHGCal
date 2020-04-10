import FWCore.ParameterSet.Config as cms

from MyHGCal.HDBClustering.hgcalHDBClusters_cfi import hgcalHDBClusters

from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import dEdX, HGCalRecHit

from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit

from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import fC_per_ele, hgceeDigitizer, hgchebackDigitizer, hfnoseDigitizer

hgcalHDBClusters.plugin.dEdXweights = cms.vdouble(dEdX.weights)
hgcalHDBClusters.plugin.fcPerMip = cms.vdouble(HGCalUncalibRecHit.HGCEEConfig.fCPerMIP)
hgcalHDBClusters.plugin.thicknessCorrection = cms.vdouble(HGCalRecHit.thicknessCorrection)
hgcalHDBClusters.plugin.fcPerEle = cms.double(fC_per_ele)
hgcalHDBClusters.plugin.noises = cms.PSet(refToPSet_ = cms.string('HGCAL_noises'))
hgcalHDBClusters.plugin.noiseMip = hgchebackDigitizer.digiCfg.noise

hgcalHDBRecoSequence = cms.Sequence( HGCalUncalibRecHit+
                                     HGCalRecHit+
                                     hgcalHDBClusters )

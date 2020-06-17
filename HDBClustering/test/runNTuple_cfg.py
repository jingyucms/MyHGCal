import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9

from myParam_cff import *
filename=filename.replace("step1", "step3").replace(".root", ".root")

process = cms.Process('ANALYSIS',Phase2C9)

process.load('Configuration.Geometry.GeometryExtended2026D41Reco_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

## HGCal raw to digis
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.hgcalRawToDigis_step=cms.Path(process.hgcalDigis)

## HGCal sim hit study
process.load('Validation.HGCalValidation.hgcSimHitStudy_cfi')

## HGCal digi study
process.load('Validation.HGCalValidation.hgcDigiStudy_cfi')

## HGCal scintillator rec hit study
process.load('MyHGCal.HGCalStuties.hgcalRecHitStudyHEScintillator_cfi')

## HGCal scintillator layer cluster study
process.load('MyHGCal.HGCalStuties.hgcalLayerClusterStudyHESintillator_cfi')

## HGCal Ntuple
#process.load("RecoLocalCalo.Configuration.hgcalLocalReco_cff")
#from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import *
#from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import *
## process.load('MyHGCal.HGCalStuties.hgcalHitNtuple_cfi')
## process.hgcalHitNtuple.fcPerMip = cms.vdouble(2.06, 3.43, 5.15)
## process.hgcalHitNtuple.thicknessCorrection = cms.vdouble(0.781, 0.775, 0.769)
## process.hgcalHitNtuple.noises = cms.PSet(
##     values = cms.vdouble(2000.0, 2400.0, 2000.0)
## )
## process.hgcalHitNtuple.dEdXweights = cms.vdouble(
##             0.0, 8.894541, 10.937907, 10.937907, 10.937907, 
##             10.937907, 10.937907, 10.937907, 10.937907, 10.937907, 
##             10.932882, 10.932882, 10.937907, 10.937907, 10.938169, 
##             10.938169, 10.938169, 10.938169, 10.938169, 10.938169, 
##             10.938169, 10.938169, 10.938169, 10.938169, 10.938169, 
##             10.938169, 10.938169, 10.938169, 32.332097, 51.574301, 
##             51.444192, 51.444192, 51.444192, 51.444192, 51.444192, 
##             51.444192, 51.444192, 51.444192, 51.444192, 51.444192, 
##             69.513118, 87.582044, 87.582044, 87.582044, 87.582044, 
##             87.582044, 87.214571, 86.888309, 86.92952, 86.92952, 
##             86.92952
##         )

## HGCal ntuple
process.load("RecoLocalCalo.Configuration.hgcalLocalReco_cff")
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import dEdX, HGCalRecHit

from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit

from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import fC_per_ele, hgceeDigitizer, hgchebackDigitizer, hfnoseDigitizer


process.load('MyHGCal.HGCalStuties.hgcalClusterNtuple_cfi')
process.hgcalClusterNtuple.dEdXweights = cms.vdouble(dEdX.weights)
process.hgcalClusterNtuple.fcPerMip = cms.vdouble(HGCalUncalibRecHit.HGCEEConfig.fCPerMIP)
process.hgcalClusterNtuple.thicknessCorrection = cms.vdouble(HGCalRecHit.thicknessCorrection)
process.hgcalClusterNtuple.fcPerEle = cms.double(fC_per_ele)
process.hgcalClusterNtuple.noises = cms.PSet(refToPSet_ = cms.string('HGCAL_noises'))
process.hgcalClusterNtuple.noiseMip = hgchebackDigitizer.digiCfg.noise


process.load('MyHGCal.HGCalStuties.hgcalHitNtuple_cfi')
process.hgcalHitNtuple.dEdXweights = cms.vdouble(dEdX.weights)
process.hgcalHitNtuple.fcPerMip = cms.vdouble(HGCalUncalibRecHit.HGCEEConfig.fCPerMIP)
process.hgcalHitNtuple.thicknessCorrection = cms.vdouble(HGCalRecHit.thicknessCorrection)
process.hgcalHitNtuple.fcPerEle = cms.double(fC_per_ele)
process.hgcalHitNtuple.noises = cms.PSet(refToPSet_ = cms.string('HGCAL_noises'))
process.hgcalHitNtuple.noiseMip = hgchebackDigitizer.digiCfg.noise

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring("file:root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/HGC/HDBSCAN/"+filename)
                            fileNames = cms.untracked.vstring("file:"+filename)
                            #skipEvents=cms.untracked.uint32(1)
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("h_"+filename),
                                   closeFileFast = cms.untracked.bool(True)
                                   )


process.simhit_analysis = cms.Path(process.hgcalSimHitStudy)
process.digi_analysis=cms.Path(process.hgcalDigiStudyEE+process.hgcalDigiStudyHEF+process.hgcalDigiStudyHEB)
process.rechit_analysis=cms.Path(process.hgcalRecHitStudyHEScintillator)
process.layercluster_scint_analysis=cms.Path(process.hgcalLayerClusterStudyHESintillator)
process.hgcal_rechit_ntuple=cms.Path(process.hgcalHitNtuple)
process.hgcal_cluster_ntuple=cms.Path(process.hgcalClusterNtuple)

process.schedule = cms.Schedule(process.hgcalRawToDigis_step,
                                #process.simhit_analysis,
                                #process.digi_analysis,
                                #process.rechit_analysis,
                                #process.layercluster_scint_analysis,
                                process.hgcal_rechit_ntuple,
                                process.hgcal_cluster_ntuple)

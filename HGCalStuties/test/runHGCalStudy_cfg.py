import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

from myParam_cff import *

#prefix="/uscms_data/d3/jingyu/HGC/layerClusering/CMSSW_10_6_0_pre4_orig/src/"
prefix=""
filename=filename.replace("step1", "step2")

#filename="step2_211_Vtx_pt20.root"
filename="step2_r160_r160_e100_p22_testGit.root"

process = cms.Process("runAnalyzer", eras.Phase2C8_timing_layer_bar)

process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')

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

## HGCal rec hit study
process.load('Validation.HGCalValidation.rechitStudy_cff')

## HGCal scintillator layer cluster study
process.load('Validation.HGCalValidation.hgcalLayerClusterStudyHESintillator_cfi')

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("file:"+prefix+filename)
                            #skipEvents=cms.untracked.uint32(1)
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2)
)

if not prefix=="":
    filename=filename.replace(".root","_orig.root")
    
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("h_"+filename),
                                   closeFileFast = cms.untracked.bool(True)
                                   )


process.simhit_analysis = cms.Path(process.hgcalSimHitStudy)
process.digi_analysis=cms.Path(process.hgcalDigiStudyEE+process.hgcalDigiStudyHEF+process.hgcalDigiStudyHEB)
process.rechit_analysis=cms.Path(process.hgcalRecHitStudyEE+process.hgcalRecHitStudyFH+process.hgcalRecHitStudyBH)
process.layercluster_scint_analysis=cms.Path(process.hgcalLayerClusterStudyHESintillator)

process.schedule = cms.Schedule(process.hgcalRawToDigis_step,
                                process.simhit_analysis,
                                process.digi_analysis,
                                process.rechit_analysis,
                                process.layercluster_scint_analysis)
                                

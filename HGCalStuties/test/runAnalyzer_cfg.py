import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

from myParam_cff import *

#prefix="root://cmseos.fnal.gov//store/user/zhangj/HGC/PGun_10_6_0_pre3/GEN-SIM/"
#prefix="root://cmseos.fnal.gov//store/user/zhangj/HGC/PGun_10_4_0/GEN-SIM/"
#filename="QGun_Pythia8_Pt500_Eta2p2-2p3_1.root"
#filename="RECO_QGun_Pythia8_Pt500_Eta2p2-2p3_1.root"
prefix=""
filename=filename.replace("step1", "step2")

filename="step2_211_Vtx_pt50.root"

#process = cms.Process("runAnalyzer")
process = cms.Process("runAnalyzer", eras.Phase2C8_timing_layer_bar)

process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

## HGCal sim hit study
#process.load('Validation.HGCalValidation.hgcSimHitStudy_cfi')

## HGCal sim hit validation 
#process.load("Validation.HGCalValidation.simhitValidation_cff")
#process.hgcalSimHitValidationEE.Verbosity=0
#process.hgcalSimHitValidationHEB.Verbosity=0
#process.hgcalSimHitValidationHEF.Verbosity=0

## HGCal raw to digis
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.hgcalRawToDigis_step=cms.Path(process.hgcalDigis)

## HGCal rec hit validation
#process.load("Validation.HGCalValidation.rechitValidation_cff")
#process.hgcalRecHitValidationEE.Verbosity=0
#process.hgcalRecHitValidationHEB.Verbosity=0
#process.hgcalRecHitValidationHEF.Verbosity=0
#process.load('Validation.HGCalValidation.rechitStudy_cff')

## HGCal digi validation
#process.load("Validation.HGCalValidation.digiValidation_cff")
#process.hgcalDigiValidationEE.Verbosity=0
#process.hgcalDigiValidationEE.SampleIndx = cms.untracked.int32(2)
#process.hgcalDigiValidationHEB.Verbosity=0
#process.hgcalDigiValidationHEF.Verbosity=0
#process.load('Validation.HGCalValidation.hgcDigiStudy_cfi')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:"+prefix+filename)
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string("h_"+filename),
#                                   closeFileFast = cms.untracked.bool(True)
#                                   )

process.load('Configuration.StandardSequences.DQMSaverAtRunEnd_cff')
process.dqmsave_step = cms.Path(process.DQMSaver)
if nparticles==1:
    process.dqmSaver.workflow="/"+"p"+str(particle[0])+"/"+"r"+str(int(rMin))+"_r"+str(int(rMax))+"_e"+str(int(eMin))+"/v6"
else:
    process.dqmSaver.workflow="/"+"p"+str(particle[0])+"/"+"r"+str(int(rMin))+"_r"+str(int(rMax))+"_e"+str(int(eMin))+"/"+"delta"+str(delta).replace(".","p")+"_v6"

process.dqmSaver.workflow="/p211/Vtx/pt50"

#process.simhit_analysis = cms.Path(process.hgcalSimHitStudy)
#process.hgcalDigiValidation=cms.Path(process.hgcalDigiValidationEE+process.hgcalDigiValidationHEB+process.hgcalDigiValidationHEF)
#process.hgcalSimHitValidation=cms.Path(process.hgcalSimHitValidationEE+process.hgcalSimHitValidationHEB+process.hgcalSimHitValidationHEF)
#process.hgcalRecHitValidation=cms.Path(process.hgcalRecHitValidationEE+process.hgcalRecHitValidationHEB+process.hgcalRecHitValidationHEF)
#process.hgcalDigiValidation=cms.Path(process.hgcalDigiStudyEE+process.hgcalDigiStudyHEF+process.hgcalDigiStudyHEB)
#process.hgcalRecHitValidation=cms.Path(process.hgcalRecHitStudyEE+process.hgcalRecHitStudyFH+process.hgcalRecHitStudyBH)

process.load('Validation.Configuration.hgcalSimValid_cff')
process.hgcalValidation_step=cms.Path(process.hgcalValidation)

process.load("Validation.HGCalValidation.HGCalPostProcessor_cff")
process.hgcalValidatorPostProcessor_step=cms.Path(process.hgcalValidatorPostProcessor)

process.schedule = cms.Schedule(process.hgcalRawToDigis_step,
                                process.hgcalValidation_step,
                                process.hgcalValidatorPostProcessor_step,
                                process.dqmsave_step)
#process.simhit_analysis)
                                #process.simhit_analysis,
                                #process.hgcalDigiValidation,
                                #process.hgcalRecHitValidation)

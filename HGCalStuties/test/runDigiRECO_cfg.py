import FWCore.ParameterSet.Config as cms
import sys

from Configuration.StandardSequences.Eras import eras

from myParam_cff import *

process = cms.Process('RECO',eras.Phase2C8_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D41Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_Fake2_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#filename="step1_eta2p5_eta2p5_pt10_pt20_p211_100.root"

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('file:'+filename),
    inputCommands = cms.untracked.vstring(
        'keep *', 
        'drop *_genParticles_*_*', 
        'drop *_genParticlesForJets_*_*', 
        'drop *_kt4GenJets_*_*', 
        'drop *_kt6GenJets_*_*', 
        'drop *_iterativeCone5GenJets_*_*', 
        'drop *_ak4GenJets_*_*', 
        'drop *_ak7GenJets_*_*', 
        'drop *_ak8GenJets_*_*', 
        'drop *_ak4GenJetsNoNu_*_*', 
        'drop *_ak8GenJetsNoNu_*_*', 
        'drop *_genCandidatesForMET_*_*', 
        'drop *_genParticlesForMETAllVisible_*_*', 
        'drop *_genMetCalo_*_*', 
        'drop *_genMetCaloAndNonPrompt_*_*', 
        'drop *_genMetTrue_*_*', 
        'drop *_genMetIC5GenJs_*_*'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:'+filename.replace("step1","step2")),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# Path and EndPath definitions
## HGCal digitization
process.digitisation_step = cms.Path(process.pdigi_valid)

## HGCal raw to digi
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.hgcalRawToDigis=cms.Path(process.hgcalDigis)

## HGCal Rec Hits
#from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import *
#from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import *

process.load("RecoLocalCalo.Configuration.hgcalLocalReco_cff")
#process.hgcalRec=cms.Path(process.HGCalUncalibRecHit+process.HGCalRecHit+process.hgcalLayerClusters)
process.hgcalRec=cms.Path(process.hgcalLocalRecoSequence)

#from RecoParticleFlow.Configuration.RecoParticleFlow_cff import *
#process.load("RecoParticleFlow.Configuration.RecoParticleFlow_cff")
#process.phase2PF=cms.Path(process._phase2_hgcal_simPFSequence)

#process.load("RecoParticleFlow.PFProducer.simPFProducer_cfi")
#process.load("RecoParticleFlow.PFTracking.hgcalTrackCollection_cfi")
#process.phase2SimPF=cms.Path(process.hgcalTrackCollection+process.simPFProducer)

#process.HGCalUncalibRecHit.HGCEEConfig.fCPerMIP = fCPerMIP_v9
#process.HGCalUncalibRecHit.HGCHEFConfig.fCPerMIP = fCPerMIP_v9
#process.HGCalRecHit.layerWeights = dEdX_weights_v10
#process.HGCalRecHit.thicknessCorrection = cms.vdouble(0.759,0.760,0.773) # v9
#process.HGCalRecHit.thicknessCorrection = cms.vdouble(0.781,0.775,0.769) # v10
#process.HGCalRecHit.HGCEE_fCPerMIP = fCPerMIP_v9
#process.HGCalRecHit.HGCHEF_fCPerMIP = fCPerMIP_v9

#from RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cff import *

## EndPath
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,
                                process.hgcalRawToDigis,
                                process.hgcalRec,
                                #process.phase2SimPF,
                                process.FEVTDEBUGHLToutput_step)

# customisation of the process.

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)

##process.MessageLogger = cms.Service(
##    "MessageLogger",
##    destinations = cms.untracked.vstring(
##        'detailedInfo',
##        'critical'
##    ),
##    detailedInfo = cms.untracked.PSet(
##        threshold = cms.untracked.string('DEBUG')
##    ),
##    debugModules = cms.untracked.vstring(
##        '*'
##    )
##)

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
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
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
    fileName = cms.untracked.string('file:root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/HGC/HDBSCAN/'+filename.replace("step1","step2")),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.input.nbPileupEvents.averageNumber = cms.double(200.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-3)
process.mix.maxBunch = cms.int32(3)
process.mix.input.fileNames = cms.untracked.vstring(['/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/4A596903-2FAB-8744-9A02-E881B18E52BD.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/ACD2DE6B-7789-774A-B0B5-2AF298A1BBD7.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/43BE832E-5381-8741-8081-B591F26167F9.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/5A4087F4-DD16-674E-8205-AD632A51D6E9.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/ACD6E194-20AD-7346-ACF2-0205D33660D0.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/8007FC3E-74FA-6C47-BE67-7BE6695AB42B.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/FFF20BB2-8D84-8A40-A041-F67CB20BDDC1.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/5A1D65A4-4E25-234C-8342-D086B3B549BC.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/AEAD0161-5950-A94E-B06A-05A6C5906ADF.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/37EBE385-1E3D-184E-9503-160030251D2E.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/64C68672-5B53-A64E-8213-D448EDB4CA2B.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/2A8CF3E5-90F0-8F41-952B-91C82949D922.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/50094427-1960-5249-A7A8-56D36A60F2ED.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/7F987C01-7AFA-3745-B3EF-523517857BE3.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/C66B229E-C956-9E46-B51E-1F91F4AC9547.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/A716825E-114C-4C49-A788-2086F5EAB282.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/E89D8696-8B47-0148-9C7C-15FDB5F86958.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/A049F31C-AEF4-1845-B58E-1B590669FE54.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/6556A9A2-5AA8-244E-9BEF-3FC09E0B5622.root', '/store/relval/CMSSW_10_6_0_patch2/RelValMinBias_14TeV/GEN-SIM/106X_upgrade2023_realistic_v3_2023D41noPU-v2/20000/3045EF71-44E0-F14F-8B84-4541E658F3BB.root'])
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

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randHelper = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randHelper.resetSeeds(zMin*rMin)

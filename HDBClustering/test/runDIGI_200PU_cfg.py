import FWCore.ParameterSet.Config as cms
import sys

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9

from myParam_cff import *

process = cms.Process('DIGI',Phase2C9)

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
    fileName = cms.untracked.string('file:root://cmseos.fnal.gov//eos/uscms/store/user/zhangj/HGC/HDBSCAN/MyHDBSCAN/'+filename.replace("step1","step2")),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.input.nbPileupEvents.averageNumber = cms.double(200.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-3)
process.mix.maxBunch = cms.int32(3)
process.mix.input.fileNames = cms.untracked.vstring(['/store/relval/CMSSW_11_0_0_pre7/RelValMinBias_14TeV/GEN-SIM/110X_mcRun4_realistic_v1_2026D41noPU-v1/20000/F9E6A55B-5159-EE42-AD1B-79553F80028C.root', '/store/relval/CMSSW_11_0_0_pre7/RelValMinBias_14TeV/GEN-SIM/110X_mcRun4_realistic_v1_2026D41noPU-v1/20000/AC566B61-92F9-1A47-97B3-6A3042C7BEA2.root', '/store/relval/CMSSW_11_0_0_pre7/RelValMinBias_14TeV/GEN-SIM/110X_mcRun4_realistic_v1_2026D41noPU-v1/20000/A61B6D35-15BC-3144-B9EB-0691A66CB93E.root', '/store/relval/CMSSW_11_0_0_pre7/RelValMinBias_14TeV/GEN-SIM/110X_mcRun4_realistic_v1_2026D41noPU-v1/20000/9CB8EC67-C1C6-F844-9685-E147E6219130.root', '/store/relval/CMSSW_11_0_0_pre7/RelValMinBias_14TeV/GEN-SIM/110X_mcRun4_realistic_v1_2026D41noPU-v1/20000/8A91E40F-729F-EA4C-9269-1359DC83E506.root', '/store/relval/CMSSW_11_0_0_pre7/RelValMinBias_14TeV/GEN-SIM/110X_mcRun4_realistic_v1_2026D41noPU-v1/20000/86877C30-2F19-2B47-8480-ADC31DDCE733.root', '/store/relval/CMSSW_11_0_0_pre7/RelValMinBias_14TeV/GEN-SIM/110X_mcRun4_realistic_v1_2026D41noPU-v1/20000/5E6A8738-4D26-7440-AAE9-46C135A1ECBA.root', '/store/relval/CMSSW_11_0_0_pre7/RelValMinBias_14TeV/GEN-SIM/110X_mcRun4_realistic_v1_2026D41noPU-v1/20000/5B524549-731D-6343-BE0C-279488A80ABB.root', '/store/relval/CMSSW_11_0_0_pre7/RelValMinBias_14TeV/GEN-SIM/110X_mcRun4_realistic_v1_2026D41noPU-v1/20000/334791F4-7A88-A747-84F4-49C9B823CC09.root', '/store/relval/CMSSW_11_0_0_pre7/RelValMinBias_14TeV/GEN-SIM/110X_mcRun4_realistic_v1_2026D41noPU-v1/20000/04016374-5FDF-AA44-B5F0-C35682F4FEFE.root'])
#process.mix.input.fileNames = cms.untracked.vstring('dbs:/RelValMinBias_14TeV/CMSSW_11_0_0_pre7-110X_mcRun4_realistic_v1_2026D41noPU-v1/GEN-SIM')
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# Path and EndPath definitions
## HGCal digitization
process.digitisation_step = cms.Path(process.pdigi_valid)

## HGCal raw to digi
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.hgcalRawToDigis=cms.Path(process.hgcalDigis)

## EndPath
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,
                                #process.hgcalRawToDigis,
                                #process.hgcalRec,
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

import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

#from myParam_cff import *

process = cms.Process('SIM',eras.Phase2C8_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D41Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D41_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

filename="step1_eta2p5_eta2p5_pt20_p15_test.root"

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Single Particle Gun'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:'+filename),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

#process.generator = cms.EDProducer("FlatRandomPtGunProducer",
#    PGunParameters = cms.PSet(
#        PartID = cms.vint32(22),
#        MinPhi = cms.double(-3.14159265359),
#        MaxPhi = cms.double(3.14159265359), ## in radians
#        MinEta = cms.double(1.9),
#        MaxEta = cms.double(1.9),
#        MinPt = cms.double(1.), # in GeV
#        MaxPt = cms.double(20.0)
#    ),
#    Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts
#    AddAntiParticle = cms.bool(False),
#)

#process.generator = cms.EDProducer("CloseByParticleGunProducer",
#    PGunParameters = cms.PSet(
#        PartID = cms.vint32(particle),
#        EnMin = cms.double(eMin),
#        EnMax = cms.double(eMax),
#        RMin = cms.double(rMin),
#        RMax = cms.double(rMax),
#        ZMin = cms.double(zMin),
#        ZMax = cms.double(zMax),
#        Delta = cms.double(delta),
#        Pointing = cms.bool(True),
#        Overlapping = cms.bool(isOverlapping),
#        RandomShoot = cms.bool(False),
#        NParticles = cms.int32(nparticles),
#        MaxEta = cms.double(etaMax),
#        MinEta = cms.double(etaMin),
#        MaxPhi = cms.double(phiMax),
#        MinPhi = cms.double(phiMin),
#    ),
#    Verbosity = cms.untracked.int32(10),
#    psethack = cms.string('single particle random energy'),
#    AddAntiParticle = cms.bool(False),
#    firstRun = cms.untracked.uint32(1)
#)

process.generator = cms.EDProducer("Pythia6PtGun",
    maxEventsToPrint = cms.untracked.int32(5),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    pythiaHepMCVerbosity = cms.untracked.bool(True),    
    PGunParameters = cms.PSet(
        ParticleID = cms.vint32(-15),
        AddAntiParticle = cms.bool(False),
        MinPhi = cms.double(0),
        MaxPhi = cms.double(0),
        MinPt = cms.double(20.0),
        MaxPt = cms.double(20.0001),
        MinEta = cms.double(2.5),
        MaxEta = cms.double(2.5)
    ),
    PythiaParameters = cms.PSet(
        pythiaTauJets = cms.vstring(
            'MDME(89,1)=0      ! no tau->electron',
            'MDME(90,1)=0      ! no tau->muon'
        ),
        pythiaUESettings = cms.vstring(
            'MSTJ(11)=3     ! Choice of the fragmentation function',
            'MSTJ(22)=2     ! Decay those unstable particles',
            'PARJ(71)=10 .  ! for which ctau  10 mm',
            'MSTP(2)=1      ! which order running alphaS',
            'MSTP(33)=0     ! no K factors in hard cross sections',
            'MSTP(51)=7     ! structure function chosen',
            'MSTP(81)=1     ! multiple parton interactions 1 is Pythia default',
            'MSTP(82)=4     ! Defines the multi-parton model',
            'MSTU(21)=1     ! Check on possible errors during program execution',
            'PARP(82)=1.9409   ! pt cutoff for multiparton interactions',
            'PARP(89)=1960. ! sqrts for which PARP82 is set',
            'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter',
            'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter',
            'PARP(90)=0.16  ! Multiple interactions: rescaling power',
            'PARP(67)=2.5    ! amount of initial-state radiation',
            'PARP(85)=1.0  ! gluon prod. mechanism in MI',
            'PARP(86)=1.0  ! gluon prod. mechanism in MI',
            'PARP(62)=1.25   ! ',
            'PARP(64)=0.2    ! ',
            'MSTP(91)=1     !',
            'PARP(91)=2.1   ! kt distribution',
            'PARP(93)=15.0  ! '
        ),
        parameterSets = cms.vstring(
            'pythiaUESettings',
            'pythiaTauJets'
        )
    )
)

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randHelper = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randHelper.resetSeeds(25)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
        getattr(process,path).insert(0, process.ProductionFilterSequence)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

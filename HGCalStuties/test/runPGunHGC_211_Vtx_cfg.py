# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: TTbar_14TeV_TuneCUETP8M1_cfi --python_file runPGun_cfg --conditions auto:phase2_realistic -n 10 --era Phase2C8_timing_layer_bar --eventcontent FEVTDEBUG --relval 9000,100 -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC14TeV --geometry Extended2023D41 --no_exec --fileout file:step1.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SIM',eras.Phase2C8_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D41_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('GGun_14TeV_TuneCUETP8M1_cfi nevts:50'),
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
    fileName = cms.untracked.string('file:step1_211_Vtx_pt20.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

##process.generator = cms.EDFilter("Pythia8PtGun",
##    maxEventsToPrint = cms.untracked.int32(1),
##    pythiaPylistVerbosity = cms.untracked.int32(1),
##    pythiaHepMCVerbosity = cms.untracked.bool(True),
##    PGunParameters = cms.PSet(
##        ParticleID = cms.vint32(1),
##        AddAntiParticle = cms.bool(True),
##        MinPhi = cms.double(-3.14159265359),
##        MaxPhi = cms.double(3.14159265359),
##        MinPt = cms.double(500.0),
##        MaxPt = cms.double(500.0),
##        MinEta = cms.double(1.7),
##        MaxEta = cms.double(1.8)
##        ),
##    PythiaParameters = cms.PSet(parameterSets = cms.vstring())       
##)

##process.generator = cms.EDProducer("CloseByParticleGunProducer",
##    PGunParameters = cms.PSet(
##        PartID = cms.vint32(211),
##        EnMin = cms.double(1000),
##        EnMax = cms.double(1000),
##        RMin = cms.double(109.119980179),
##        RMax = cms.double(109.119980179),
##        #RMin = cms.double(53.0643),
##        #RMax = cms.double(53.0643),
##        ZMin = cms.double(321.05),
##        ZMax = cms.double(321.05),
##        Delta = cms.double(2.5),
##        Pointing = cms.bool(True),
##        Overlapping = cms.bool(False),
##        #MaxEta = cms.double(2.5),
##        #MinEta = cms.double(2.5),
##        MaxEta = cms.double(1.8),
##        MinEta = cms.double(1.8),
##        MaxPhi = cms.double(3.14159265359),
##        MinPhi = cms.double(-3.14159265359),
##    ),
##    Verbosity = cms.untracked.int32(10),
##    psethack = cms.string('single gamma random energy'),
##    AddAntiParticle = cms.bool(True),
##    firstRun = cms.untracked.uint32(1)
##)

pT=20.0

process.generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(211),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359), ## in radians
        MinEta = cms.double(1.7),
        MaxEta = cms.double(1.7),
        MinPt = cms.double(pT), # in GeV
        MaxPt = cms.double(pT)
    ),
    Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts
    AddAntiParticle = cms.bool(False),
)

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randHelper = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randHelper.resetSeeds(1)

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

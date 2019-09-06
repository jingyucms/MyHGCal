# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: TTbar_14TeV_TuneCUETP8M1_cfi --python_file runPGun_cfg --conditions auto:phase2_realistic -n 10 --era Phase2C8_timing_layer_bar --eventcontent FEVTDEBUG --relval 9000,100 -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC14TeV --geometry Extended2023D41 --no_exec --fileout file:step1.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from myParam_cff import *

process = cms.Process('SIM',eras.Phase2C9_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D46Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D46_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(nevt)
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
    fileName = cms.untracked.string('file:'+filename),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'phase2_realistic_T15', '')

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

process.generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(particle[0]),
        MinPhi = cms.double(phiMin),
        MaxPhi = cms.double(phiMax), ## in radians
        MinEta = cms.double(etaMin),
        MaxEta = cms.double(etaMax),
        MinPt = cms.double(ptMin), # in GeV
        MaxPt = cms.double(ptMax)
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

# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --conditions 131X_mcRun4_realistic_v9 -n 2 --era Phase2C17I13M9 --eventcontent FEVTDEBUGHLT -s RAW2DIGI,L1,L1P2GT --datatier GEN-SIM-DIGI-RAW-MINIAOD --fileout file:test.root --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,Configuration/DataProcessing/Utils.addMonitoring,L1Trigger/Configuration/customisePhase2.addHcalTriggerPrimitives,L1Trigger/Configuration/customisePhase2FEVTDEBUGHLT.customisePhase2FEVTDEBUGHLT --geometry Extended2026D95 --nThreads 8 --filein /store/mc/Phase2Spring23DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_L1TFix_Trk1GeV_131X_mcRun4_realistic_v9-v1/50000/005bc30b-cf79-4b3b-9ec1-a80e13072afd.root --mc --inputCommands=keep *, drop l1tPFJets_*_*_*, drop l1tTrackerMuons_l1tTkMuonsGmt_*_* --outputCommands=drop l1tPFJets_*_*_*, drop l1tTrackerMuons_l1tTkMuonsGmt_*_*
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('GMT',Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D95Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.SimPhase2L1GlobalTriggerEmulator_cff')
process.load('L1Trigger.Configuration.Phase2GTMenus.SeedDefinitions.prototypeSeeds')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
'''
process.MessageLogger = cms.Service("MessageLogger",
       destinations   = cms.untracked.vstring(
                                             'detailedInfo'
                                               ,'critical'
                                               ,'cerr'
                    ),
       critical       = cms.untracked.PSet(
                      threshold = cms.untracked.string('ERROR') 
        ),
       detailedInfo   = cms.untracked.PSet(
                      threshold = cms.untracked.string('INFO') 
       ),
       cerr           = cms.untracked.PSet(
                       threshold  = cms.untracked.string('WARNING') 
        )
)
 
 
process.myAnalysisModule = cms.EDAnalyzer('ModuleThatIssuesMessages')
process.p = cms.Path(process.myAnalysisModule)
'''

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
from fixed_DT_samples import *
#import numpy as np
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring(
#'/store/user/bachtis/TEST140X_V1/DYToLL_M-50_TuneCP5_14TeV-pythia8/PHASEII_DY200a/240417_092910/0000/test_175.root'
#np.flip(dy)
dy
),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop l1tPFJets_*_*_*',
        'drop l1tTrackerMuons_l1tTkMuonsGmt_*_*'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:2'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-MINIAOD'),
        filterName = cms.untracked.string('')
    ),
    #fileName = cms.untracked.string('file:test.root'),
    fileName = cms.untracked.string('file:/uscmst1b_scratch/lpc1/3DayLifetime/hancelin/dummy.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '131X_mcRun4_realistic_v9', '')
process.FEVTDEBUGHLToutput.outputCommands.append('drop l1tPFJets_*_*_*')
process.FEVTDEBUGHLToutput.outputCommands.append('drop l1tTrackerMuons_l1tTkMuonsGmt_*_*')
process.FEVTDEBUGHLToutput.outputCommands.append('keep *_*_*_GMT')

# Path and EndPath definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.L1simulation_step = cms.Path(process.SimL1Emulator)

process.L1simulation_step = cms.Path(process.l1tStubsGmt*process.l1tKMTFMuonsGmt*process.l1tFwdMuonsGmt*process.l1tSAMuonsGmt*process.l1tTkMuonsGmt)

#process.Phase2L1GTProducer = cms.Path(process.l1tGTProducerSequence)
#process.Phase2L1GTAlgoBlockProducer = cms.Path(process.l1tGTAlgoBlockProducerSequence)
#process.pDoubleEGEle37_24 = cms.Path(process.DoubleEGEle3724)
#process.pDoubleIsoTkPho22_12 = cms.Path(process.DoubleIsoTkPho2212)
#process.pDoublePuppiTau52_52 = cms.Path(process.DoublePuppiTau5252)
#process.pDoubleTkEle25_12 = cms.Path(process.DoubleTkEle2512)
#process.pDoubleTkMuon15_7 = cms.Path(process.DoubleTkMuon157)
#process.pIsoTkEleEGEle22_12 = cms.Path(process.IsoTkEleEGEle2212)
#process.pPuppiHT400 = cms.Path(process.PuppiHT400)
#process.pPuppiHT450 = cms.Path(process.PuppiHT450)
#process.pPuppiMET200 = cms.Path(process.PuppiMET200)
#process.pQuadJet70_55_40_40 = cms.Path(process.QuadJet70554040)
#process.pSingleEGEle51 = cms.Path(process.SingleEGEle51)
#process.pSingleIsoTkEle28 = cms.Path(process.SingleIsoTkEle28)
#process.pSingleIsoTkPho36 = cms.Path(process.SingleIsoTkPho36)
#process.pSinglePuppiJet230 = cms.Path(process.SinglePuppiJet230)
#process.pSingleTkEle36 = cms.Path(process.SingleTkEle36)
#process.pSingleTkMuon22 = cms.Path(process.SingleTkMuon22)
#process.pTripleTkMuon5_3_3 = cms.Path(process.TripleTkMuon533)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads = 1
process.options.numberOfStreams = 0

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000 

#call to customisation function customise_aging_1000 imported from SLHCUpgradeSimulations.Configuration.aging
process = customise_aging_1000(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customisePhase2
from L1Trigger.Configuration.customisePhase2 import addHcalTriggerPrimitives 

#call to customisation function addHcalTriggerPrimitives imported from L1Trigger.Configuration.customisePhase2
process = addHcalTriggerPrimitives(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customisePhase2FEVTDEBUGHLT
from L1Trigger.Configuration.customisePhase2FEVTDEBUGHLT import customisePhase2FEVTDEBUGHLT 

#call to customisation function customisePhase2FEVTDEBUGHLT imported from L1Trigger.Configuration.customisePhase2FEVTDEBUGHLT
process = customisePhase2FEVTDEBUGHLT(process)

# End of customisation functions


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

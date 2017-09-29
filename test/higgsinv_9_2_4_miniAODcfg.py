import FWCore.ParameterSet.Config as cms
process = cms.Process("MAIN")
import sys

################################################################                                                                                         
# Required additional modules !!CHECK ALL UP TO DATE
# Electron conversion veto calculator:
#   in CMSSW src: git clone https://github.com/ajgilbert/ICAnalysis-ElectronConversionCalculator.git ICAnalysis/ElectronConversionCalculator
# Electron cut based ID:
#   in CMSSW src: git cms-merge-topic ikrav:egm_id_phys14
################################################################                                                                                           

################################################################                                                                                         
# Define some strings                                                                                                                          
################################################################                                                                                           
electronLabel = cms.InputTag("slimmedElectrons")
muonLabel     = cms.InputTag("slimmedMuons")
photonLabel   = cms.InputTag("slimmedPhotons")
tauLabel      = cms.InputTag("slimmedTaus")
vtxLabel      = cms.InputTag("offlineSlimmedPrimaryVertices")

################################################################                                                                                         
# Setup and Read Options                                                                                                                          
################################################################                                                                                           
import FWCore.ParameterSet.VarParsing as parser
opts = parser.VarParsing ('analysis')


opts.register('file',
              0, parser.VarParsing.multiplicity.singleton,
              #'root://xrootd.unl.edu//store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/150/00000/34A57FB8-D819-E611-B0A4-02163E0144EE.root',parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.string, "input file")


opts.register('globalTag', '92X_dataRun2_Prompt_v9', parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.string, "global tag") #to be frequently updated from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Prompt_reconstruction_Global_Tag

opts.register('isData', 1, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Process as data?")

opts.register('release', '92XMINIAOD', parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.string, "Release label")

#opts.register('isNLO', False, parser.VarParsing.multiplicity.singleton,
#    parser.VarParsing.varType.bool, "Store sign of weight?")

opts.register('doHT', 0, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Store HT and number of outgoing partons?")

opts.register('isReHLT', 0, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Process as reHLT sample?")

opts.parseArguments()
infile      = opts.file
if not infile: infile = "file:/tmp/file.root"
isData      = opts.isData
tag         = opts.globalTag
release     = opts.release
if not isData:
  doHT      = opts.doHT
  isReHLT   = opts.isReHLT
else:
  doHT      = 0
  isReHLT   = 0

v_doBXloop  = False  ## False = not do the loop on all the possible BX; True = do it!
##v_doBXloop  = True  ## False = not do the loop on all the possible BX; True = do it!


if not release in ['80XMINIAOD', '92XMINIAOD']:
  print 'Release not recognised, exiting!'
  sys.exit(1)
print 'release     : '+release
print 'isData      : '+str(isData)
print 'globalTag   : '+str(tag)
print 'doHT        : '+str(doHT)
print 'isReHLT     : '+str(isReHLT)

################################################################                                                                                          
# Standard setup                                                                                                                                           
################################################################                                                                                          
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.TFileService = cms.Service("TFileService",
  fileName      = cms.string("EventTree.root"),
  closeFileFast = cms.untracked.bool(True)
)

################################################################                                                                                        
# Message Logging, summary, and number of events                                                                                                          
################################################################                                                                                          
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(100)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options   = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True)
)

################################################################ 
# Input files and global tags
################################################################ 

process.load("CondCore.CondDB.CondDB_cfi")
#from CondCore.DBCommon.CondDBSetup_cfi import *

if not isData:                         
  process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
       ##                      CHECK WHAT TO USE
                              infile
                              #'file:3881F5BF-3825-E611-92B1-0025905D1E0A.root'  ##SOMETHING SIMILAR IF LOCAL
                                )
                              )
else:
  process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
       ##                      CHECK WHAT TO USE
                              infile
                                )
                              )

process.GlobalTag.globaltag = cms.string(tag)

#'root://xrootd-cms.infn.it///store/mc/Phys14DR/VBF_HToInv_M-125_13TeV_powheg-pythia6/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/58D548B0-AB6F-E411-B468-3417EBE34D1A.root')                                                                                         
#load JEC not in global tag
# process.jec = cms.ESSource("PoolDBESSource",
#       DBParameters = cms.PSet(
#         messageLevel = cms.untracked.int32(0)
#         ),
#       timetype = cms.string('runnumber'),
#       toGet = cms.VPSet(
#       cms.PSet(
#             record = cms.string('JetCorrectionsRecord'),
#             tag    = cms.string('JetCorrectorParametersCollection_Summer15_25nsV5_DATA_AK4PFchs'),
#             label  = cms.untracked.string('AK4PFchs')
#             ),
#       #..................................................
#       ## here you add as many jet types as you need
#       ## note that the tag name is specific for the particular sqlite file 
#       ), 
#       connect = cms.string('sqlite:Summer15_25nsV5_DATA.db')
#      # uncomment above tag lines and this comment to use MC JEC
#      # connect = cms.string('sqlite:Summer12_V7_MC.db')
# )
# ## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
# process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

 
import UserCode.ICHiggsTauTau.default_producers_cfi as producers
import UserCode.ICHiggsTauTau.default_selectors_cfi as selectors

#process.load("PhysicsTools.PatAlgos.patSequences_cff")
################################################################
# Re-do PFTau reconstruction
################################################################
#process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")
#switchToPFTauHPS(process)


################################################################
# Declare sequences
################################################################
process.ic80XSequence = cms.Sequence()
process.icMiniAODSequence = cms.Sequence()


################################################################
# Object selection
################################################################
#process.selectedVertices = cms.EDFilter("VertexRefSelector",
#  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
#  cut = cms.string("ndof >= 4 & abs(z) <= 24 & abs(position.Rho) <= 2")
#)

process.selectedElectrons = cms.EDFilter("PATElectronRefSelector",
  src = electronLabel,
  cut = cms.string("pt > 9.5 & abs(eta) < 2.6")
)

#process.selectedElectrons = cms.EDFilter("PATElectronSelector", 
  #src = cms.InputTag("slimmedElectrons"),
  #cut = cms.string("pt > 5 && abs(superCluster.eta)<2.5")
#)

process.selectedPFMuons = cms.EDFilter("PATMuonRefSelector",
  src = muonLabel,
  cut = cms.string("pt > 3.0 & abs(eta) < 2.6")
)

process.selectedPFTaus = cms.EDFilter("PATTauRefSelector",
  src = tauLabel,
  cut = cms.string('pt > 18.0 & abs(eta) < 2.6 & tauID("decayModeFindingNewDMs") > 0.5')
)

process.load('Configuration.StandardSequences.Services_cff')

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedElectrons = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                                                      engineName = cms.untracked.string('TRandom3'),
                                                                                    ),
                                                   calibratedPhotons = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                                                    engineName = cms.untracked.string('TRandom3'),
                                                                                  ),
                                                  )
process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
process.load('EgammaAnalysis.ElectronTools.calibratedPhotonsRun2_cfi')


#process.selectedPF = cms.EDFilter("PATPackedCandidateSelector",#
#  src = cms.InputTag("packedPFCandidates"),
#  cut = cms.string("fromPV > 1&&charge!=0&&pt>0.9")
#)

process.icSelectionSequence = cms.Sequence(
  #process.selectedPF+
  #process.selectedVertices+
  process.selectedElectrons+
  process.selectedPFMuons+
  process.selectedPFTaus
  #+process.calibratedElectrons
  #+process.calibratedPhotons 
)


################################################################                                                                                            
# PF sequence for lepton isolation                                             
#!!andrew custom code to redo iso, look into taking default miniaod values by making value map
#used as input for CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff
################################################################                                                                                          

process.load("CommonTools.ParticleFlow.pfParticleSelection_cff")

if release in ['80XMINIAOD', '92XMINIAOD']:
  process.pfPileUp = cms.EDFilter("CandPtrSelector",
                                  src = cms.InputTag("packedPFCandidates"),
                                  cut = cms.string("fromPV <= 1")
                                  )
  process.pfNoPileUp = cms.EDFilter("CandPtrSelector",
                                    src = cms.InputTag("packedPFCandidates"),
                                    cut = cms.string("fromPV > 1")
                                    )
  process.pfAllNeutralHadrons = cms.EDFilter("CandPtrSelector",
                                             src = cms.InputTag("pfNoPileUp"),
                                             cut = cms.string("abs(pdgId) = 111 | abs(pdgId) = 130 | " \
                                                              "abs(pdgId) = 310 | abs(pdgId) = 2112")
                                             )
  process.pfAllChargedHadrons= cms.EDFilter("CandPtrSelector",
                                            src = cms.InputTag("pfNoPileUp"),
                                            cut = cms.string("abs(pdgId) = 211 | abs(pdgId) = 321 | " \
                                                             "abs(pdgId) = 999211 | abs(pdgId) = 2212")
                                            )
  process.pfAllPhotons= cms.EDFilter("CandPtrSelector",
                                     src = cms.InputTag("pfNoPileUp"),
                                     cut = cms.string("abs(pdgId) = 22")
                                     )
  process.pfAllChargedParticles= cms.EDFilter("CandPtrSelector",
                                              src = cms.InputTag("pfNoPileUp"),
                                              cut = cms.string("abs(pdgId) = 211 | abs(pdgId) = 321 | " \
                                                               "abs(pdgId) = 999211 | abs(pdgId) = 2212 | " \
                                                               "abs(pdgId) = 11 | abs(pdgId) = 13")
                                              )
  process.pfPileUpAllChargedParticles= cms.EDFilter("CandPtrSelector",
                                                    src = cms.InputTag("pfPileUp"),
                                                    cut = cms.string("abs(pdgId) = 211 | abs(pdgId) = 321 | " \
                                                                     "abs(pdgId) = 999211 | abs(pdgId) = 2212 | " \
                                                                     "abs(pdgId) = 11 | abs(pdgId) = 13")
                                                    )
  process.pfAllNeutralHadronsAndPhotons = cms.EDFilter("CandPtrSelector",
                                                       src = cms.InputTag("pfNoPileUp"),
                                                       cut = cms.string("abs(pdgId) = 111 | abs(pdgId) = 130 | " \
                                                                        "abs(pdgId) = 310 | abs(pdgId) = 2112 | abs(pdgId) = 22")
                                                       )

  process.pfParticleSelectionSequence = cms.Sequence(
    process.pfPileUp+
    process.pfNoPileUp+
    process.pfAllNeutralHadrons+
    process.pfAllChargedHadrons+
    process.pfAllPhotons+
    process.pfAllChargedParticles+
    process.pfPileUpAllChargedParticles+
    process.pfAllNeutralHadronsAndPhotons
    )

################################################################                                                                                          
# Vertices                                                                                                                                                 
################################################################                                                                                           
process.icVertexProducer = producers.icVertexProducer.clone(
  branch  = cms.string("vertices"),
  input = cms.InputTag("offlineSlimmedPrimaryVertices"),
  firstVertexOnly = cms.bool(True)
)

process.icGenVertexProducer = producers.icGenVertexProducer.clone(
  input = cms.InputTag("prunedGenParticles")
)

process.icVertexSequence = cms.Sequence(
  process.icVertexProducer+
  process.icGenVertexProducer
)

if isData :
  process.icVertexSequence.remove(process.icGenVertexProducer)


################################################################
# Electrons
################################################################
process.icElectronSequence = cms.Sequence()

#ICElectronConversionCalculator NOT final, but at least have a running version for now
process.icElectronConversionCalculator = cms.EDProducer('ICElectronConversionCalculator',
    input       = electronLabel,
    beamspot    = cms.InputTag("offlineBeamSpot"),
    conversions = cms.InputTag("reducedEgamma:reducedConversions")
    )

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
process.electronMVAValueMapProducer.src = electronLabel
process.electronMVAValueMapProducer.srcMiniAOD = electronLabel

process.icElectronSequence+=cms.Sequence(
#   process.egmGsfElectronIDSequence+
   process.electronMVAValueMapProducer
   )

#Electron PF iso sequence:
process.load("CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff")
#  process.elPFIsoValueCharged03PFIdPFIso    = process.icElectronCharged.clone()
#process.elPFIsoValueCharged03PFIdPFIso    = process.elPFIsoValueCharged03PFId.clone()
process.elPFIsoValueChargedAll03PFIdPFIso = process.elPFIsoValueChargedAll03PFId.clone()
#process.elPFIsoValueGamma03PFIdPFIso      = process.elPFIsoValueGamma03PFId.clone()
#process.elPFIsoValueNeutral03PFIdPFIso    = process.elPFIsoValueNeutral03PFId.clone()
#process.elPFIsoValuePU03PFIdPFIso         = process.elPFIsoValuePU03PFId.clone()
#  process.elPFIsoDepositCharged.src     = electronLabel
process.elPFIsoDepositChargedAll.src  = electronLabel
#  process.elPFIsoDepositNeutral.src     = electronLabel
#  process.elPFIsoDepositGamma.src       = electronLabel
#  process.elPFIsoDepositPU.src          = electronLabel
process.elPFIsoValueCharged03PFIdPFIso = cms.EDProducer('ICRecoElectronIsolation',
  input        = electronLabel,
  deltaR       = cms.double(0.3),
  iso_type = cms.string("charged_iso") 
)    
process.elPFIsoValueGamma03PFIdPFIso = cms.EDProducer('ICRecoElectronIsolation',
  input        = electronLabel,
  deltaR       = cms.double(0.3),
  iso_type = cms.string("photon_iso") 
)    
process.elPFIsoValueNeutral03PFIdPFIso = cms.EDProducer('ICRecoElectronIsolation',
  input        = electronLabel,
  deltaR       = cms.double(0.3),
  iso_type = cms.string("neutral_iso") 
)    
process.elPFIsoValuePU03PFIdPFIso = cms.EDProducer('ICRecoElectronIsolation',
  input        = electronLabel,
  deltaR       = cms.double(0.3),
  iso_type = cms.string("pu_iso") 
)    

process.electronPFIsolationValuesSequence = cms.Sequence(
      process.elPFIsoValueCharged03PFIdPFIso+
      process.elPFIsoValueChargedAll03PFIdPFIso+
      process.elPFIsoValueGamma03PFIdPFIso+
      process.elPFIsoValueNeutral03PFIdPFIso+
      process.elPFIsoValuePU03PFIdPFIso
      )

if release in ['80XMINIAOD', '92XMINIAOD']:
  process.elPFIsoDepositGamma.ExtractorPSet.ComponentName = cms.string("CandViewExtractor")

process.icElectronSequence += cms.Sequence(
#      process.electronPFIsolationDepositsSequence+
    process.elPFIsoDepositChargedAll+
    process.electronPFIsolationValuesSequence
    )

process.elPFIsoValueChargedAll03PFIdPFIso.deposits[0].vetos = (
    cms.vstring('EcalEndcaps:ConeVeto(0.015)','EcalBarrel:ConeVeto(0.01)'))

if release in ['80XMINIAOD', '92XMINIAOD']:
  process.elPFIsoValueCharged04PFIdPFIso = cms.EDProducer('ICElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.4),
    iso_type = cms.string("charged_iso") 
  )    
  process.elPFIsoValueGamma04PFIdPFIso = cms.EDProducer('ICElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.4),
    iso_type = cms.string("photon_iso") 
  )    
  process.elPFIsoValueNeutral04PFIdPFIso = cms.EDProducer('ICElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.4),
    iso_type = cms.string("neutral_iso") 
  )    
  process.elPFIsoValuePU04PFIdPFIso = cms.EDProducer('ICElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.4),
    iso_type = cms.string("pu_iso") 
  )    

  process.elEcalPFClusterIso = cms.EDProducer('ICElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.3), 
    iso_type     = cms.string("ecal_pf_cluster_iso")
  )

  process.elHcalPFClusterIso = cms.EDProducer('ICElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.3), 
    iso_type     = cms.string("hcal_pf_cluster_iso")
  )

process.elPFIsoValueChargedAll04PFIdPFIso = process.elPFIsoValueChargedAll04PFId.clone()
process.elPFIsoValueChargedAll04PFIdPFIso.deposits[0].vetos = (
  cms.vstring('EcalEndcaps:ConeVeto(0.015)','EcalBarrel:ConeVeto(0.01)'))

process.electronPFIsolationValuesSequence +=cms.Sequence(
  process.elPFIsoValueCharged04PFIdPFIso+
  process.elPFIsoValueGamma04PFIdPFIso+
  process.elPFIsoValuePU04PFIdPFIso+
  process.elPFIsoValueNeutral04PFIdPFIso+
  process.elPFIsoValueChargedAll04PFIdPFIso+
  process.elEcalPFClusterIso+
  process.elHcalPFClusterIso
)


process.icElectronProducer = producers.icElectronProducer.clone(
  branch                    = cms.string("electrons"),
  input                     = cms.InputTag("selectedElectrons"),
  includeConversionMatches  = cms.bool(True),
  inputConversionMatches    = cms.InputTag("icElectronConversionCalculator"),
  includeVertexIP           = cms.bool(True),
  inputVertices             = vtxLabel,
  includeBeamspotIP         = cms.bool(True),
  inputBeamspot             = cms.InputTag("offlineBeamSpot"),
  includeFloats = cms.PSet(
     mvaNonTrigSpring15    = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
     mvaTrigSpring15       = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values")
  ),
  includeClusterIso        = cms.bool(True),
  includePFIso03           = cms.bool(True),
  includePFIso04           = cms.bool(True)
)


#if release in ['80XMINIAOD', '92XMINIAOD']:
#  process.icElectronProducer.includeClusterIso = cms.bool(False)

process.icElectronSequence += cms.Sequence(
  process.icElectronConversionCalculator+
  process.icElectronProducer
)


################################################################ 
# Muons
################################################################

process.icMuonSequence = cms.Sequence()

#redo isolation calculation
process.load("CommonTools.ParticleFlow.Isolation.pfMuonIsolation_cff")
process.load("CommonTools.ParticleFlow.deltaBetaWeights_cff")


process.muPFIsoDepositCharged.src                         = muonLabel
process.muPFIsoDepositChargedAll.src                      = muonLabel
process.muPFIsoDepositNeutral.src                         = muonLabel
process.muPFIsoDepositNeutral.ExtractorPSet.inputCandView = cms.InputTag("pfWeightedNeutralHadrons")
process.muPFIsoDepositGamma.src                           = muonLabel
process.muPFIsoDepositGamma.ExtractorPSet.inputCandView   = cms.InputTag("pfWeightedPhotons")
process.muPFIsoDepositPU.src                              = muonLabel



process.muPFIsoValueChargedAll03PFIso  = process.muPFIsoValueChargedAll03.clone()
process.muPFIsoValueNeutral04PFWeights = process.muPFIsoValueNeutral04.clone()
process.muPFIsoValueGamma04PFWeights   = process.muPFIsoValueGamma04.clone()
process.muPFIsoValueNeutral03PFWeights = process.muPFIsoValueNeutral03.clone()
process.muPFIsoValueGamma03PFWeights   = process.muPFIsoValueGamma03.clone()
process.muPFIsoValueCharged03PFIso     = cms.EDProducer('ICMuonIsolation',
  input        = muonLabel,
  deltaR       = cms.double(0.3),
  iso_type = cms.string("charged_iso") 
)    
process.muPFIsoValueGamma03PFIso       = cms.EDProducer('ICMuonIsolation',
  input        = muonLabel,
  deltaR       = cms.double(0.3),
  iso_type = cms.string("photon_iso") 
)    
process.muPFIsoValueNeutral03PFIso     = cms.EDProducer('ICMuonIsolation',
  input        = muonLabel,
  deltaR       = cms.double(0.3),
  iso_type = cms.string("neutral_iso") 
)    
process.muPFIsoValuePU03PFIso          = cms.EDProducer('ICMuonIsolation',
  input        = muonLabel,
  deltaR       = cms.double(0.3),
  iso_type = cms.string("pu_iso") 
)

process.muonPFIsolationValuesSequence  = cms.Sequence(
   
   process.muPFIsoValueCharged03PFIso+
   process.muPFIsoValueChargedAll03PFIso+
   process.muPFIsoValueGamma03PFIso+
   process.muPFIsoValueNeutral03PFIso+
   process.muPFIsoValuePU03PFIso+
   process.muPFIsoValueNeutral04PFWeights+
   process.muPFIsoValueGamma04PFWeights+
   process.muPFIsoValueNeutral03PFWeights+
   process.muPFIsoValueGamma03PFWeights
   )


process.muPFIsoValueCharged04PFIso     = cms.EDProducer('ICMuonIsolation',
  input        = muonLabel,
  deltaR       = cms.double(0.4),
  iso_type = cms.string("charged_iso") 
)    
process.muPFIsoValueGamma04PFIso       = cms.EDProducer('ICMuonIsolation',
  input        = muonLabel,
  deltaR       = cms.double(0.4),
  iso_type = cms.string("photon_iso") 
)    
process.muPFIsoValueNeutral04PFIso     = cms.EDProducer('ICMuonIsolation',
  input        = muonLabel,
  deltaR       = cms.double(0.4),
  iso_type = cms.string("neutral_iso") 
)    
process.muPFIsoValuePU04PFIso          = cms.EDProducer('ICMuonIsolation',
  input        = muonLabel,
  deltaR       = cms.double(0.4),
  iso_type = cms.string("pu_iso") 
)    
process.muPFIsoValueChargedAll04PFIso  = process.muPFIsoValueChargedAll04.clone()
process.muonPFIsolationValuesSequence +=cms.Sequence(
  process.muPFIsoValueCharged04PFIso+
  process.muPFIsoValueGamma04PFIso+
  process.muPFIsoValuePU04PFIso+
  process.muPFIsoValueNeutral04PFIso+
  process.muPFIsoValueChargedAll04PFIso
)


process.icMuonSequence += cms.Sequence(
    process.pfDeltaBetaWeightingSequence+
    process.muPFIsoDepositChargedAll+
    process.muPFIsoDepositNeutral+
    process.muPFIsoDepositGamma+
    process.muonPFIsolationValuesSequence
    )

#make and store ic muon object
process.icMuonProducer = producers.icMuonProducer.clone(
  branch                    = cms.string("muons"),
  input                     = cms.InputTag("selectedPFMuons"),
  isPF                      = cms.bool(False),
  includeVertexIP           = cms.bool(True),
  inputVertices             = vtxLabel,
  includeBeamspotIP         = cms.bool(True),
  inputBeamspot             = cms.InputTag("offlineBeamSpot"),
  includeDoubles = cms.PSet(
   neutral_pfweighted_iso_03 = cms.InputTag("muPFIsoValueNeutral03PFWeights"),
   neutral_pfweighted_iso_04 = cms.InputTag("muPFIsoValueNeutral04PFWeights"),
   gamma_pfweighted_iso_03 = cms.InputTag("muPFIsoValueGamma03PFWeights"),
   gamma_pfweighted_iso_04 = cms.InputTag("muPFIsoValueGamma04PFWeights")
  ),
  requestTracks           = cms.bool(True),
  includePFIso03           = cms.bool(True),
  includePFIso04           = cms.bool(True),
)

if release in ['80XMINIAOD', '92XMINIAOD']:
  process.icMuonProducer.isPF = cms.bool(False)

process.icMuonSequence += cms.Sequence(
  process.icMuonProducer
)

################################################################                                                                                            
# Taus                                                                                                                                                   
################################################################                                                                                            
import UserCode.ICHiggsTauTau.tau_discriminators_cfi as tauIDs

process.icTauProducer = producers.icTauProducer.clone(
  input                   = cms.InputTag("selectedPFTaus"),
  inputVertices           = vtxLabel,
  includeVertexIP         = cms.bool(True),
  requestTracks           = cms.bool(True),
  includeTotalCharged     = cms.bool(False),
  totalChargedLabel       = cms.string('totalCharged'),
  tauIDs = tauIDs.dynamicStripIds

)

if release in ['80XMINIAOD', '92XMINIAOD']:
  process.icTauProducer = cms.EDProducer("ICPFTauFromPatProducer",
    branch                  = cms.string("taus"),
    input                   = cms.InputTag("selectedPFTaus"),
    inputVertices           = vtxLabel,
    includeVertexIP         = cms.bool(True),
    requestTracks           = cms.bool(False),
    includeTotalCharged     = cms.bool(False),
    totalChargedLabel       = cms.string('totalCharged'),
    requestPFCandidates   = cms.bool(False),
    inputPFCandidates     = cms.InputTag("pfCandidates"),
    isSlimmed             = cms.bool(False),
    tauIDs = cms.PSet()
)

process.icTauSequence = cms.Sequence(
  process.icTauProducer
)

##make and store taus
#process.icTauProducer = cms.EDProducer("ICPFTauFromPatProducer",
#  branch                  = cms.string("taus"),
#  input                   = cms.InputTag("selectedPFTaus"),
#  inputVertices           = vtxLabel,
#  includeVertexIP         = cms.bool(True),
#  requestTracks           = cms.bool(False),
#  tauIDs = cms.PSet()
#)
#
#process.icTauSequence = cms.Sequence(
#  process.icTauProducer
#)

################################################################
# L1 Taus
################################################################

if isData or isReHLT:
  process.icL1TauProducer = cms.EDProducer("ICL1TObjectProducer<l1t::Tau>",
    branch = cms.string("l1tTau"),
    input  = cms.InputTag("caloStage2Digis","Tau","RECO"),
    doBXloop = cms.bool(v_doBXloop)
  )
  if isReHLT:
    process.icL1TauProducer.input = cms.InputTag("caloStage2Digis","Tau","HLT2")


##################################################################                                                                                          
#  Photons
#################################################################
process.icPhotonSequence = cms.Sequence()

#Calculate photon conversion veto
process.icPhotonElectronConversionVetoCalculator = cms.EDProducer("ICPhotonElectronConversionVetoFromPatCalculator",
                                                        input = photonLabel
)

#produce photons using 0.3 radius isolation taken straight from slimmedPhotons
process.icPhotonProducer = producers.icPhotonProducer.clone(
                                          input = photonLabel,
                                          branch = cms.string("photons"),
                                          includeElectronVeto=cms.bool(True),
                                          inputElectronVeto=cms.InputTag("icPhotonElectronConversionVetoCalculator"),
                                          includeHadTowOverEm=cms.bool(True),
                                          includePFIso03=cms.bool(True),
                                          includeIsoFromPat=cms.bool(True) 
)

process.icPhotonSequence += cms.Sequence(
  process.icPhotonElectronConversionVetoCalculator+
  process.icPhotonProducer
)

################################################################## 
#  Jets and rho
#################################################################

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets

#Make ak4 CHS jets:#rebuild ak4 chs jets as in  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Jets

process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

  #Reapply JECs:
if not isData:
  updateJetCollection(
    process,
    jetSource = cms.InputTag("slimmedJets"),
    labelName = "UpdatedJEC",
    jetCorrections = ("AK4PFchs", cms.vstring(['L1FastJet','L2Relative','L3Absolute']), 'None')
    )
else:
  updateJetCollection(
    process,
    jetSource = cms.InputTag("slimmedJets"),
    labelName = "UpdatedJEC",
    jetCorrections = ("AK4PFchs", cms.vstring(['L1FastJet','L2Relative','L3Absolute','L2L3Residual']), 'None')
    )
  

#Also make non-chs jets:
#process.ak4PFJets = ak4PFJets.clone(src='packedPFCandidates',doAreaFastjet=True)
 
#Get slimmedJets direct from miniAOD
process.selectedSlimmedJetsAK4 = cms.EDFilter("PATJetRefSelector",
                                              #src = cms.InputTag("slimmedJets"),
                                              src = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
                                              cut = cms.string("pt > 15")
                                              )

#Get slimmedJetsPuppi direct from miniAOD
process.selectedSlimmedJetsPuppiAK4 = cms.EDFilter("PATJetRefSelector",
                                              src = cms.InputTag("slimmedJetsPuppi"),
                                              cut = cms.string("pt > 15")
                                              )

# select hadrons and partons for the jet flavour
#process.selectedHadronsAndPartons = cms.EDProducer('HadronAndPartonSelector',
#    src = cms.InputTag("generator"),
#    particles = cms.InputTag("prunedGenParticles"),
#    partonMode = cms.string("Auto")
#)

#process.pfchsJetFlavourAssociation = cms.EDProducer("JetFlavourClustering",
#    jets                     = cms.InputTag("ak4PFJetsCHS"),
#    bHadrons                 = cms.InputTag("selectedHadronsAndPartons","bHadrons"),
#    cHadrons                 = cms.InputTag("selectedHadronsAndPartons","cHadrons"),
#    partons                  = cms.InputTag("selectedHadronsAndPartons","algorithmicPartons"),
#    leptons                  = cms.InputTag("selectedHadronsAndPartons","leptons"),
#    jetAlgorithm             = cms.string("AntiKt"),
#    rParam                   = cms.double(0.4),
#    ghostRescaling           = cms.double(1e-18),
#    hadronFlavourHasPriority = cms.bool(False)
#)

#process.icPFchsJetFlavourCalculator = cms.EDProducer('ICJetFlavourCalculator',
#    input       = cms.InputTag("ak4PFJetsCHS"),
#    flavourMap  = cms.InputTag("pfchsJetFlavourAssociation")
#  )


# Jet energy corrections for reclustered
#CHS

 # Jet energy corrections
 # ----------------------
process.ak4PFL1FastjetCHS = cms.EDProducer("L1FastjetCorrectorProducer",
#    srcRho = cms.InputTag("kt6PFJets", "rho"),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll"),
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1FastJet')
)
process.ak4PFL2RelativeCHS = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2Relative')
)
process.ak4PFL3AbsoluteCHS = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L3Absolute')
)
process.ak4PFResidualCHS = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2L3Residual')
)

#Corrections applied to miniaod slimmedJets
pfchsJECS = cms.PSet(
  L1FastJet  = cms.string("ak4PFL1FastjetCHS"),
  L2Relative = cms.string("ak4PFL2RelativeCHS"),
  L3Absolute = cms.string("ak4PFL3AbsoluteCHS")
)
if isData: pfchsJECS = cms.PSet(
  L1FastJet  = cms.string("ak4PFL1FastjetCHS"),
  L2Relative = cms.string("ak4PFL2RelativeCHS"),
  L3Absolute = cms.string("ak4PFL3AbsoluteCHS"),
  L2L3Residual = cms.string("ak4PFResidualCHS")
)


# b-tagging for reclustered
process.load("RecoJets.JetAssociationProducers.ak4JTA_cff")
from RecoJets.JetAssociationProducers.ak4JTA_cff import ak4JetTracksAssociatorAtVertex
process.load("RecoBTag.Configuration.RecoBTag_cff")
import RecoBTag.Configuration.RecoBTag_cff as btag
#process.jetTracksAssociatorAtVertexAK4PFCHS = ak4JetTracksAssociatorAtVertex.clone(
#  jets   = cms.InputTag("ak4PFJetsCHS"),
#  tracks = cms.InputTag("unpackedTracksAndVertices"),
#  pvSrc  = cms.InputTag("unpackedTracksAndVertices")
#  )

#process.impactParameterTagInfosAK4PFCHS = btag.pfImpactParameterTagInfos.clone(
#  #jetTracks = cms.InputTag('jetTracksAssociatorAtVertexAK4PF'),
#  computeGhostTrack=cms.bool(False),
#  candidates=cms.InputTag("packedPFCandidates"),
#  primaryVertex = cms.InputTag("unpackedTracksAndVertices")
#  )

#process.secondaryVertexTagInfosAK4PFCHS = btag.pfSecondaryVertexTagInfos.clone(
#  trackIPTagInfos = cms.InputTag('impactParameterTagInfosAK4PFCHS')
#  )

##Btag discriminators need to be properly updated

#process.simpleSecondaryVertexHighEffBJetTagsAK4PFCHS = btag.pfSimpleSecondaryVertexHighEffBJetTags.clone (
#  tagInfos = cms.VInputTag('secondaryVertexTagInfosAK4PFCHS')
#  )
#process.simpleSecondaryVertexHighPurBJetTagsAK4PFCHS = btag.pfSimpleSecondaryVertexHighPurBJetTags.clone (
#  tagInfos = cms.VInputTag('secondaryVertexTagInfosAK4PFCHS')
#  )
#process.combinedSecondaryVertexBJetTagsAK4PFCHS      = btag.pfCombinedSecondaryVertexV2BJetTags.clone (
#  tagInfos = cms.VInputTag('impactParameterTagInfosAK4PFCHS', 'secondaryVertexTagInfosAK4PFCHS')
#  )

#process.btaggingSequenceAK4PFCHS = cms.Sequence(
#  process.jetTracksAssociatorAtVertexAK4PFCHS
#  +process.impactParameterTagInfosAK4PFCHS
#  +process.secondaryVertexTagInfosAK4PFCHS
#  +process.simpleSecondaryVertexHighEffBJetTagsAK4PFCHS
#  +process.simpleSecondaryVertexHighPurBJetTagsAK4PFCHS
#  +process.combinedSecondaryVertexBJetTagsAK4PFCHS
# )

#process.btaggingSequenceAK4PFCHS += cms.Sequence(
#  process.pfImpactParameterTagInfos+
#  process.pfSecondaryVertexTagInfos+
#  process.pfCombinedSecondaryVertexV2BJetTags
#)

if release in ['80XMINIAOD', '92XMINIAOD']:
  process.pfImpactParameterTagInfos.primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices")
  process.pfImpactParameterTagInfos.candidates = cms.InputTag("packedPFCandidates")




# Pileup ID
# ---------

process.load("RecoJets.JetProducers.PileupJetID_cfi")

if release in ['80XMINIAOD', '92XMINIAOD']:
  process.pileupJetIdCalculator.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
  process.pileupJetIdCalculator.jets     = cms.InputTag("slimmedJets")
#  process.pileupJetIdCalculator.rho      = cms.InputTag("fixedGridRhoFastjetAll")
  process.pileupJetIdEvaluator.jets      = cms.InputTag("slimmedJets")
#  process.pileupJetIdEvaluator.rho       = cms.InputTag("fixedGridRhoFastjetAll")


process.pileupJetIdUpdated = process.pileupJetId.clone(
  jets=cms.InputTag("slimmedJets"),
  inputIsCorrected=True,
  applyJec=True,
  vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
  )

#print process.pileupJetIdUpdated.dumpConfig()
process.updatedPatJetsUpdatedJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']

#process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
#process.patJetCorrFactorsReapplyJEC = process.updatedPatJetCorrFactors.clone(
  #src = cms.InputTag("slimmedJets"),
  #levels = ['L1FastJet', 'L2Relative', 'L3Absolute'] )
#process.updatedJets = process.updatedPatJets.clone(
  #jetSource = cms.InputTag("slimmedJets"),
  #jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
  #)
#process.updatedJets.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']

#Produce and store jets taken straight from miniAOD
process.icPFJetProducerFromPat = producers.icPFJetFromPatProducer.clone(
  branch                    = cms.string("pfJetsPFlow"),
  input                     = cms.InputTag("selectedSlimmedJetsAK4"),
  #input                     = cms.InputTag("slimmedJets"),
  srcConfig = cms.PSet(
    isSlimmed               = cms.bool(True),
    slimmedPileupIDLabel    = cms.string('pileupJetIdUpdated:fullDiscriminant'),
    #slimmedPileupIDLabel    = cms.string('pileupJetId:fullDiscriminant'),
    includeJetFlavour       = cms.bool(True),
    includeJECs             = cms.bool(True),
    inputSVInfo             = cms.InputTag(""),
    requestSVInfo           = cms.bool(False)
    ),
  destConfig = cms.PSet(
    includePileupID         = cms.bool(True),
    inputPileupID           = cms.InputTag("puJetMva", "fullDiscriminant"),
    includeTrackBasedVars   = cms.bool(False),
    inputTracks             = cms.InputTag("unpackedTracksAndVertices"),
    inputVertices           = cms.InputTag("unpackedTracksAndVertices"),
    requestTracks           = cms.bool(False)
    )
  )

process.icPFJetProducerFromPatPuppi = producers.icPFJetFromPatProducer.clone(
  branch                    = cms.string("ak4SlimmedJetsPuppi"),
  input                     = cms.InputTag("selectedSlimmedJetsPuppiAK4"),
  srcConfig = cms.PSet(
    isSlimmed               = cms.bool(True),
    slimmedPileupIDLabel    = cms.string('pileupJetId:fullDiscriminant'),
    includeJetFlavour       = cms.bool(True),
    includeJECs             = cms.bool(True),
    inputSVInfo             = cms.InputTag(""),
    requestSVInfo           = cms.bool(False)
    ),
  destConfig = cms.PSet(
    includePileupID         = cms.bool(False),
    inputPileupID           = cms.InputTag("puJetMva", "fullDiscriminant"),
    includeTrackBasedVars   = cms.bool(False),
    inputTracks             = cms.InputTag("unpackedTracksAndVertices"),#!!check this and line below
    inputVertices           = cms.InputTag("unpackedTracksAndVertices"),#!!
    requestTracks           = cms.bool(False)
    )
  )

process.icPFJetSequence = cms.Sequence()
process.icPFJetSequence += cms.Sequence(
  process.pileupJetIdUpdated+
  process.patJetCorrFactorsUpdatedJEC+
  process.updatedPatJetsUpdatedJEC+
  #process.patJetCorrFactorsReapplyJEC+
  #process.updatedJets+
  process.selectedUpdatedPatJetsUpdatedJEC+
  process.selectedSlimmedJetsAK4+
# WHY?  process.unpackedTracksAndVertices+
  #process.ak4PFL1FastjetCHS+
  #process.ak4PFL2RelativeCHS+
  #process.ak4PFL3AbsoluteCHS+  
  #process.ak4PFResidualCHS+
  #process.ak4PFJetsCHS+
  #process.btaggingSequenceAK4PFCHS+
  #process.puJetMvaCHS
  #process.patJetsReapplyJEC+
  #process.selectedReJECSlimmedJetsAK4+
  #process.icPFJetProducer+ #Not from slimmed jets!
  process.icPFJetProducerFromPat+
  process.pileupJetIdCalculator+
  process.pileupJetIdEvaluator+
  process.selectedSlimmedJetsPuppiAK4+
  process.icPFJetProducerFromPatPuppi
)

##################################################################                                                                                          
#  MET                                                                                                                                                    
#!!currently type 1 met only, no type0pc+type1 recipe and also haven't implemented mva met
#################################################################





# ################################################################
# # MVA MET and PF MET
# ################################################################

#Make uncorrected pfMet from packedPF
from RecoMET.METProducers.PFMET_cfi import pfMet

process.icPfMet = pfMet.clone(src = "packedPFCandidates")
process.icPfMet.calculateSignificance = False # this can't be easily implemented on packed PF candidates at the moment

process.icuncorrectedPfMetProducer = producers.icMetProducer.clone(
                                                    input = cms.InputTag("icPfMet"),
                                                    branch = cms.string("uncorrectedpfMet"),
                                                    includeCustomID = cms.bool(False),
                                                    inputCustomID = cms.InputTag(""),
                                                    )


#setup met significance calculator for
process.load("RecoMET/METProducers.METSignificance_cfi")
process.load("RecoMET/METProducers.METSignificanceParams_cfi")

#!!clone met sig producer with different name and input for met uncertainties

#get genmet
process.icGenMetProducer = producers.icMetFromPatProducer.clone(
  input = cms.InputTag("slimmedMETs"),
  branch = cms.string("genMetCollection"),
  includeCustomID = cms.bool(False),
  inputCustomID = cms.InputTag(""),
  doGenMet = cms.bool(True)
  )

# Use the newly corrected jets we produced above
#process.basicJetsForMet.src = cms.InputTag('updatedPatJetsUpdatedJEC')

#if not isData:
#    process.patPFMet.addGenMET = cms.bool(False)

process.ictype1PfMetProducer = producers.icMetFromPatProducer.clone(
  input = cms.InputTag("slimmedMETs"),
  #input = cms.InputTag("patPFMetT1"),
  branch = cms.string("pfMetType1Collection"),
  includeCustomID = cms.bool(False),
  inputCustomID = cms.InputTag(""),
  includeExternalMetsigMethod2 = cms.bool(True),
  includeMetCorrections = cms.bool(True),
  includeMetUncertainties = cms.bool(True)
  )

#process.ictype1PfMetMuEGCleanProducer = producers.icMetFromPatProducer.clone(
  #input = cms.InputTag("slimmedMETsMuEGClean"),
  ##input = cms.InputTag("patPFMetT1"),
  #branch = cms.string("pfMetType1MuEGCleanCollection"),
  #includeCustomID = cms.bool(False),
  #inputCustomID = cms.InputTag(""),
  #includeExternalMetsigMethod2 = cms.bool(True),
  #includeMetCorrections = cms.bool(True),
  #includeMetUncertainties = cms.bool(True)
  #)

#puppi met, need to update significance....

#from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
#makePuppiesFromMiniAOD( process );

#runMetCorAndUncFromMiniAOD(process,
#                           isData=bool(isData),
#                           metType="Puppi",
#                           pfCandColl=cms.InputTag("puppiForMET"),
#                           recoMetFromPFCs=True,
#                           reclusterJets=True,
#                           jetFlavor="AK4PFPuppi",
#                           postfix="Puppi"
#                           )

process.icPfMetPuppiProducer = producers.icMetFromPatProducer.clone(
  #input = cms.InputTag("patPFMetT1Puppi"),
  input = cms.InputTag("slimmedMETsPuppi"),
  branch = cms.string("pfMetPuppiCollection"),
  includeCustomID = cms.bool(False),
  inputCustomID = cms.InputTag(""),
  includeMetCorrections = cms.bool(False),
  includeMetUncertainties = cms.bool(False)
  )

# process.ictype1PfMetProducermetsigoutofbox = producers.icMetProducer.clone(
#                                                     input = cms.InputTag("slimmedMETs"),
#                                                     branch = cms.string("pfMetType1metsigoutofbox"),
#                                                     includeCustomID = cms.bool(False),
#                                                     inputCustomID = cms.InputTag(""),
#                                                     includeExternalMetsig = cms.bool(False)
#                                                     )

## Following lines are for default MET for Type1 corrections.
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

# If you only want to re-correct for JEC and get the proper uncertainties for the default MET
runMetCorAndUncFromMiniAOD(process,
                           isData=bool(isData))

process.load('RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff')
process.badGlobalMuonTaggerMAOD.muons       = muonLabel
process.badGlobalMuonTaggerMAOD.vtx         = vtxLabel
process.badGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)
process.cloneGlobalMuonTaggerMAOD.muons       = muonLabel
process.cloneGlobalMuonTaggerMAOD.vtx         = vtxLabel
process.cloneGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)

#from PhysicsTools.PatUtils.tools.muonRecoMitigation import muonRecoMitigation

#muonRecoMitigation(process = process,
#                   pfCandCollection = "packedPFCandidates", #input PF Candidate Collection
#                   runOnMiniAOD = True, #To determine if you are running on AOD or MiniAOD
#                   selection="", #You can use a custom selection for your bad muons. Leave empty if you would like to use the bad muon recipe definition.
#                   muonCollection="", #The muon collection name where your custom selection will be applied to. Leave empty if you would like to use the bad muon recipe definition.
#                   cleanCollName="cleanMuonsPFCandidates", #output pf candidate collection ame
#                   cleaningScheme="computeAllApplyClone", #Options are: "all", "computeAllApplyBad","computeAllApplyClone". Decides which (or both) bad muon collections to be used for MET cleaning coming from the bad muon recipe.
#                   postfix="" #Use if you would like to add a post fix to your muon / pf collections
#                   )

#runMetCorAndUncFromMiniAOD(process,
#                           isData=bool(isData),
#                           pfCandColl="cleanMuonsPFCandidates",
#                           recoMetFromPFCs=True,
#                           postfix="MuClean"
#                           )

#process.mucorMET = cms.Sequence(process.badGlobalMuonTaggerMAOD *
#                                process.cloneGlobalMuonTaggerMAOD *
#                                #process.badMuons * # If you are using cleaning mode "all", uncomment this line
#                                process.cleanMuonsPFCandidates *
#                                process.fullPatMetSequenceMuClean
#                                )


# Now you are creating the e/g corrected MET on top of the bad muon corrected MET (on re-miniaod)
from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG
corMETFromMuonAndEG(process,
                    pfCandCollection="", #not needed
                    electronCollection="slimmedElectronsBeforeGSFix",
                    photonCollection="slimmedPhotonsBeforeGSFix",
                    corElectronCollection="slimmedElectrons",
                    corPhotonCollection="slimmedPhotons",
                    allMETEGCorrected=True,
                    muCorrection=False,
                    eGCorrection=True,
                    runOnMiniAOD=True,
                    postfix="MuEGClean"
                    )
process.slimmedMETsMuEGClean = process.slimmedMETs.clone()
process.slimmedMETsMuEGClean.src = cms.InputTag("patPFMetT1MuEGClean")
process.slimmedMETsMuEGClean.rawVariation =  cms.InputTag("patPFMetRawMuEGClean")
process.slimmedMETsMuEGClean.t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
del process.slimmedMETsMuEGClean.caloMET
# If you are running in the scheduled mode:
process.egcorrMET = cms.Sequence(
  process.cleanedPhotonsMuEGClean+
  process.cleanedCorPhotonsMuEGClean+
  process.matchedPhotonsMuEGClean+ 
  process.matchedElectronsMuEGClean +
  process.corMETPhotonMuEGClean+
  process.corMETElectronMuEGClean+
  process.patPFMetT1MuEGClean+
  process.patPFMetRawMuEGClean+
  process.patPFMetT1SmearMuEGClean+
  process.patPFMetT1TxyMuEGClean+
  process.patPFMetTxyMuEGClean+
  process.patPFMetT1JetEnUpMuEGClean+
  process.patPFMetT1JetResUpMuEGClean+
  process.patPFMetT1SmearJetResUpMuEGClean+
  process.patPFMetT1ElectronEnUpMuEGClean+
  process.patPFMetT1PhotonEnUpMuEGClean+
  process.patPFMetT1MuonEnUpMuEGClean+
  process.patPFMetT1TauEnUpMuEGClean+
  process.patPFMetT1UnclusteredEnUpMuEGClean+
  process.patPFMetT1JetEnDownMuEGClean+
  process.patPFMetT1JetResDownMuEGClean+
  process.patPFMetT1SmearJetResDownMuEGClean+
  process.patPFMetT1ElectronEnDownMuEGClean+
  process.patPFMetT1PhotonEnDownMuEGClean+
  process.patPFMetT1MuonEnDownMuEGClean+
  process.patPFMetT1TauEnDownMuEGClean+
  process.patPFMetT1UnclusteredEnDownMuEGClean+
  process.slimmedMETsMuEGClean)

#process.icMetSequenceCorrected = cms.Sequence(
#  process.mucorMET  *
#  process.fullPatMetSequence * # If you are re-correctign the default MET
#  process.egcorrMET
#)

process.icMetSequence = cms.Sequence(
  process.METSignificance+
  process.icPfMet+
  process.icuncorrectedPfMetProducer+
  process.ictype1PfMetProducer+
  process.icPfMetPuppiProducer
  #+process.icMvaMetSequence
  #process.ictype1PfMetProducermetsigoutofbox
  #process.correctionTermsPfMetType1Type2+ #!!needs particle flow, need to find appropriate bit and change to packed version
  #process.correctionTermsPfMetType0PFCandidate + #!!currently causing errors
  #process.pfMetT0pcT1+
  #process.icPfMetT0pcT1Producer
)

#if isData:
  #process.icMetSequence += cms.Sequence(process.ictype1PfMetMuEGCleanProducer)

#!!MET UNCERTAINTIES
#!!MET FILTERS


################################################################                                                                                            
#!! Simulation only: GenParticles, GenJets, PileupInfo                                                                                                        
################################################################                                                                                            
process.icGenSequence = cms.Sequence()

#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.icPrunedGenParticles = cms.EDProducer("GenParticlePruner",
#  src = cms.InputTag("prunedGenParticles","","PAT"),
#  select = cms.vstring(
#    "drop  *",
#    "keep status == 3 || status == 22 || status == 23",  # all status 3
#    "keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15",  # all charged leptons
#    "keep abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16",  # all neutrinos
#    "keep++ abs(pdgId) == 15",  # keep full tau decay chain
#    "keep (4 <= abs(pdgId) <= 6)", # keep heavy flavour quarks !!note extended to top quarks as well
#    "keep (21 <= abs(pdgId) <= 25)", # keep bosons
#    "keep (400 <= abs(pdgId) < 600) || (4000 <= abs(pdgId) < 6000)", # keep b and c hadrons
#    "keep abs(pdgId) = 10411 || abs(pdgId) = 10421 || abs(pdgId) = 10413 || abs(pdgId) = 10423 || abs(pdgId) = 20413 || abs(pdgId) = 20423 || abs(pdgId) = 10431 || abs(pdgId) = 10433 || abs(pdgId) = 20433", # additional c hadrons for jet fragmentation studies
#    "keep abs(pdgId) = 10511 || abs(pdgId) = 10521 || abs(pdgId) = 10513 || abs(pdgId) = 10523 || abs(pdgId) = 20513 || abs(pdgId) = 20523 || abs(pdgId) = 10531 || abs(pdgId) = 10533 || abs(pdgId) = 20533 || abs(pdgId) = 10541 || abs(pdgId) = 10543 || abs(pdgId) = 20543" # additional b hadrons for jet fragmentation studies
#  )
#)

#process.prunedGenParticlesTaus = cms.EDProducer("ICGenParticlePruner53X",
#  src = cms.InputTag("prunedGenParticles","","PAT"),
#  select = cms.vstring(
#    "drop  *",
#    "keep++ abs(pdgId) == 15",  # keep full tau decay chain
#  )
#)

#Store ALL pruned gen particles
process.icGenParticleProducer = producers.icGenParticleProducer.clone(
  input   = cms.InputTag("prunedGenParticles","","PAT"),
  includeMothers = cms.bool(True),
  includeDaughters = cms.bool(True),
  includeStatusFlags = cms.bool(True)
)

process.icGenParticleProducerFromLHEParticles = producers.icGenParticleFromLHEParticlesProducer.clone(
##  input   = cms.InputTag("externalLHEProducer","","LHE")
)


#process.icGenParticleTauProducer = producers.icGenParticleProducer.clone(
#  input   = cms.InputTag("prunedGenParticlesTaus"),
#  branch = cms.string("genParticlesTaus"),
#  includeMothers = cms.bool(True),
#  includeDaughters = cms.bool(True),
#  includeStatusFlags = cms.bool(True)
#)

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.genParticlesForJets.ignoreParticleIDs = cms.vuint32(
  1000022, 2000012, 2000014,
  2000016, 1000039, 5000039,
  4000012, 9900012, 9900014,
  9900016, 39, 12, 14, 16
)
if release in ['80XMINIAOD', '92XMINIAOD']:
  # But of course this won't work because genParticlesForJets(InputGenJetsParticleSelector)
  # requires a vector<GenParticle> input. There's no alternative filter for the PackedGenParticle
  # type at the moment. Probably we could make our own generic cut-string selector, but
  # not in this package
  process.genParticlesForJets.src = cms.InputTag("packedGenParticles")

process.load("RecoJets.JetProducers.ak4GenJets_cfi")
process.ak4GenJetsNoNuBSM  =  process.ak4GenJets.clone()
process.ak4GenJetsNoNuBSM.src=cms.InputTag("packedGenParticles") #This still contains nus in 72, should be fixed in 74

process.selectedGenJets = cms.EDFilter("GenJetRefSelector",
  src = cms.InputTag("ak4GenJetsNoNuBSM"),
  cut = cms.string("pt > 10.0")
)

process.icGenJetProducer = producers.icGenJetProducer.clone(
  branch = cms.string("genJetsReclustered"),
  input   = cms.InputTag("selectedGenJets"),
  inputGenParticles = cms.InputTag("prunedGenParticles","","PAT"),
  requestGenParticles = cms.bool(False),
  isSlimmed  = cms.bool(True)#set true because based on gen jets made from packedcandidates
)

process.icGenJetProducerFromSlimmed = producers.icGenJetProducer.clone(
  branch = cms.string("genJets"),
  input = cms.InputTag("slimmedGenJets"),
  inputGenParticles=cms.InputTag("genParticles"),
  requestGenParticles = cms.bool(False),
  isSlimmed = cms.bool(True)
  ) 

process.icPileupInfoProducer = producers.icPileupInfoProducer.clone(
  input   = cms.InputTag("slimmedAddPileupInfo")
  )


if not isData:
  process.icGenSequence += (
    #process.icPrunedGenParticles+
    #process.prunedGenParticlesTaus+
    #process.icGenParticleTauProducer+
    process.icGenParticleProducer
  )
  process.icGenSequence += (
      process.ak4GenJetsNoNuBSM+
      process.selectedGenJets+
      process.icGenJetProducer+
      process.icGenJetProducerFromSlimmed+
      process.icGenMetProducer+
      process.icPileupInfoProducer
    )

if doHT:
    process.icGenSequence += (
      process.icGenParticleProducerFromLHEParticles
)
##################################################################                                                                                          
#!! Trigger                                                                                                                                                 
#!! Need to get L1 extra
################################################################## 

from PhysicsTools.PatAlgos.tools.trigTools import *
process.icTriggerSequence = cms.Sequence()
process.icTriggerObjectSequence = cms.Sequence()

process.patTriggerPath = cms.Path()
if release in ['80XMINIAOD', '92XMINIAOD']:
  switchOnTrigger(process, path = 'patTriggerPath',  outputModule = '')

process.icTriggerPathProducer = producers.icTriggerPathProducer.clone(
  branch = cms.string("triggerPaths"),
  inputIsStandAlone = cms.bool(True),
  input = cms.InputTag("TriggerResults", "", "HLT"),
  inputPrescales=cms.InputTag("patTrigger")
  )

if isReHLT:
  process.icTriggerPathProducer.input = cms.InputTag("TriggerResults", "", "HLT2")


process.icTriggerSequence += cms.Sequence(
  #process.patTrigger+
  #process.patTriggerEvent+
  process.icTriggerPathProducer
  )

## It produces failures in the production of DATA_ntuples


if not isReHLT:
  #process.icPFMET170NoiseCleanedObjectProducer = producers.icTriggerObjectProducer.clone(
    #branch = cms.string("triggerObjectsPFMET170NoiseCleaned"),
    #input   = cms.InputTag("patTriggerEvent"),
    #hltPath = cms.string("HLT_PFMET170_NoiseCleaned"),
    #inputIsStandAlone = cms.bool(False),
    #storeOnlyIfFired = cms.bool(True)
  #)

  process.icDiPFJet40DEta3p5MJJ600PFMETNoMu140ObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsDiPFJet40DEta3p5MJJ600PFMETNoMu140"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140"),
    inputIsStandAlone = cms.bool(False),
    storeOnlyIfFired = cms.bool(True)
    )

  process.icDiPFJet40DEta3p5MJJ600PFMETNoMu80ObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsDiPFJet40DEta3p5MJJ600PFMETNoMu80"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80"),
    inputIsStandAlone = cms.bool(False),
    storeOnlyIfFired = cms.bool(True)
    )

  process.icPFMETNoMu90PFMHTNoMu90IDTightObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsPFMETNoMu90PFMHTNoMu90IDTight"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight"),
    inputIsStandAlone = cms.bool(False),
    storeOnlyIfFired = cms.bool(True)
    )

  process.icPFMETNoMu100PFMHTNoMu100IDTightObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsPFMETNoMu100PFMHTNoMu100IDTight"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight"),
    inputIsStandAlone = cms.bool(False),
    storeOnlyIfFired = cms.bool(True)
    )

  process.icPFMETNoMu110PFMHTNoMu110IDTightObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsPFMETNoMu110PFMHTNoMu110IDTight"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight"),
    inputIsStandAlone = cms.bool(False),
    storeOnlyIfFired = cms.bool(True)
    )

  process.icPFMETNoMu120PFMHTNoMu120IDTightObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsPFMETNoMu120PFMHTNoMu120IDTight"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"),
    inputIsStandAlone = cms.bool(False),
    storeOnlyIfFired = cms.bool(True)
    )

  process.icIsoMu27ObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsIsoMu27"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_IsoMu27"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

  process.icPFMET120PFMHT120IDTightObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsPFMET120PFMHT120IDTight"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_PFMET120_PFMHT120_IDTight"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

  process.icDiJet11035Mjj650PFMET110ObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsDiJet11035Mjj650PFMET110"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_DiJet110_35_Mjj650_PFMET110"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

if isReHLT:
  process.icDiPFJet40DEta3p5MJJ600PFMETNoMu140ObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsDiPFJet40DEta3p5MJJ600PFMETNoMu140"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140"),
    inputIsStandAlone = cms.bool(False),
    storeOnlyIfFired = cms.bool(True),
    inputTriggerResults = cms.InputTag("TriggerResults", "", "HLT2")
    )

  process.icDiPFJet40DEta3p5MJJ600PFMETNoMu80ObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsDiPFJet40DEta3p5MJJ600PFMETNoMu80"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80"),
    inputIsStandAlone = cms.bool(False),
    storeOnlyIfFired = cms.bool(True),
    inputTriggerResults = cms.InputTag("TriggerResults", "", "HLT2")
    )

  process.icPFMETNoMu90PFMHTNoMu90IDTightObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsPFMETNoMu90PFMHTNoMu90IDTight"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight"),
    inputIsStandAlone = cms.bool(False),
    storeOnlyIfFired = cms.bool(True),
    inputTriggerResults = cms.InputTag("TriggerResults", "", "HLT2")
    )

  process.icPFMETNoMu100PFMHTNoMu100IDTightObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsPFMETNoMu100PFMHTNoMu100IDTight"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight"),
    inputIsStandAlone = cms.bool(False),
    storeOnlyIfFired = cms.bool(True),
    inputTriggerResults = cms.InputTag("TriggerResults", "", "HLT2")
    )

  process.icPFMETNoMu110PFMHTNoMu110IDTightObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsPFMETNoMu110PFMHTNoMu110IDTight"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight"),
    inputIsStandAlone = cms.bool(False),
    storeOnlyIfFired = cms.bool(True),
    inputTriggerResults = cms.InputTag("TriggerResults", "", "HLT2")
    )

  process.icPFMETNoMu120PFMHTNoMu120IDTightObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsPFMETNoMu120PFMHTNoMu120IDTight"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"),
    inputIsStandAlone = cms.bool(False),
    storeOnlyIfFired = cms.bool(True),
    inputTriggerResults = cms.InputTag("TriggerResults", "", "HLT2")
    )

  process.icIsoMu27ObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsIsoMu27"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_IsoMu27"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False),
    inputTriggerResults = cms.InputTag("TriggerResults", "", "HLT2")
    )

  process.icPFMET120PFMHT120IDTightObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsPFMET120PFMHT120IDTight"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_PFMET120_PFMHT120_IDTight"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False),
    inputTriggerResults = cms.InputTag("TriggerResults", "", "HLT2")
    )

  process.icDiJet11035Mjj650PFMET110ObjectProducer = producers.icTriggerObjectProducer.clone(
    branch = cms.string("triggerObjectsDiJet11035Mjj650PFMET110"),
    input   = cms.InputTag("patTriggerEvent"),
    hltPath = cms.string("HLT_DiJet110_35_Mjj650_PFMET110"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False),
    inputTriggerResults = cms.InputTag("TriggerResults", "", "HLT2")
    )


process.icTriggerObjectSequence += cms.Sequence(
  #process.icPFMET170NoiseCleanedObjectProducer+
  process.icDiPFJet40DEta3p5MJJ600PFMETNoMu80ObjectProducer+
  process.icDiPFJet40DEta3p5MJJ600PFMETNoMu140ObjectProducer+
  process.icPFMETNoMu90PFMHTNoMu90IDTightObjectProducer+
  process.icPFMETNoMu100PFMHTNoMu100IDTightObjectProducer+
  process.icPFMETNoMu110PFMHTNoMu110IDTightObjectProducer+
  process.icPFMETNoMu120PFMHTNoMu120IDTightObjectProducer
  )

if release in ['80XMINIAOD', '92XMINIAOD']:
  for name in process.icTriggerObjectSequence.moduleNames():
    mod = getattr(process, name)
    mod.inputIsStandAlone = cms.bool(True)
    mod.input = cms.InputTag("selectedPatTrigger")
  if isReHLT:
    for name in process.icTriggerObjectSequence.moduleNames():
      mod = getattr(process, name)
      mod.inputTriggerResults = cms.InputTag("TriggerResults", "","HLT2")

## Need to unpack filterLabels on slimmedPatTrigger then make selectedPatTrigger
process.patTriggerUnpacker = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
   patTriggerObjectsStandAlone = cms.InputTag("slimmedPatTrigger"),
   triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
   unpackFilterLabels = cms.bool(True)
)
process.selectedPatTrigger = cms.EDFilter(
 'PATTriggerObjectStandAloneSelector',
  cut = cms.string('!filterLabels.empty()'),
  src = cms.InputTag('patTriggerUnpacker')
)
process.icTriggerSequence += cms.Sequence(
    process.patTriggerUnpacker +
    process.selectedPatTrigger
    )

#L1 Extra information

if isData or isReHLT:
  process.icL1EtSumProducer = cms.EDProducer('ICL1TObjectProducer<l1t::EtSum>',
    branch = cms.string("l1tEtSum"),
    input = cms.InputTag("caloStage2Digis","EtSum","RECO"),
    doBXloop = cms.bool(v_doBXloop)
  )
  if isReHLT:
    process.icL1EtSumProducer.input = cms.InputTag("caloStage2Digis","EtSum","HLT2")

  process.icL1MuonsProducer = cms.EDProducer('ICL1TObjectProducer<l1t::Muon>',
    branch = cms.string("l1tMuon"),
    input = cms.InputTag("gmtStage2Digis","Muon","RECO"),
    doBXloop = cms.bool(v_doBXloop)
  )
  if isReHLT:
    process.icL1MuonsProducer.input = cms.InputTag("gmtStage2Digis","Muon","HLT2")

  process.icL1EGammaProducer = cms.EDProducer('ICL1TObjectProducer<l1t::EGamma>',
    branch = cms.string("l1tEGamma"),
    input = cms.InputTag("caloStage2Digis","EGamma","RECO"),
    doBXloop = cms.bool(v_doBXloop)
  )
  if isReHLT:
    process.icL1EGammaProducer.input = cms.InputTag("caloStage2Digis","EGamma","HLT2")

  process.icL1JetProducer = cms.EDProducer('ICL1TObjectProducer<l1t::Jet>',
    branch = cms.string("l1tJet"),
    input = cms.InputTag("caloStage2Digis","Jet","RECO"),
    doBXloop = cms.bool(v_doBXloop)
  )
  if isReHLT:
    process.icL1JetProducer.input = cms.InputTag("caloStage2Digis","Jet","HLT2")

  process.icL1Sequence = cms.Sequence(
      process.icL1EtSumProducer+
      process.icL1MuonsProducer+
      process.icL1EGammaProducer+
      process.icL1JetProducer
      )



################################################################
# EventInfo
################################################################
#Load the MET filters here
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')

process.BadPFMuonFilter.muons        = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.taggingMode  = cms.bool(True)
process.BadChargedCandidateFilter.muons        = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.taggingMode  = cms.bool(True)


# So in this case the 100 events you are checking simply do not contain any bad muons.
# Because these modules select bad events, the event fails the filters if it does not have any bad muons.
# This is why you arre using ~ to invert the behaviour of the filters.
# If you look at  https://github.com/gpetruc/cmssw/blob/badMuonFilters_80X_v2/RecoMET/METFilters/plugins/badGlobalMuonTaggerMAOD.cc#L118-L125 you see that the opposite is true when you switch on tagging mode: if tagging mode is set to true then events always pass the filter, and so if you use ~ all of the events will fail the filter.
# So you should redefine noBadGlobalMuons as (process.badGlobalMuonTaggerMAOD + process.cloneGlobalMuonTaggerMAOD).
# See also https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2786/2/1.html.

process.icEventInfoProducer = producers.icEventInfoProducer.clone(
  includeJetRho       = cms.bool(True),
  includeHT           = cms.bool(False),
  lheProducer         = cms.InputTag("externalLHEProducer"),
  inputJetRho         = cms.InputTag("fixedGridRhoFastjetAll"),
  includeLeptonRho    = cms.bool(False),
  inputLeptonRho      = cms.InputTag("fixedGridRhoFastjetAll"),
  includeVertexCount  = cms.bool(True),
  inputVertices       = vtxLabel,
  includeCSCFilter    = cms.bool(False),
  inputCSCFilter      = cms.InputTag("BeamHaloSummary"),
  includeFiltersFromTrig = cms.bool(True),
  inputfiltersfromtrig = cms.InputTag("TriggerResults","","RECO"),
  filters             = cms.PSet(
    badChargedHadronFilter = cms.InputTag("BadChargedCandidateFilter"),
    badMuonFilter          = cms.InputTag("BadPFMuonFilter"),
    Flag_badMuons          = cms.InputTag("badGlobalMuonTaggerMAOD"),
    Flag_duplicateMuons    = cms.InputTag("cloneGlobalMuonTaggerMAOD")
    ),
  filtersfromtrig     = cms.vstring("*")
)

if not isData:
  process.icEventInfoProducer.inputfiltersfromtrig = cms.InputTag("TriggerResults","","PAT")

if isData:
  process.icEventInfoProducer.filters=cms.PSet(
    badChargedHadronFilter = cms.InputTag("BadChargedCandidateFilter"),
    badMuonFilter          = cms.InputTag("BadPFMuonFilter")
  )



#if doHT:
 # process.icEventInfoProducer.includeHT = cms.bool(True)

 
process.icEventInfoSequence = cms.Sequence(
  process.BadPFMuonFilter+
  process.BadChargedCandidateFilter+ 
  process.badGlobalMuonTaggerMAOD+
  process.cloneGlobalMuonTaggerMAOD+
  process.icEventInfoProducer
  )

if isData:
  process.icEventInfoSequence.remove(process.badGlobalMuonTaggerMAOD)
  process.icEventInfoSequence.remove(process.cloneGlobalMuonTaggerMAOD)

################################################################
# Event
################################################################
process.icEventProducer = producers.icEventProducer.clone()

process.icProdSeq = cms.Sequence()

#if isData:
#  process.icProdSeq += cms.Sequence(process.icMetSequenceCorrected)

process.icProdSeq += cms.Sequence(
  process.ic80XSequence+
  process.icMiniAODSequence+
  process.icSelectionSequence+
  process.pfParticleSelectionSequence+
  process.icVertexSequence+
  process.icElectronSequence+
  process.icMuonSequence+
  process.icTauSequence+
  process.icPhotonSequence
)

if isData or isReHLT:
  process.icProdSeq += process.icL1TauProducer

process.icProdSeq += cms.Sequence(
  process.icPFJetSequence+
  process.icGenSequence+
  process.icMetSequence+
  process.icTriggerSequence+
  process.icTriggerObjectSequence+
  process.icEventInfoSequence
)

if isData or isReHLT:
  process.icProdSeq += process.icL1Sequence

process.icProdSeq += process.icEventProducer

process.p = cms.Path(process.icProdSeq)

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string("edmdump.root")
#                               )

#process.outpath = cms.EndPath(process.out)


#process.schedule = cms.Schedule(process.patTriggerPath, process.p)
process.schedule = cms.Schedule(process.p)

#make an edm output ntuple with everything in it
#process.schedule = cms.Schedule(process.p,process.outpath)



#Uncomment below for a dump of the config
#print process.dumpPython()

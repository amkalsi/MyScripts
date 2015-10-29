import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMuon.MuonIdentification.me0MuonReco_cff')
process.load('Validation.RecoMuon.associators_cff')
process.load('Validation.RecoMuon.selectors_cff')
process.load('Validation.RecoMuon.MuonTrackValidator_cfi')
process.load('Validation.RecoMuon.RecoMuonValidator_cfi')
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V6::All', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        'root://xrootd.unl.edu//store/user/amkalsi/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola_HGCALGS_PU140_MinBias_ME0_RECO_TimingResolution1ps/b0b349462b87ba52be1798b06a8d86fb/out_reco_1005_1_YLW.root'
        )
                            )
from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
from SimTracker.TrackAssociation.TrackAssociatorByPosition_cff import *
#need a propagator in case analysis is made a radius!=0
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import *

from SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi import *
process.TrackAssociatorByChi2ESProducer = TrackAssociatorByChi2ESProducer.clone(chi2cut = 100.0,ComponentName = 'TrackAssociatorByChi2')
process.me0muon = cms.EDProducer("ME0MuonTrackCollProducer",
#                         me0MuonTag = cms.InputTag("me0SegmentMatching"),
#                         selectionTags = cms.vstring('All'),
                                 muonQ = cms.untracked.bool(False),
                                 # true for tight muons, false for loose muons 
                                 )
#--------------------

process.me0muonColl_seq = cms.Sequence(
                             process.me0muon
                             )
process.demo = cms.EDAnalyzer("ME0TimingAnalysis",
                      
                              UseAssociators = cms.bool(True),
                              associators = cms.vstring('TrackAssociatorByChi2'),
                              label = cms.VInputTag('me0muon'),  
     # input TrackingParticle collections

                              )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("OutTree.root")
                                   )


process.recoMuonValidation = cms.Sequence(#probeTracks_seq*
    #(selectedVertices * selectedFirstPrimaryVertex) * 
    #bestMuonTuneP_seq*
    #muonColl_seq*trackColl_seq*extractedMuonTracks_seq*bestMuon_seq*trackerMuon_seq*
    process.me0muonColl_seq
    #((process.muonValidation_seq))
    )

process.p = cms.Path(process.recoMuonValidation*process.demo)

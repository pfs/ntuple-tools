import FWCore.ParameterSet.Config as cms

import sys
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')
options.register('inputFiles', None,               VarParsing.multiplicity.singleton, VarParsing.varType.string, "input directory")
options.register('outputFile', 'ROIanalysis.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, "output filey")
options.parseArguments()

process = cms.Process("HGCROI")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoLocalCalo.HGCalRecProducers.HGCalLocalRecoSequence_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
from FastSimulation.Event.ParticleFilter_cfi import *
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import *

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                            )

try:
    import os
    inDir=options.inputFiles
    if not '.root' in inDir:
        allFiles=[os.path.join(inDir,f) for f in os.listdir('/eos/cms/%s'%inDir)]
        process.source.fileNames.extend(allFiles)
        print len(process.source.fileNames),'files found in',inDir 
    else:
        print 'Processing only',inDir
        inDir=inDir.replace('/eos/cms','')
        process.source.fileNames.extend([inDir])

except:
    print 'could not parse input directory, using default'
    process.source.fileNames.extend(['/store/cmst3/user/psilva/HGCal/H125gg_EE/CMSSW_9_3_2/DIGI_PU00_0p0/RECO/Events_24.root'])


process.ana = cms.EDAnalyzer("HGCROIAnalyzer",
                             dEdXWeights = dEdX_weights,
                             thicknessCorrection = cms.vdouble(1.132,1.092,1.084),
                             byClosest = cms.bool(False)
                             )
                             
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outputFile)
                                   )
process.p = cms.Path(process.ana)

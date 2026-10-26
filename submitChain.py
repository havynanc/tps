import os


data = {
#    'MuonGUN':'/DoubleMuon_FlatPt-1To100-gun/Phase2Fall22DRMiniAOD-PU200_125X_mcRun4_realistic_v2-v1/GEN-SIM-DIGI-RAW-MINIAOD',
#    'TTBar':'/TT_TuneCP5_14TeV-powheg-pythia8/Phase2Spring24DIGIRECOMiniAOD-PU200_AllTP_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD',
     'DY200_23':'/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2Spring23DIGIRECOMiniAOD-PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/GEN-SIM-DIGI-RAW-MINIAOD',
#     'DY200_24':'/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2Spring24DIGIRECOMiniAOD-PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD',
#     'MinBias200_23_partialFalse':'/MinBias_TuneCP5_14TeV-pythia8/Phase2Spring23DIGIRECOMiniAOD-PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/GEN-SIM-DIGI-RAW-MINIAOD',
#    'Displaced':'/DisplacedMuons_Pt-30To100_Dxy-0To3000-gun/Phase2Fall22DRMiniAOD-PU200_125X_mcRun4_realistic_v2-v1/GEN-SIM-DIGI-RAW-MINIAOD',   
     'Jpsi_23':'/JPsiToMuMu_Pt-0to100_TuneCP5-pythia8-gun/Phase2Spring23DIGIRECOMiniAOD-PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/GEN-SIM-DIGI-RAW-MINIAOD',
#     'Jpsi_24':'/JPsiToMuMu_Pt-0To100-pythia8-gun/Phase2Spring24DIGIRECOMiniAOD-PU200_Trk1GeV_140X_mcRun4_realistic_v4-v2/GEN-SIM-DIGI-RAW-MINIAOD'
}

#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 1
for tag,dataset in  data.items():
    FILE="""
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'skim_{tag}'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runOnlyGMT.py'
config.Data.inputDataset = '{dataset}'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/hancelin/new_tflayer_stubs'
config.Data.publication = True
config.Data.ignoreLocality = True
config.Data.outputDatasetTag = 'PHASEII_{tag}'
config.Data.partialDataset = False
config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.whitelist = ['T2_US_*']
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 3000
""".format(tag=tag,dataset=dataset)
    f=open("crab_{tag}.py".format(tag=tag),"w")
    print(FILE)
    f.write(FILE)
    f.close()
    os.system("crab submit -c crab_{PT}.py".format(PT=tag))



from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.workArea='Mar6'

config.section_('JobType')
config.JobType.psetName = '/afs/cern.ch/work/a/agilbert/CMSSW_TEST/CMSSW_5_3_11/src/UserCode/ICHiggsTauTau/test/higgstautau_new_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['EventTree.root']
config.JobType.pyCfgParams = ['isData=0', 'release=53X', 'globalTag=START53_V22::All']

config.section_('Data')
config.Data.unitsPerJob = 20000
config.Data.totalUnits = 2000000
#config.Data.unitsPerJob = 1
config.Data.splitting = 'EventAwareLumiBased'
config.Data.publication = False
config.Data.outputDatasetTag = 'Mar6Pilot2'

config.section_('User')
config.User.voGroup = 'dcms'

config.section_('Site')
config.Site.storageSite = 'T2_DE_DESY'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print hte.headers

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    tasks=list()

    tasks.append(('QCD_Pt_20_MuEnrichedPt_15','/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM'))

    for task in tasks:
        print task[0]
        config.General.requestName = task[0]
        config.Data.inputDataset = task[1]
        submit(config)




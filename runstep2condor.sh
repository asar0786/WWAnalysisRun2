#!/bin/bash
echo "Starting job on " `date`
echo "Running on: `uname -a`"
echo "System software: `cat /etc/redhat-release`"
source /cvmfs/cms.cern.ch/cmsset_default.sh
### copy the input root files if they are needed only if you require local reading
xrdcp -s root://cmseos.fnal.gov//store/user/asahmed/SecondStep/WWTree_CommonNtuple_For1and2Lepton_2019_05_27_05h39/CMSSW_8_0_26_patch1.tgz  .
tar -xf CMSSW_8_0_26_patch1.tgz
rm CMSSW_8_0_26_patch1.tgz
cd CMSSW_8_0_26_patch1/src/WWAnalysis/WWAnalysisRun2
echo "====> List files : " 
ls -alh
echo "====> Remove any file with name similar to WWTree*.root... " 
rm WWTree*.root
scramv1 b ProjectRename
eval `scram runtime -sh`
echo "====> List files : " 
ls -alh
python python/produceWWNtuples.py -i /store/user/lnujj/WpWm_aQGC_Ntuples_Ram/FirstStepOutput/BaconNtuples/ $*
echo "====> List files : " 
ls -alh
echo "====> List root files : " 
ls WWTree*.root
echo "====> copying WWTree*.root file to stores area..." 
xrdcp -f WWTree*.root root://cmseos.fnal.gov//store/user/asahmed/SecondStep/WWTree_CommonNtuple_For1and2Lepton_2019_05_27_05h39
rm WWTree*.root
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_8_0_26_patch1

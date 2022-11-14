#!/usr/bin/env bash

# make sure that you are using O2 from CVMFS
# via for e.g. 
# > /cvmfs/alice.cern.ch/bin/alienv enter VO_ALICE@O2Physics::nightly-20221107-1

#baseDir="/Users/andry/cernbox"
baseDir="/afs/cern.ch/user/a/arakotoz/mycernbox"
myDir="${baseDir}/alice/mft/pilotbeam/505713/test-reco"
cd $myDir
pwd

## Ensure necessary files are present in the node work directory
#ccdbBaseDir="${myDir}/ccdb"
#geomFilenameAligned="snapshot.root"
#ccdbGeomAlignedPath="GLO/Config/GeometryAligned"

#if [ ! -e "${ccdbBaseDir}/${ccdbGeomAlignedPath}/${geomFilenameAligned}" ]; then
#    echo "Cannot find ${ccdbBaseDir}/${ccdbGeomAlignedPath}/${geomFilenameAligned}"
#    exit 1
#fi

inputfile="small-ctf-alien-file-list.txt"

shmSize=16000000000

severity="info"
logConfig="--severity ${severity} --timeframes-rate-limit 3 --timeframes-rate-limit-ipcid 0 " #--infologger-mode \"stdout\""

#readCmd="o2-ctf-reader-workflow --copy-cmd \"alien_cp ?src file://?dst\" --remote-regex \"^alien:///alice/data/.+\"  --ctf-input ${inputfile} --delay 1 --loop 0 --onlyDet MFT --shm-segment-size ${shmSize} ${logConfig} --allow-missing-detectors --condition-remap file://${ccdbBaseDir}=${ccdbGeomAlignedPath} "
readCmd="o2-ctf-reader-workflow --copy-cmd \"alien_cp ?src file://?dst\" --remote-regex \"^alien:///alice/data/.+\"  --ctf-input ${inputfile} --delay 1 --loop 0 --onlyDet MFT --shm-segment-size ${shmSize} ${logConfig} --allow-missing-detectors -b "

recoOptions="MFTTracking.FullClusterScan=true;MFTTracking.LTFclsRCut=0.2;"
#recoOptions="MFTTracking.FullClusterScan=true;MFTTracking.LTFclsRCut=0.2;MFTTracking.trackmodel=2;MFTAlpideParam.roFrameLengthInBC=198;"
#recoCmd="o2-mft-reco-workflow --shm-segment-size ${shmSize} ${logConfig} --nThreads 2 --clusters-from-upstream --mft-track-writer --mft-cluster-writer --disable-mc --pipeline mft-tracker:1 --run-assessment --configKeyValues \""${recoOptions}"\" --condition-remap file://${ccdbBaseDir}=${ccdbGeomAlignedPath} "
recoCmd="o2-mft-reco-workflow --shm-segment-size ${shmSize} ${logConfig} --nThreads 2 --clusters-from-upstream --mft-cluster-writer --disable-mc --pipeline mft-tracker:1 --run-assessment --configKeyValues \""${recoOptions}"\" -b "

# Concatenate workflow
runCmd=" $readCmd "
runCmd+=" | $recoCmd"
runCmd+=" | o2-dpl-run ${logConfig} --shm-segment-size ${shmSize} -b --run > ctf2cltrack.log"

# List input files and command line
echo "======================="
echo "Input file : ${inputfile}"
echo "number of files :"
cat  ${inputfile} | wc -l
echo "======================="
echo "Running reconstruction command: "
echo "${runCmd} "
startTime=$(date +"%Y %m %d %H:%M:%S")
echo "Start "${startTime}
echo "======================="

## Runs reconstruction comand stored in $runCmd
eval "${runCmd}"

endTime=$(date +"%Y %m %d %H:%M:%S")
echo "End "${endTime}
echo "======================="
#cd ${myDir}

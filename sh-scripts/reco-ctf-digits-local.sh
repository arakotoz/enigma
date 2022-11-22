#!/usr/bin/env bash

# make sure that you are using O2 from CVMFS
# via for e.g. 
# > /cvmfs/alice.cern.ch/bin/alienv enter VO_ALICE@O2Physics::nightly-20221017-1
## to execute, do:
# > . ~/cernbox/alice/enigma/sh-scripts/reco-ctf-digits-local.sh
## or:
# > source ~/cernbox/alice/enigma/sh-scripts/reco-ctf-digits-local.sh
#

MODE=$(echo $HOME | grep afs) # if empty, then using laptop, else running on lxplus
baseDir="/Users/andry/cernbox"
inputfile="small-ctf-local-file-list.txt"

echo "======================="

if [ ! -z "${MODE}" ]
then
    echo "MODE = lxplus"
    baseDir="/afs/cern.ch/user/a/arakotoz/mycernbox"
    inputfile="small-ctf-alien-file-list.txt"
else
    echo "MODE = laptop"
    echo "LD_LIBRARY_PATH="$LD_LIBRARY_PATH
fi

myDir="${baseDir}/alice/mft/LHC22h/520509/test-reco"
cd $myDir
echo "======================="
pwd

## Ensure necessary files are present in the node work directory
ccdbBaseDir="${myDir}/ccdb"
geomFilenameAligned="snapshot.root"
ccdbGeomAlignedPath="GLO/Config/GeometryAligned"

if [ ! -e "${ccdbBaseDir}/${ccdbGeomAlignedPath}/${geomFilenameAligned}" ]; then
    echo "Cannot find ${ccdbBaseDir}/${ccdbGeomAlignedPath}/${geomFilenameAligned}"
    exit 1
fi

magfieldFilename="snapshot.root"
ccdbMagfieldPath="GLO/Config/GRPMagField"

if [ ! -e "${ccdbBaseDir}/${ccdbMagfieldPath}/${magfieldFilename}" ]; then
    echo "Cannot find ${ccdbBaseDir}/${ccdbMagfieldPath}/${magfieldFilename}"
    exit 1
fi

noiseFilename="snapshot.root"
ccdbNoisePath="MFT/Calib/NoiseMap"

if [ ! -e "${ccdbBaseDir}/${ccdbNoisePath}/${noiseFilename}" ]; then
    echo "Cannot find ${ccdbBaseDir}/${ccdbNoisePath}/${noiseFilename}"
    exit 1
fi


shmSize=16000000000

severity="info"
logConfig="--severity ${severity} --timeframes-rate-limit 3 --timeframes-rate-limit-ipcid 0 " #--infologger-mode \"stdout\""

readCmd="o2-ctf-reader-workflow --copy-cmd no-copy --ctf-input ${inputfile} --delay 1 --loop 0 --onlyDet MFT --mft-digits --shm-segment-size ${shmSize} ${logConfig} --allow-missing-detectors --condition-remap file://${ccdbBaseDir}=${ccdbGeomAlignedPath},${ccdbMagfieldPath},${ccdbNoisePath} -b "

## laptop
recoOptions="MFTTracking.FullClusterScan=true;MFTTracking.LTFclsRCut=0.2;MFTTracking.trackmodel=2;MFTClustererParam.maxBCDiffToMaskBias=10;" # MFTClustererParam.maxBCDiffToSquashBias=-10;MFTAlpideParam.roFrameLengthInBC=396;"
threadOptions="--nThreads 2"
if [ ! -z "${MODE}" ]
then
    ## lxplus
    recoOptions="MFTTracking.forceZeroField=true;MFTTracking.FullClusterScan=true;MFTTracking.LTFclsRCut=0.2;MFTTracking.trackmodel=2;MFTClustererParam.maxBCDiffToMaskBias=10;"
    threadOptions=""
fi

recoCmd="o2-mft-reco-workflow --shm-segment-size ${shmSize} ${logConfig} ${threadOptions} --digits-from-upstream --mft-cluster-writer --disable-mc --pipeline mft-tracker:1 --run-assessment --configKeyValues \""${recoOptions}"\" --condition-remap file://${ccdbBaseDir}=${ccdbGeomAlignedPath},${ccdbMagfieldPath},${ccdbNoisePath} -b "

## Concatenate workflow
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
echo "======================="
startTime=$(date +"%Y %m %d %H:%M:%S")
echo "Start "${startTime}
echo "======================="

## Runs reconstruction comand stored in $runCmd
eval "${runCmd}" # "echo | " is a hack (to provide input stream to O2 workflows?)

endTime=$(date +"%Y %m %d %H:%M:%S")
echo "End "${endTime}
echo "======================="

if [ -z "${MODE}" ]
then 
    ## laptop
    cd ${myDir}
fi

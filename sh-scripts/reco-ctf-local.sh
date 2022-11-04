#!/usr/bin/env bash

# make sure that you are using O2 

myDir="/Users/andry/cernbox/alice/mft/pilotbeam/505713/test-reco"
cd $myDir
pwd
LD_LIBRARY_PATH=/usr/local/lib:/usr/lib
XROOTD_DIR=$(which xrootd-config | sed 's/bin\/xrootd-config//')
LD_LIBRARY_PATH=${XROOTD_DIR}/lib:/Users/andry/alice/sw/osx_x86-64/JAliEn-ROOT/latest/lib:${LD_LIBRARY_PATH}
echo $LD_LIBRARY_PATH
echo "======================="

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

inputfile="small-ctf-local-file-list.txt"

severity="info"
logConfig="--severity ${severity} --fairmq-rate-logging 0 " #--infologger-mode \"stdout\""
shmSize=16000000000

readCmd="o2-ctf-reader-workflow --copy-cmd no-copy --ctf-input ${inputfile} --delay 8 --onlyDet MFT --shm-segment-size ${shmSize} ${logConfig} --allow-missing-detectors --condition-remap file://${ccdbBaseDir}=${ccdbGeomAlignedPath},${ccdbMagfieldPath} -b "

recoOptions="MFTTracking.forceZeroField=true;MFTTracking.FullClusterScan=true;MFTTracking.LTFclsRCut=0.2;MFTTracking.trackmodel=2;MFTAlpideParam.roFrameLengthInBC=198;"
recoCmd="o2-mft-reco-workflow --shm-segment-size ${shmSize} ${logConfig} --clusters-from-upstream --mft-track-writer --mft-cluster-writer --disable-mc --pipeline mft-tracker:1 --run-assessment --configKeyValues \""${recoOptions}"\" --condition-remap file://${ccdbBaseDir}=${ccdbGeomAlignedPath},${ccdbMagfieldPath} "

# Concatenate workflow
runCmd=" $readCmd "
runCmd+=" | $recoCmd"
runCmd+=" | o2-dpl-run  ${logConfig} --shm-segment-size ${shmSize} -b --run > ctf2cltrack.log"
#runCmd+=" -b --run "

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
echo | eval "${runCmd}" # "echo | " is a hack (to provide input stream to O2 workflows?)

endTime=$(date +"%Y %m %d %H:%M:%S")
echo "End "${endTime}
echo "======================="
exit "$?"

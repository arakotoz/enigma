#!/usr/bin/env bash

currentDir=`pwd`

echo "* *****************************************************"
echo "* PATH: ${PATH}"
echo "* LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}"
echo "* Workdir: ${currentDir}"
echo "* *****************************************************"

## Ensure necessary files are present in the node work directory

geomFilenameAligned="snapshot.root"
if [ ! -e "${geomFilenameAligned}" ]; then
    echo "Cannot find $geomFilenameAligned"
    exit 1
fi

ccdbGeomAlignedPath="GLO/Config/GeometryAligned"
localCcdbBaseDir="${currentDir}/ccdb"
mkdir -p "${localCcdbBaseDir}/${ccdbGeomAlignedPath}"
cp -pfv "${geomFilenameAligned} ${localCcdbBaseDir}/${ccdbGeomAlignedPath}/."

if [ ! -e "${localCcdbBaseDir}/${ccdbGeomAlignedPath}/${geomFilenameAligned}" ]; then
    echo "Cannot find ${localCcdbBaseDir}/${ccdbGeomAlignedPath}/${geomFilenameAligned}"
    exit 1
fi

function makeFileList() { # Extract input file names from the XML file
    if [ "$3" -eq 1 ]; then
        grep -oE "turl=\"[a-zA-Z0-9/_.:]+\"" "$1" | cut -d '"' -f 2 >"$2"
    else
        grep -oE "turl=\"[a-zA-Z0-9/_.:]+\"" "$1" | cut -d '"' -f 2 | xargs -L 1 basename >"$2"
    fi
}

nodownload=1

inFilename="${1-wn.xml}"
tmpInFile="wn.txt"

makeFileList "${inFilename}" "${tmpInFile}" "${nodownload}"

shmSize=16000000000

severity="info"
logConfig="--severity ${severity} --resources-monitoring 50 --resources-monitoring-dump-interval 50 --early-forward-policy noraw --fairmq-rate-logging 0 --timeframes-rate-limit 1 --timeframes-rate-limit-ipcid 0 " #--infologger-mode \"stdout\""

#readCmd="o2-ctf-reader-workflow --copy-cmd no-copy --ctf-input ${inputfile} --delay 8 --loop 0 --onlyDet MFT --shm-segment-id 0 --shm-segment-size ${shmSize} ${logConfig} --allow-missing-detectors --condition-remap file://${localCcdbBaseDir}=${ccdbGeomAlignedPath} "
readCmd="o2-ctf-reader-workflow --copy-cmd \"alien_cp ?src file:?dst\" --remote-regex \"^alien:///alice/data/.+\" --ctf-input ${inputfile} --delay 8 --loop 0 --onlyDet MFT --shm-segment-id 0 --shm-segment-size ${shmSize} ${logConfig} --allow-missing-detectors --condition-remap file://${localCcdbBaseDir}=${ccdbGeomAlignedPath} "

recoOptions="MFTTracking.FullClusterScan=true;MFTTracking.LTFclsRCut=0.2;MFTTracking.trackmodel=2;MFTAlpideParam.roFrameLengthInBC=198;"
recoCmd="o2-mft-reco-workflow --shm-segment-id 0 --shm-segment-size ${shmSize} ${logConfig} --nThreads 1 --clusters-from-upstream --mft-track-writer --mft-cluster-writer --disable-mc --pipeline mft-tracker:1 --run-assessment --configKeyValues \""${recoOptions}"\" --condition-remap file://${localCcdbBaseDir}=${ccdbGeomAlignedPath} "

# Concatenate workflow
runCmd=" $readCmd "
runCmd+=" | $recoCmd"
runCmd+=" | o2-dpl-run  ${logConfig} --shm-segment-size ${shmSize} -b --run"
runCmd+=" | o2-dpl-run ${logConfig} --shm-segment-id 0 --shm-segment-size ${shmSize} -b --run > ctf2cltrack.log"

# List input files and command line
echo "======================="
echo "InputFiles on ${tmpInFile} :"
cat  ${tmpInFile}
echo "number of files :"
cat  ${tmpInFile} | wc -l
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

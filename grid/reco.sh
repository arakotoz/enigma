#!/usr/bin/env bash

## Ensure necessary files are present in the node work directory
geomFilename="o2sim_geometry.root"
geomFilenameAligned="o2sim_geometry-aligned.root"
ctfDict="ctf_dictionary.root"
mftDict="MFTdictionary.bin"
grpFile="o2sim_grp.root"

if [ ! -e "$geomFilename" ]; then
    echo "Cannot find $geomFilename"
    exit 1
fi

if [ ! -e "$geomFilenameAligned" ]; then
    echo "Cannot find $geomFilenameAligned"
    exit 1
fi

if [ ! -e "$ctfDict" ]; then
    echo "Cannot find $ctfDict"
    exit 1
fi


if [ ! -e "$mftDict" ]; then
    echo "Cannot find $mftDict"
    exit 1
fi

if [ ! -e "$grpFile" ]; then
    echo "Cannot find $grpFile"
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

severity="info"
logConfig="--severity ${severity} " #--infologger-mode \"stdout\""
shmSize=$(( 8 << 30 ))

readCmd="o2-ctf-reader-workflow --copy-cmd no-copy --ctf-input ${tmpInFile} --delay .1 --onlyDet MFT --shm-segment-size ${shmSize} ${logConfig} -b "

recoOptions="MFTTracking.forceZeroField=true; MFTTracking.FullClusterScan=true; MFTTracking.MinTrackPointsLTF=5; MFTTracking.MinTrackStationsLTF=4; MFTTracking.LTFclsRCut=0.2; MFTTracking.LTFseed2BinWin=3; MFTTracking.MinTrackStationsCA=4; MFTTracking.ROADclsRCut=0.04; MFTTracking.MinTrackPointsCA=5;MFTAlpideParam.roFrameLengthInBC=198;"
recoCmd="o2-mft-reco-workflow ${logConfig} --clusters-from-upstream --mft-cluster-writer --disable-mc --configKeyValues \""${recoOptions}"\""

assessCmd="o2-mft-assessment-workflow ${logConfig} --disable-mc -b"

# Concatenate workflow
runCmd=" $readCmd "
runCmd+=" | $recoCmd"
#runCmd+=" | $assessCmd"
runCmd+=" | o2-dpl-run  ${logConfig} -b --run"

# List input files and command line
echo "======================="
echo "InputFiles on ${tmpInFile} :"
cat  ${tmpInFile}
echo "======================="
echo "Running reconstruction command: "
echo "${runCmd} "
startTime=$(date +"%Y %m %d %H:%M:%S")
echo "Start "${startTime}
echo "======================="

## Runs reconstruction comand stored in $runCmd
echo | eval "${runCmd}" # "echo | " is a hack (to provide input stream to O2 workflows?)
exit "$?"

#!/usr/bin/env bash

# make sure that you are using O2 from CVMFS
# via for e.g. 
# > /cvmfs/alice.cern.ch/bin/alienv enter VO_ALICE@O2PDPSuite::nightly-20211026-1
# or 
# > /cvmfs/alice.cern.ch/bin/alienv enter VO_ALICE@O2Physics::nightly-20211026-1

myDir="/Users/andry/cernbox/alice/mft/pilotbeam/505713/test-reco"
#myDir="~/mycernbox/alice/mft/pilotbeam/505713/test-reco"
cd $myDir
pwd

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

xmlList="505713.xml"
inputfile="small-ctf-file-list.txt"

#makeFileList "${xmlList}" "${inputfile}" "${nodownload}"

severity="info"
logConfig="--severity ${severity} --fairmq-rate-logging 0 " #--infologger-mode \"stdout\""
shmSize=16000000000

readCmd="o2-ctf-reader-workflow --copy-cmd no-copy --ctf-input ${inputfile} --delay 8 --onlyDet MFT --shm-segment-size ${shmSize} ${logConfig} --allow-missing-detectors -b "
recoOptions="MFTTracking.forceZeroField=true;MFTTracking.FullClusterScan=true;MFTTracking.LTFclsRCut=0.2;MFTAlpideParam.roFrameLengthInBC=198;"
recoOptions+="NameConf.mDirGeom=${myDir};"
recoOptions+="NameConf.mDirCollContext=${myDir};"
recoOptions+="NameConf.mDirGRP=${myDir};"
recoCmd="o2-mft-reco-workflow --shm-segment-size ${shmSize} ${logConfig} --clusters-from-upstream --mft-track-writer --mft-cluster-writer --disable-mc --pipeline mft-tracker:1 --run-assessment --configKeyValues \""${recoOptions}"\""

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
exit "$?"

echo "======================="

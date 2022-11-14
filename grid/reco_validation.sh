#!/bin/bash
##################################################
validateout=$(dirname "$0")
validatetime=$(date)
error=0
if [ -z "$validateout" ]; then
   validateout="."
fi

cd "$validateout" || exit 1
validateworkdir=$(pwd)

echo "*******************************************************"
echo "* Executing validation script           *"
echo "* Time:    $validatetime "
echo "* Dir:     $validateout"
echo "* Workdir: $validateworkdir"
echo "* ----------------------------------------------------*"
ls -la ./
echo "* ----------------------------------------------------*"

##################################################

if [ ! -f stderr ]; then
   error=1
   echo "* ########## Job not validated - no stderr  ###"
   echo "Error = $error"
fi
parArch=$(grep -Ei "Cannot Build the PAR Archive" stderr)
segViol=$(grep -Ei "Segmentation violation" stderr)
segFault=$(grep -Ei "Segmentation fault" stderr)
glibcErr=$(grep -Ei '\*\*\* glibc detected \*\*\*' stderr)

if [ "$parArch" != "" ]; then
   error=1
   echo "* ########## Job not validated - PAR archive not built  ###"
   echo "$parArch"
   echo "Error = $error"
fi
if [ "$segViol" != "" ]; then
   error=1
   echo "* ########## Job not validated - Segment. violation  ###"
   echo "$segViol"
   echo "Error = $error"
fi
if [ "$segFault" != "" ]; then
   error=1
   echo "* ########## Job not validated - Segment. fault  ###"
   echo "$segFault"
   echo "Error = $error"
fi
if [ "$glibcErr" != "" ]; then
   error=1
   echo "* ########## Job not validated - *** glibc detected ***  ###"
   echo "$glibcErr"
   echo "Error = $error"
fi
if ! [ -f mfttracks.root ]; then
   error=1
   echo "Output file mfttracks.root not found. Job FAILED !"
   echo "Output file mfttracks.root not found. Job FAILED !" >>stderr
fi
if ! [ -f mftclusters.root ]; then
   error=1
   echo "Output file mftclusters.root not found. Job FAILED !"
   echo "Output file mftclusters.root not found. Job FAILED !" >>stderr
fi
if ! [ -f MFTAssessment.root ]; then
   error=1
   echo "Output file MFTAssessment.root not found. Job FAILED !"
   echo "Output file MFTAssessment.root not found. Job FAILED !" >>stderr
fi

if [ $error = 0 ]; then
   echo "* ----------------   Job Validated  ------------------*"
   echo "* === Logs std* will be deleted === "
   rm -f std*
fi
echo "* ----------------------------------------------------*"
echo "*******************************************************"
cd - || exit 1
exit $error

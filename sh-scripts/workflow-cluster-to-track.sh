#!/usr/bin/env bash

# make sure that you are using O2 from CVMFS
# via for e.g. 
# > /cvmfs/alice.cern.ch/bin/alienv enter VO_ALICE@O2PDPSuite::nightly-20211026-1
# or 
# > /cvmfs/alice.cern.ch/bin/alienv enter VO_ALICE@O2Physics::nightly-20211026-1

function groupfiles () {
    # create a subgroup of files for reco
    # from a list of files, obtained for e.g. with:
    #      alien_ls -c /alice/data/2021/OCT/504980/raw/*.tf > tflist
    # or
    #      alien_ls -c /alice/data/2021/OCT/504980/raw/*.root > ctflist
    # use:
    #      FilesPerGroup=15 groupfiles ctflist
    # alternatively, use the bash program split instead:
    #      split -l 15 -d ctflist group_
    FilesPerGroup=${FilesPerGroup:-"50"}
    GroupCounter=0
    FileCounter=0
    rm group_*
    for file in $(cat $@)
    do
	echo ${file} >> group_${GroupCounter}
	((FileCounter++))

	if [ ${FileCounter} -eq ${FilesPerGroup} ]
	then
            ((GroupCounter++))
            FileCounter=0
	fi          
    done
}
export -f groupfiles


function go_cluster_to_track() {
    # read cluster file(s), run reco workflow on these clusters, output tracks ROOT files
    # use:
    #     go_cluster_to_track file_with_clusterfilelist
    inputfile=${@}
    MINCLUSTERS=${MINCLUSTERS:-"4"}
    ROADRADIUS=${ROADRADIUS:-"0.5"}
    prefix=${@}_${MINCLUSTERS}cls_${ROADRADIUS}road

    o2-mft-cluster-reader-workflow --shm-segment-size 15000000000 --mft-cluster-infile ${inputfile}  | o2-mft-reco-workflow --mft-track-writer "--outfile mfttracks_${prefix}.root" --mft-cluster-writer "--outfile /dev/null" --shm-segment-size 15000000000 --clusters-from-upstream --disable-mc --configKeyValues "MFTTracking.FullClusterScan=true;  MFTTracking.MinTrackPointsLTF=${MINCLUSTERS}; MFTTracking.MinTrackStationsLTF=1; MFTTracking.LTFclsRCut=${ROADRADIUS}; MFTTracking.LTFseed2BinWin=10; MFTTracking.MinTrackStationsCA=1; MFTTracking.ROADclsRCut=${ROADRADIUS}; MFTTracking.MinTrackPointsCA=${MINCLUSTERS};" -b --run > cluster2track_${prefix}.log 
}
export -f go_cluster_to_track


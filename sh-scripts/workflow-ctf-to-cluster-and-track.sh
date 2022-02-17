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

function go_ctf_to_cluster_track() {
    # read CTF file(s), run reco workflow on these files, output tracks and cluster ROOT files
    # with a list of CTF files obtained from for e.g.
    #      alien_ls -c /alice/data/2021/OCT/504980/raw/*.root > ctflist
    # use:
    #      go_ctf_to_cluster_track ctflist
    inputfile=${@}
    MINCLUSTERS=${MINCLUSTERS:-"4"}
    ROADRADIUS=${ROADRADIUS:-"0.1"}
    prefix=${@}_${MINCLUSTERS}cls_${ROADRADIUS}road

    o2-ctf-reader-workflow --copy-cmd "alien_cp ?src ?dst" --remote-regex "/alice/data/2021/.+" --onlyDet MFT --delay 0.1 --ctf-input ${inputfile}  | o2-mft-reco-workflow --mft-track-writer "--outfile mfttracks_${prefix}.root" --mft-cluster-writer "--outfile mftclusters_${prefix}.root" --shm-segment-size 15000000000 --clusters-from-upstream --disable-mc --configKeyValues "MFTTracking.FullClusterScan=true;  MFTTracking.MinTrackPointsLTF=${MINCLUSTERS}; MFTTracking.MinTrackStationsLTF=1; MFTTracking.LTFclsRCut=${ROADRADIUS}; MFTTracking.LTFseed2BinWin=6; MFTTracking.MinTrackStationsCA=1; MFTTracking.ROADclsRCut=${ROADRADIUS}; MFTTracking.MinTrackPointsCA=${MINCLUSTERS};" -b --run > reco_${prefix}.log && mv ${inputfile} done_${inputfile}
}

export -f go_ctf_to_cluster_track



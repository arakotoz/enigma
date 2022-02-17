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
    # the command used here was:
    #      split -l 30 -a 3 -d raw_0100_tf a_raw_0100_tf_
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


function go-tf-to-cluster-track() {
    # read TF file(s), run reco workflow on these files, output tracks and cluster ROOT files
    # with a list of TF files obtained from for e.g.
    #      alien_ls -c /alice/data/2021/OCT/504980/raw/*.tf > tflist
    # use:
    #      go-tf-to-cluster-track tflist
    inputfile=${@}
    prefix=${@}
    RECO_OPTIONS="MFTTracking.FullClusterScan=false; MFTTracking.MinTrackPointsLTF=5; MFTTracking.MinTrackStationsLTF=4; MFTTracking.LTFclsRCut=0.2; MFTTracking.LTFseed2BinWin=3; MFTTracking.MinTrackStationsCA=4; MFTTracking.ROADclsRCut=0.04; MFTTracking.MinTrackPointsCA=5;MFTAlpideParam.roFrameLengthInBC=198;"
    OUTDIR=${prefix}-outdir
    LOG=${OUTDIR}/tf-to-cluster-track.log

    mkdir -p ${OUTDIR}/raw
    mkdir -p finished
    echo 
    echo "-----------------------------"
    echo "go-tf-to-cluster-track() ${inputfile}"
    echo "go-tf-to-cluster-track() ${inputfile}" > ${LOG}
    echo
    START_TIME=$(date +"%Y %m %d %H:%M:%S")
    echo "Start "${START_TIME}
    echo "Start "${START_TIME} >> ${LOG}
    echo >> ${LOG}
    echo "o2-raw-tf-reader-workflow -b --copy-cmd \"alien_cp ?src file://?dst\" --remote-regex \"^/alice/data/.+\" --shm-segment-size 15000000000 --input-data "${inputfile}" --onlyDet MFT |  o2-itsmft-stf-decoder-workflow -b --runmft --shm-segment-size 15000000000 --nthreads 4 --digits --no-clusters --no-cluster-patterns | o2-mft-reco-workflow -b --run --mft-track-writer \"--outfile "${OUTDIR}"/mfttracks.root\" --mft-cluster-writer \"--outfile "${OUTDIR}"/mftclusters.root\" --shm-segment-size 15000000000 --digits-from-upstream --disable-mc --configKeyValues \""${RECO_OPTIONS}"\"" >> ${LOG}
    echo >> ${LOG}
    o2-raw-tf-reader-workflow -b --copy-cmd "alien_cp ?src file://?dst" --remote-regex "^/alice/data/2021/.+" --shm-segment-size 15000000000 --input-data ${inputfile} --onlyDet MFT |  o2-itsmft-stf-decoder-workflow  -b --runmft --shm-segment-size 15000000000 --nthreads 4 --digits --no-clusters --no-cluster-patterns | o2-mft-reco-workflow -b --run --mft-track-writer "--outfile ${OUTDIR}/mfttracks.root" --mft-cluster-writer "--outfile ${OUTDIR}/mftclusters.root" --shm-segment-size 15000000000 --digits-from-upstream --disable-mc --configKeyValues ${RECO_OPTIONS} >> ${LOG}
    echo
    mv -v ${inputfile} finished/${inputfile}
    mv rawdump*.raw ${OUTDIR}/raw/.
    echo
    echo >> ${LOG}
    STOP_TIME=$(date +"%Y %m %d %H:%M:%S")
    echo "Stop "${STOP_TIME}
    echo "Stop "${STOP_TIME} >> ${LOG}
    echo "-----------------------------"
}
export -f go-tf-to-cluster-track

function go-tf-to-cluster-track-b0() {
    # read TF file(s), run reco workflow with forced B=0 on these files, output tracks and cluster ROOT files
    # with a list of TF files obtained from for e.g.
    #      alien_ls -c /alice/data/2021/OCT/504980/raw/*.tf > tflist
    # use:
    #      go-tf-to-cluster-track-b0 tflist
    inputfile=${@}
    prefix=${@}
    RECO_OPTIONS="MFTTracking.forceZeroField=true; MFTTracking.FullClusterScan=false; MFTTracking.MinTrackPointsLTF=5; MFTTracking.MinTrackStationsLTF=4; MFTTracking.LTFclsRCut=0.2; MFTTracking.LTFseed2BinWin=3; MFTTracking.MinTrackStationsCA=4; MFTTracking.ROADclsRCut=0.04; MFTTracking.MinTrackPointsCA=5;MFTAlpideParam.roFrameLengthInBC=198;"
    OUTDIR=${prefix}-outdir
    LOG=${OUTDIR}/tf-to-cluster-track-b0.log
    SHMSIZE=$(( 8 << 30 ))
    NTHREADS=3

    mkdir -p ${OUTDIR}/raw
    mkdir -p finished
    echo 
    echo "-----------------------------"
    echo "go-tf-to-cluster-track-b0() ${inputfile}"
    echo "go-tf-to-cluster-track-b0() ${inputfile}" > ${LOG}
    echo
    START_TIME=$(date +"%Y %m %d %H:%M:%S")
    echo "Start "${START_TIME}
    echo "Start "${START_TIME} >> ${LOG}
    echo >> ${LOG}
    echo "o2-raw-tf-reader-workflow -b --copy-cmd \"alien_cp ?src file://?dst\" --remote-regex \"^/alice/data/.+\" --shm-segment-size "${SHMSIZE}" --input-data "${inputfile}" --onlyDet MFT |  o2-itsmft-stf-decoder-workflow -b --runmft --shm-segment-size "${SHMSIZE}" --nthreads "${NTHREADS}" --digits --no-clusters --no-cluster-patterns | o2-mft-reco-workflow -b --run --mft-track-writer \"--outfile "${OUTDIR}"/mfttracks.root\" --mft-cluster-writer \"--outfile "${OUTDIR}"/mftclusters.root\" --shm-segment-size "${SHMSIZE}" --digits-from-upstream --disable-mc --configKeyValues \""${RECO_OPTIONS}"\"" >> ${LOG}
    echo >> ${LOG}
    o2-raw-tf-reader-workflow -b --copy-cmd "alien_cp ?src file://?dst" --remote-regex "^/alice/data/2021/.+" --shm-segment-size ${SHMSIZE} --input-data ${inputfile} --onlyDet MFT |  o2-itsmft-stf-decoder-workflow  -b --runmft --shm-segment-size ${SHMSIZE} --nthreads ${NTHREADS} --digits --no-clusters --no-cluster-patterns | o2-mft-reco-workflow -b --run --mft-track-writer "--outfile ${OUTDIR}/mfttracks.root" --mft-cluster-writer "--outfile ${OUTDIR}/mftclusters.root" --shm-segment-size ${SHMSIZE} --digits-from-upstream --disable-mc --configKeyValues ${RECO_OPTIONS} >> ${LOG}
    echo
    mv -v ${inputfile} finished/${inputfile}
    mv rawdump*.raw ${OUTDIR}/raw/.
    echo
    echo >> ${LOG}
    STOP_TIME=$(date +"%Y %m %d %H:%M:%S")
    echo "Stop "${STOP_TIME}
    echo "Stop "${STOP_TIME} >> ${LOG}
    echo "-----------------------------"
}
export -f go-tf-to-cluster-track-b0
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

function go_tf_to_digits() {
    # read TF file(s), run decoder on these files, output digit ROOT files
    # with a list of TF files obtained from for e.g.
    #      alien_ls -c /alice/data/2021/OCT/504980/raw/*.tf > tflist
    # use:
    #      go_tf_to_digits tflist
    inputfile=${@}
    prefix=${@}

    o2-raw-tf-reader-workflow --copy-cmd "alien_cp ?src ?dst" --remote-regex "/alice/data/2021/.+"  --shm-segment-size 15000000000 --input-data ${inputfile} --onlyDet MFT |  o2-itsmft-stf-decoder-workflow  --shm-segment-size 15000000000 --nthreads 8 --runmft --digits --no-clusters --no-cluster-patterns | o2-itsmft-digit-writer-workflow  --disable-mc --outfile mftdigits_${prefix}.root -b --runmft > decode_and_digits_TF${prefix}.log && mv ${inputfile} done_${inputfile}
}
export -f go_tf_to_digits



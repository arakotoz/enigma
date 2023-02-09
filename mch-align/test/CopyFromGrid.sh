#!/bin/bash
TARGET_DIR=$1
LOCAL_DIR=$2
SET=$3
RUN=$4
PASS=$5
CTF=$6

NEW_DIR=$SET"_"$RUN"_"$PASS"_"$CTF

mkdir -p $LOCAL_DIR/$NEW_DIR
cd $LOCAL_DIR/$NEW_DIR

echo "Start copying from $TARGET_DIR ..."
CTF_LIST=`alien_ls $TARGET_DIR | grep o2_ctf_*`
echo "Found jobs:"
echo $CTF_LIST

for CTF in $CTF_LIST
do
	mkdir -p $CTF
	echo "Processing $CTF"
	for ITEM in log_archive.zip mchtracks.root QC.root stderr.log stdout.log
	do
		echo "Copying $ITEM ..."
		alien_cp alien:$TARGET_DIR/$CTF/$ITEM file:$LOCAL_DIR/$NEW_DIR/$CTF/$ITEM
		echo "$ITEM copied."
	done

done

hadd -f mchtracks.root o2_ctf_*/mchtracks.root
#!/bin/bash

CHR="$1"
MAJOR_NAME="_$2"
MAJOR_CODE="$3"
MINOR_CODE="$4"

#---------------#
#  input files  #
#---------------#
RMSK_IN="/scratch/dpham4/PI/data/rmsk.txt.gz"
CHIP_SEQ_IN="/home/mcb/blanchem/wgEncodeRegTfbsClustered/${MAJOR_CODE:0:6}$MAJOR_NAME.bed"
MAJOR_IN="/scratch/blanchem/$CHR/sites/sites.$MAJOR_CODE.gz"
MINOR_IN="/scratch/blanchem/$CHR/sites/sites.$MINOR_CODE.gz"

#----------------#
#  output files  #
#----------------#
RMSK_OUT="/scratch/dpham4/PI/data/$CHR/rmsk.txt"
CHIP_SEQ_OUT="/scratch/dpham4/PI/data/$CHR/chip_seq_${MAJOR_CODE:0:6}.txt"
MAJOR_OUT="/scratch/dpham4/PI/data/$CHR/$MAJOR_CODE.txt"
MINOR_OUT="/scratch/dpham4/PI/data/$CHR/$MINOR_CODE.txt"

#--------#
#  main  #
#--------#
if [ ! -f $RMSK_OUT ]
then
    zcat $RMSK_IN | grep -w $CHR > $RMSK_OUT
    if [ "$?" != "0" ]
    then
        rm -f $RMSK_OUT
        exit 1
    fi
fi

if [ ! -f $CHIP_SEQ_OUT ]
then
    grep -w $CHR $CHIP_SEQ_IN > $CHIP_SEQ_OUT
    if [ "$?" != "0" ]
    then
        rm -f $CHIP_SEQ_OUT
        exit 1
    fi
fi

if [ ! -f $MAJOR_OUT ]
then
    zcat $MAJOR_IN | grep -w '^0' > $MAJOR_OUT
    if [ "$?" != "0" ]
    then
        rm -f $MAJOR_OUT
        exit 1
    fi
fi

if [ ! -f $MINOR_OUT ]
then
    zcat $MINOR_IN | grep -w '^0' > $MINOR_OUT
    if [ "$?" != "0" ]
    then
        rm -f $MINOR_OUT
        exit 1
    fi
fi
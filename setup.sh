#!/bin/bash

CHR="$1"
TF1_NAME="_$2"
TF1_CODE="$3"
TF2_CODE="$4"

#---------------#
#  input files  #
#---------------#
RMSK_IN="/scratch/dpham4/PI/data/rmsk.txt.gz"
CHIP_IN="/home/mcb/blanchem/wgEncodeRegTfbsClustered/${TF1_CODE:0:6}$TF1_NAME.bed"
TF1_IN="/scratch/blanchem/$CHR/sites/sites.$TF1_CODE.gz"
TF2_IN="/scratch/blanchem/$CHR/sites/sites.$TF2_CODE.gz"

#----------------#
#  output files  #
#----------------#
RMSK_OUT="/scratch/dpham4/PI/data/$CHR/rmsk.txt"
CHIP_OUT="/scratch/dpham4/PI/data/$CHR/chip_seq_${TF1_CODE:0:6}.txt"
TF1_OUT="/scratch/dpham4/PI/data/$CHR/$TF1_CODE.txt"
TF2_OUT="/scratch/dpham4/PI/data/$CHR/$TF2_CODE.txt"

#--------#
#  main  #
#--------#
if [ ! -d "/scratch/dpham4/PI/data/$CHR" ]
then
    mkdir "/scratch/dpham4/PI/data/$CHR"
fi

if [ ! -f $RMSK_OUT ]
then
    zcat $RMSK_IN | grep -w $CHR > $RMSK_OUT
    python min_rmsk.py $CHR
    if [ "$?" != "0" ]
    then
        rm -f $RMSK_OUT
        exit 1
    fi
fi

if [ ! -f $CHIP_OUT ]
then
    grep -w $CHR $CHIP_IN > $CHIP_OUT
    if [ "$?" != "0" ]
    then
        rm -f $CHIP_OUT
        exit 2
    fi
fi

if [ ! -f $TF1_OUT ]
then
    zcat $TF1_IN | grep -w '^0' > $TF1_OUT
    if [ "$?" != "0" ]
    then
        rm -f $TF1_OUT
        exit 3
    fi
fi

if [ ! -f $TF2_OUT ]
then
    zcat $TF2_IN | grep -w '^0' > $TF2_OUT
    if [ "$?" != "0" ]
    then
        rm -f $TF2_OUT
        exit 4
    fi
fi
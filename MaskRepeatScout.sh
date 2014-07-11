#!/bin/bash
############################################
# Nathan Figueroa - Miami University
# April 01 2013
#
# Simple script for automating all of 
# RepeatScout's procedures
#
# Takes one parameter: Sequence file name
#
############################################

SEQUENCE=$1
OUTPUT=$2
PAS=$3
RMDIR=$4

FREQ=${OUTPUT//".fa"/".freq.fa"}
FILTERED=${OUTPUT//".fa"/".filtered.fa"}
FINAL=${OUTPUT//".fa"/".final.fa"}

#RepeatMasker $SEQUENCE -lib $FILTERED -pa $PAS -dir $RMDIR
RepeatMasker $SEQUENCE -lib $OUTPUT -pa $PAS -dir $RMDIR
cat $OUTPUT | perl filter-stage-2.prl --cat=$RMDIR/$SEQUENCE.out --thresh=10 > $FINAL
RepeatMasker $SEQUENCE -lib $FINAL -pa $PAS -dir $RMDIR

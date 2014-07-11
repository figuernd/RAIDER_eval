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
L=$2
OUTPUT=$3
MIN=$4

FREQ=${OUTPUT//".fa"/".freq.fa"}
FILTERED=${OUTPUT//".fa"/".filtered.fa"}

./build_lmer_table -min $MIN -l $L -sequence $SEQUENCE -freq $FREQ
./RepeatScout -sequence $SEQUENCE -freq $FREQ -output $OUTPUT
#cat $OUTPUT | perl filter-stage-1.prl > $FILTERED


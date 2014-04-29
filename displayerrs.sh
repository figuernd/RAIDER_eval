#!/bin/bash
echo "Simulated Sequences:"
head -5 *sim.*.fa.batch.e*
echo
echo "RAIDER:"
head -5 raider.e*
echo
echo "Consensus Sequences:"
head -5 *.sim.*.consensus.fa.batch.e*
echo
echo "REPEATMASKER:"
head -5 repeatmasker.e*
echo 
echo "STATS:"
head -5 stats.e*

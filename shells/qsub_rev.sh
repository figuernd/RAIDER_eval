#!/bin/bash 
c=2
rep=5
inc=5
cnt=15
while [ $cnt -lt 100 ]
        do
            qsub -v c=$c,cnt=$cnt,reps=$rep rev.batch
            sleep 10m
        cnt=$(($cnt+$inc))
done


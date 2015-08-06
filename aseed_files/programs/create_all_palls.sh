#!/bin/bash

# create_all_palls.sh: Wrapper to create seed files of all combinations of weight
# and length between specified min and max weight (minw, maxw) and specified min
# and max length (minl, maxl). Arguments: minw maxw minl maxl
# by Carly Schaeffer


aseed_dir='../aseeds/'
dir=${aseed_dir}L_${minl}_${maxl}_W_${minw}_${maxw}
currl=${minl}
currw=${minw}

mkdir ${dir}

echo "min weight= ${minw} max weight= ${maxw}"
echo "min length= ${minl} max length= ${maxl}"

while [[ $currl -le $maxl ]]; do
    currw=${minw}
    echo "currl= ${currl}"
    while [[ $currw -le $maxw && $currw -le $currl ]]; do
        echo "currw= ${currw}"
        python3 create_all_pals.py -l ${currl} -w ${currw} ${dir}/L${currl}_W${currw}
        let currw=currw+1
    done
    let currl=currl+1
done

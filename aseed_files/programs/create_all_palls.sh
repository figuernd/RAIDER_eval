#!/bin/bash


#pass in minw maxw minl maxl

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


#!/bin/bash

cd ~/karro/RAIDER_eval/
##comp=comparisons_overlap3




seed=$(head -n 1 ${sfile})
weight=$(grep -o "1" <<< "$seed" | wc -l)
L=${#seed} 
basedata=${data##*/}
id=$(basename $sfile .txt)

fname=${dir}/F${freq}_${id}



##data=hg18.chr22.fa


if [ -z $M ]; then
    mem=""
else
    mem="--mem"
    fname=${fname}_MEM
fi

echo "Output directory: ${fname}"
echo "Data: ${basedata}"
echo "Seed file: ${sfile}"

if [ ${type} == "chrom_sim" ]; then
    fname=${dir}/sim/F${freq}_${id}
    if [ -z $T ]; then
        t_str=""
        T=1
    else
        t_str="-t ${T}"
    fi
    if [ -z $ST ]; then
        st_str=""
        ST=0
    else
        st_str="--st ${ST}"
    fi
    if [ -z $RNG ]; then
        rng_str=""
        RNG=None
    else
        rng_str="--rng_seed ${RNG}"
    fi
    fname=${fname}_T${T}_ST${ST}_RNG${RNG}
    
    if [ -z $FAM ]; then
        fam_str=""
    else
        fam_str="--family_file ${FAM}"
        fname=${fname}_FAM
    fi


    ./RAIDER_eval.py -R --R2 --RS --nuke -r $fname -f ${freq} --sf ${sfile} --age 1 $mem --timing_jobs chrom_sim -k 0 ${t_str} ${st_str} ${rng_str} ${fam_str} ${data} 
    mv ${fname}/SOURCE_DATA/ sim_data/SIM_T${T}_ST${ST}_RNG${RNG}_FAM${FAM}
else
    fname=${dir}/seq/F${freq}_${id}
    ./RAIDER_eval.py -R --R2 --RS --nuke -r $fname -f ${freq} --sf ${sfile} --age 1 $mem --timing_jobs seq_files ${data} 
fi

##dir=${fname} ./cleanup.sh

./doSortedAnalysis.R -f ${fname}/stats.txt --formula=~-tpr -v -a -m

#!/bin/bash

# rev_sf.sh : shell script that automates the process of naming RAIDER_eval output
#  directory, running the RAIDER_eval program, cleaning up the output directory, and
#  creating a formatted statistics output file.
# Arguments:
#   $dir - base directory for output dir, to be modified
#   $sfile OR $seed - depending on which, adjusts naming of output dir to indicate
#   $freq
#   $data
#   $type (default is seq_files)
# By Carly Schaeffer


cd ~/karro/RAIDER_eval/ # need to change this to your RAIDER_eval directory

if [ -n "$sfile" ]; then
    id=$(basename $sfile .txt)
    seed_str="--sf ${sfile}"
else
    L=${#seed}
    weight=$(grep -o "1" <<< "$seed" | wc -l)
    id=L${L}_W${weight}_S${seed}
    seed_str="-s ${seed}"
fi

echo $id
basedata=${data##*/}
fname=${dir}/F${freq}_${id}

module load python-3.3.3

if [ -z $opts ]; then
  opts=""
fi
#if [ -v $PS ]; then
#    opts="${opts} --ps"
#    opts_dir_str=${opts_dir_str}_PROSPLIT
#fi


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


    ./RAIDER_eval.py --R2 -R --nuke -r $fname -f ${freq} ${seed_str} --age 1  --timing_jobs chrom_sim ${t_str} ${st_str} ${rng_str} ${fam_str} ${data}
    mv ${fname}/SOURCE_DATA/ sim_data/SIM_T${T}_ST${ST}_RNG${RNG}_FAM${FAM}

else
    fname=${dir}/seq/F${freq}_${id}
    ./RAIDER_eval.py --R2 --RA --nuke ${opts} --no -r $fname -f ${freq} ${seed_str} --timing_jobs seq_files ${data}
fi

dir=${fname} ./cleanup.sh

./doSortedAnalysis.R -f ${fname}/stats.txt -v -a -m -c

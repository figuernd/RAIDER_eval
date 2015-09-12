rm nohup.out
nohup python testing_pipeline.py -r BEST --sf BEST.txt --mn -f 2 --pa 2 -R --R2 --RS sim_data/hg38.chr22.fa sim_data/hg38.chr11.fa hg38.chr1.fa ce10.chrV.fa mm10.chr19.fa > nohup.out&

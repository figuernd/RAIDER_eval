rm nohup.out
nohup python testing_pipeline.py -r SEEDS12 --sf seeds12.txt -f 2 3 --mn --pa 2 --R2 sim_data/hg38.chr22.fa > nohup.out&

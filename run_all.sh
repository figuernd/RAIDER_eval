rm nohup.out
nohup python3.3 testing_pipeline.py -r SEEDS1 --sf seeds1.txt -f 2 --mn --pa 2 --R2 --RS sim_data/ce10.chrV.fa sim_data/hg38.chr22.fa > nohup.out&

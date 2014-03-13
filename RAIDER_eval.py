#!/usr/bin/python3.3
import sys
import subprocess
import os.path
import os.path
import argparse



def parse_params(args):
    parser = argparse.ArgumentParser(description = "Evaluate RAIDER against RepeatScout")
    parser.add_argument('-f', type = int, help = "e.r. occurance threshold", default = 2)
    parser.add_argument('-s', '--seed', help = "spaced seed string", default = "1111011110111101111")
    parser.add_argument('seq_file', help = "File containing genomic sequence")
    parser.add_argument('output_dir', help = "Output directory")
    return parser.parse_args(args)
    


if __name__:
    args = parse_params(sys.argv[1:])
    if not os.path.exists('./%s' % (args.output_dir)):
        print("HERE")
        os.makedirs('./%s' % (args.output_dir))

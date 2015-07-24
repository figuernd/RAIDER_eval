#!/software/python/3.3.3/bin/python3.3

import itertools
import argparse


def generate_seed_file(seed_file, length, weight): 
    f = open(seed_file, 'w')
    for bits in itertools.combinations(range((length - 2)//2), (weight - 2)//2):
        s = ['0'] * length
        s[0] = '1'
        s[length - 1] = '1'
        for bit in bits:
            s[bit + 1] = '1'
            s[length - bit - 2] = '1'
        f.write(''.join(s) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Generate All Palindromic Seeds With Specified Length and Weight.")
    parser.add_argument('-l', '--length', type = int, help = "Seed length", default = 23)
    parser.add_argument('-w', '--weight', type = int, help = "Seed weight", default = 13)
    parser.add_argument("seed_file", help = "Seed output file")
    args = parser.parse_args()
    generate_seed_file(args.seed_file, args.length, args.weight)

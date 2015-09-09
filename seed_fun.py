import random
import re

def compact_seed(s):
    while True:
        r = re.search("(111+)", s)
        if r:
            w = r.group(1)
            s = re.sub("1"*len(w), "1^{" + str(len(w)) + "}", s)
        else:
            break

    while True:
        r = re.search("(000+)", s)
        if r:
            w = r.group(1)
            s = re.sub("0"*len(w), "0^{" + str(len(w)) + "}", s)
        else:
            break
    return s

def weight(key):
    return sum([1 for s in key if s == "1"])

def check_weights(file):
    for line in open(file):
        seed = line.rstrip()
        w = weight(seed)
        print(len(seed), w, float(w)/len(seed))

def gen_seed(L):
    s = ""
    t = 1
    for i in L:
        s += str(t) + "^{" + str(i) + "}"
        t = 1 - t

    return s + "\n"

def random_seed(w,p):
    """generate a random seed of weight w, with each character
    being a 0 with probability p."""
    c = 1
    s = "1"
    while c < w-1:
        if random.uniform(0,1) < p:
            s = s + "0"
        else:
            s = s + "1"
            c = c + 1
    return s + "1"


def random_seed2(w,l):
    """Generate a random seed of length l and weight w"""
    L = ["1"]*(w-2) + ["0"]*(l-w)
    random.shuffle(L)
    return "1" + "".join(L) + "1"

def random_seed3(w,l):
    """Generate a random palindomic seed o length l and weight w"""
    L = ["1"]*(w//2 - 1) + ["0"]*((l-w)//2)
    random.shuffle(L)
    return "1" + "".join(L + L[::-1]) + "1"

def generate_list(w, p_lower, p_upper, p_step, n):
    debug = 0
    M = {}

    p = p_lower
    while p < p_upper:
        l = int(round(float(w)/p))
        for k in range(n):
            num_tries = 0
            while True:
                s = random_seed3(w,l)
                if s not in M:
                    M[s] = True
                    break
                num_tries += 1
                assert(num_tries < 20)
            debug += 1
        p += p_step
    print(debug)
    return sorted(M.keys(), key=len)
        
def main():
    fp = open("seeds11.txt", "w")
    L = generate_list(40, 0.80, 0.801, 0.025, 500)
    fp.write("\n".join([compact_seed(s) for s in L]))
    fp.close()

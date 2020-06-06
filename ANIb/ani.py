import argparse
import os
from multiprocessing import Pool


def parse_args():
    parser = argparse.ArgumentParser(description='ANI all vs all computation')
    parser.add_argument('genome_list', type=str, help='A file containing all the genome to compare. One genome per line.')
    parser.add_argument('--cores', '-c', type=int, default=8, help='The number of cores to use')
    
    args = parser.parse_args()

    return args

def parse_list(file):
    root = "/".join(file.split("/")[:-1])
    files = []
    with open(file) as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            files.append(f"{root}/{line}")

    return files

def compute_ANI(genome1, genome2):
    exec_path = '/'.join(__file__.split('/')[:-1])
    if len(exec_path) == 0:
        exec_path = '.'
    stream = os.popen(f'ruby {exec_path}/ani.rb --seq1 {genome1} --seq2 {genome2}')
    output = stream.read()

    last_line = output.split("\n")[2].strip()
    last_line = last_line.split(": ")[1]
    last_line = last_line.split(' ')[0]
    ANI = float(last_line[:-1])

    return ANI

def compute_all_ani(lst):
    for idx, file in enumerate(lst):
        for file2 in lst[idx+1:]:
            ANI = compute_ANI(file, file2)
            out1 = file.split('/')[-1]
            out2 = file2.split('/')[-1]
            print(f"{out1}\t{out2}\t{ANI}")

def compute_parallel_ani(lst, cores):
    pairs = []
    for idx, file in enumerate(lst):
        for file2 in lst[idx+1:]:
            pairs.append((file, file2))

    p = Pool(cores)
    results = p.starmap(compute_ANI, pairs)

    for idx, pair in enumerate(pairs):
        out1 = pair[0].split('/')[-1]
        out2 = pair[1].split('/')[-1]
        print(f"{out1}\t{out2}\t{results[idx]}")

def main():
    args = parse_args()
    f_lst = parse_list(args.genome_list)
    compute_parallel_ani(f_lst, args.cores)


if __name__ == "__main__":
    main()

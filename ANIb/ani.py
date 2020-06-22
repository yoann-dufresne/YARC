import argparse
import os
from multiprocessing import Pool


def parse_args():
    parser = argparse.ArgumentParser(description='ANI all vs all computation')
    parser.add_argument('genome_list', type=str, help='A file containing all the genome to compare. One genome per line.')
    parser.add_argument('--cores', '-c', type=int, default=8, help='The number of cores to use')
    parser.add_argument('--queries', '-q', help='A file list containing the names of the genomes to compare with the whole genome list. If not specified, all vs all is performed.')
    
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

def compute_parallel_ani(file_pairs, cores):
    p = Pool(cores)
    results = p.starmap(compute_ANI, file_pairs)

    for idx, pair in enumerate(file_pairs):
        out1 = pair[0].split('/')[-1]
        out2 = pair[1].split('/')[-1]
        print(f"{out1}\t{out2}\t{results[idx]}")

def main():
    args = parse_args()
    f_lst = parse_list(args.genome_list)
    pairs = []
    if args.queries is None:
        compute_parallel_ani(combinations(f_lst,2), args.cores)
    else:
        q_lst = parse_list(args.queries)
        compute_parallel_ani(product(f_lst, q_lst), args.cores)


if __name__ == "__main__":
    main()

import paf_parser as pp
import sys
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from pyfaidx import Fasta


def main():
    args = parse_arguments()
    
    with open(args.paf) as paf:
        extract_insertions(paf, args.reads)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--paf', "-p", required=True, help="Pass a paf file to parse the alignment using cs string")
    parser.add_argument('--reads', "-r", required=True, help="A fasta file containing the reads mapped onto a reference?.")
    return parser.parse_args()


def extract_insertions(paf, fasta):
    reads = Fasta(fasta)
    
    for idx, match in enumerate(pp.parse(paf)):
        if idx % 10000 == 0:
            print(idx, file=sys.stderr)
            
        insert_sizes = get_max_insertion_positions(match["cs"])
        if len(insert_sizes) > 0:
            read_name = match["query"]["name"]

            for _, _, start, end in insert_sizes:
                subsequence = None

                # Offset the match position regarding orientation
                if match["orientation"] == "+":
                    start += match["query"]["start"]
                    end += match["query"]["start"]
                    subsequence = reads[read_name][start:end]
                else:
                    save = start
                    start = match["query"]["end"] - end
                    end = match["query"]["end"] - save
                    subsequence = Seq(str(reads[read_name][start:end]))
                    subsequence = subsequence.reverse_complement()
                    subsequence = str(subsequence)

                # Output fasta
                print(f">{read_name};offset={start};size={end-start}\n{subsequence}")



def get_max_insertion_positions(sequence, min_size=1000):
    insertions = []

    align_idx = 0
    current_align_idx = 0
    current_idx = 0
    current_ins = 0
    insertion_open = False

    for idx, tupl in enumerate(sequence):
        operation, value = tupl
            
        # Extends or close a current insertion detected
        if insertion_open:
            if operation == '*':
                current_ins += 1
            elif operation == '+':
                current_ins += len(value)
            else:
                insertion_open = False
                if current_ins >= min_size:
                    insertions.append((current_ins, current_idx, current_align_idx, align_idx))
        # Try to open an insertion
        else:
            if operation == '*':
                current_ins = 1
                current_idx = idx
                current_align_idx = align_idx
                insertion_open = True
            elif operation == '+':
                current_ins = len(value)
                current_idx = idx
                current_align_idx = align_idx
                insertion_open = True

        # Update the position in the sequence
        if operation == ":":
            align_idx += value
        elif operation == "*":
            align_idx += 1
        elif operation == "+":
            align_idx += len(value)

    if insertion_open and current_ins >= min_size:
        insertions.append((current_ins, current_idx, current_align_idx, align_idx))

    return insertions


if __name__ == "__main__":
    main()

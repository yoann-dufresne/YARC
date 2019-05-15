#!/usr/bin/env python3

import sys
import argparse
import csv
import random

from Bio import SeqIO
from Bio.Seq import Seq
from metropolis_hasting import MetropolisHasting1D


def parse_arguments():
    # Parse the command line
    parser = argparse.ArgumentParser(description="Process a reads FASTA file to integrate errors based on a PacBio model")
    parser.add_argument('--length_distribution', '-l', required=True, help="CSV file containing a read length distribution. The columns X and Y must be present where X is the read length and Y the number of reads with this length.")
    parser.add_argument('--depth', '-d', type=int, required=True, help="Sequencing depth. Without subsampling, this sequencing depth is perfect (ie each nucleotide of the ref is present exactly d time in the reads).")
    parser.add_argument('--sub-sample', '-s', type=float, default=1.0, help="Sub-sample ratio. Depending of this ratio, a bigger sequencing depth will be computed. From this bigger pool of reads, the software randomly sub sample to generate an average depth d.")
    parser.add_argument('--outprefix', '-o', help="The file prefix where the reads will be saved. If not specified, outputed on stdout (and stderr for paired end)")
    parser.add_argument('--reference', '-r', help="The fasta file containing the reference used to generate reads")
    parser.add_argument("--paired", "-p", action='store_true', help="If set, generate a pair of read files instead of one. See the gap argument for details on gap size.")
    parser.add_argument("--gap_distribution", '-g', help="The gap distribution file used for paired end reads. If not set, the gap will be 0.")
    parser.add_argument('--circular', '-c', action='store_true', help="Assume that all the input sequences are circular chromosomes or plasmids.")
    parser.add_argument('--seed', type=int, help="Seed value to initialize the random generator. Used for reproducibility.")

    args = parser.parse_args()

    if args.sub_sample <= 0:
        print("sub sample ratio must be positive and non-null", file=sys.stderr)
        exit(1)
    elif args.sub_sample > 1:
        print("sub sample ratio can't be over 1", file=sys.stderr)
        exit(1)

    return args


def parse_csv(filename):
    """ Parse the CSV file for a distribution and return a list of tuples (x, y).

    Args:
        filename: The CSV filename

    Return:
        A list of tuples (x, y) regarding the columns X and Y in the csv file.
    """

    tuples = []

    with open(filename, newline='') as csv_file:
        reader = csv.DictReader(csv_file, delimiter=',')

        for row in reader:
            t = (float(row["X"]), float(row["Y"]))
            if t[1] != 0:
                tuples.append(t)

    return tuples


def init_random_generators(len_distrib, seed=None):
    """ Init random generator regarding the distribution of read lengths

    Args:
        len_distrib: A CSV filename containing the distributions for the length.
        In this CSV file, the columns X and Y must be present.

    Return:
        A namespace containing three Metropolis hasting generators regarding the three different
        distributions.
    """
    if seed:
        random.seed(seed)

    # Generate the Metropolis hasting for insertions
    distrib = parse_csv(len_distrib)
    return MetropolisHasting1D(distrib, 10000, seed=random.randint(0, 1000000000))


def generate(filename, len_generator, depth, paired=False, circular=False, subsample_ratio=1, outprefix=None):
    """ Read the reference fasta file and write the out_file with reads. The lengths are
    randomized using the length generator.

    Args:
        filename: The input reference FASTA filename. The sequences will be splited into reads.
        If None, stdin is used instead.

        len_generator: A random generator for read length.

        depth: the output sequencing depth (ie, the average coverage on each nucleotide from the
        reference when the outputed reads are aligned)

        circular: Boolean to assume if the input sequences are 

        out_file: Filename for the outputed fasta. All the modified reads will be outputed in this
        file. The filename is optional. If the filename is not specified, then the outfile will
        be stdout

    Return: None
    """
    print("TODO: Fix problem with floating point depth", file=sys.stderr)

    full_depth = depth / subsample_ratio

    # Get the input name
    reader = None
    if not filename:
        reader = sys.stdin
    else:
        reader = open(filename)

    writer = None
    writer2 = None
    if not paired:
        # Get the outfile name
        if not outprefix:
            writer = sys.stdout
        else:
            writer = open(outprefix + ".fasta", "w")
    else:
        if not outprefix:
            writer = sys.stdout
            writer2 = sys.stderr
        else:
            writer = open(outprefix + '_R1.fasta', 'w')
            writer2 = open(outprefix + '_R2.fasta', 'w')

    # split sequences
    random_lengths = []
    read_idx = 0
    for record in SeqIO.parse(reader, "fasta"):
        str_seq = str(record.seq)

        # Depth
        for d in range(round(full_depth)):
            # Circular to linear part
            split_idx = 0
            if circular:
                # select the split idx
                split_idx = random.randint(0, len(str_seq))
                # reorganize the genome
                str_seq = "{}{}".format(str_seq[split_idx:], str_seq[:split_idx])

            # current idx on reference for next split
            left_idx = 0
            right_idx = len(str_seq)

            # Make reads from left and right
            while left_idx <= right_idx:
                # Length computation
                if len(random_lengths) == 0:
                    random_lengths = len_generator.next_values(10000)
                    random.shuffle(random_lengths)
                length, _ = random_lengths.pop(0)
                length = round(length)
                total_length = length
                if paired:
                    total_length = total_length * 2 #TODO : + gap

                # From left or right ?
                left = (random.random() <= 0.5)

                origin = split_idx
                if left:
                    sequence = str_seq[left_idx:min(right_idx+1, left_idx +total_length)]
                    origin += left_idx
                    left_idx += total_length
                else:
                    sequence = str_seq[max(left_idx, right_idx-total_length+1):right_idx+1]
                    origin += max(left_idx, right_idx-total_length+1)
                    right_idx -= total_length

                origin %= len(str_seq)

                # sub-sampling
                if random.random() > subsample_ratio:
                    continue

                # Reverse complement with 50% chance
                reversed = False
                if random.random() <= 0.5:
                    sequence = str(Seq(sequence).reverse_complement())

                if paired:
                    # Extract the beginning of the paired sequence
                    seq_1 = sequence[:min(length, len(sequence))]
                    # Extract the end of the paired sequence
                    seq_2 = sequence[max(0, len(sequence)-length):]
                    seq_2 = str(Seq(seq_2).reverse_complement())

                    # Output the reads
                    writer.write(">{}_{}_R1;origin={};\n{}\n".format(
                        record.description, read_idx, origin, seq_1))
                    # Output the reads
                    writer2.write(">{}_{}_R2;origin={};\n{}\n".format(
                        record.description, read_idx, origin, seq_2))
                else:
                    # Output the read
                    writer.write(">{}_{};origin={};\n{}\n".format(
                        record.description, read_idx, origin, sequence))

                read_idx += 1


    # Close the outfile if needed
    if not outprefix:
        writer.close()
        if paired:
            writer2.close()
    if not filename:
        reader.close()


def main():
    args = parse_arguments()
    rnd_gens = init_random_generators(
        args.length_distribution
    )
    #        filename,       len_generator, depth, subsample_ratio=1, outprefix=None
    generate(args.reference, rnd_gens, args.depth, circular=args.circular, paired=args.paired, subsample_ratio=args.sub_sample, outprefix=args.outprefix)


if __name__ == "__main__":
    main()

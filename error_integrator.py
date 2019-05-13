#!/usr/bin/env python3

import sys
import argparse
import csv
import random

from Bio import SeqIO
from metropolis_hasting import MetropolisHasting1D


def parse_arguments():
    # Parse the command line
    parser = argparse.ArgumentParser(description="Process FASTA files to integrate errors based on a PacBio model")
    parser.add_argument('--insertions_distribution', '-i', required=True, help="CSV file containing an insertion distribution. The columns X and Y must be present where X is the insertion rate and Y the number of reads with this rate.")
    parser.add_argument('--deletions_distribution', '-d', required=True, help="CSV file containing an deletion distribution. The columns X and Y must be present where X is the deletion rate and Y the number of reads with this rate.")
    parser.add_argument('--substitutions_distribution', '-s', required=True, help="CSV file containing an substitution distribution. The columns X and Y must be present where X is the substitution rate and Y the number of reads with this rate.")
    parser.add_argument('--outfile', '-o', help="The file where the modified reads will be saved. If not specified, outputed on stdout")
    parser.add_argument('--fasta', '-f', help="The fasta file containing all the reads to be modified. If not specified, use stdin as input")
    parser.add_argument('--seed', type=int, help="Seed value to initialize the random generator. Used for reproducibility.")

    args = parser.parse_args()

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


def init_random_generators(ins_distrib, del_distrib, sub_distrib, seed=None):
    """ Init random generator and add them into a Namespace

    Args:
        The three arguments must be CSV filenames containing the distributions for the errors.
        In these CSV files, the columns X and Y must be present.

    Return:
        A namespace containing three Metropolis hasting generators regarding the three different
        distributions.
    """
    if seed:
        random.seed(seed)

    # Generate the Metropolis hasting for insertions
    distrib = parse_csv(ins_distrib)
    ins_rnd = MetropolisHasting1D(distrib, 1000, seed=random.randint(0, 10000), burning_steps=10000)

    # Generate the Metropolis hasting for deletions
    distrib = parse_csv(del_distrib)
    del_rnd = MetropolisHasting1D(distrib, 1000, seed=random.randint(0, 10000), burning_steps=10000)

    # Generate the Metropolis hasting for substitutions
    distrib = parse_csv(sub_distrib)
    sub_rnd = MetropolisHasting1D(distrib, 1000, seed=random.randint(0, 10000), burning_steps=10000)

    return argparse.Namespace(ins=ins_rnd, dele=del_rnd, sub=sub_rnd)


def introduce_errors(filename, error_generators, outfile=None):
    """ Read all the sequence in the fasta file and write the out_file including errors in the
    reads. The error numbers are randomized using the error generators. The header of each output
    read willl contain a section errors between two ';'. In this section you'll find the number of
    error included for each type.

    Args:
        filename: The input FASTA filename. The sequences from this file will be modified to
        innclude errors. If none, stdin is used instead.

        error_generators: A namespace including 3 random generators (one for each error type).

        out_file: Filename for the outputed fasta. All the modified reads will be outputed in this
        file. The filename is optional. If the filename is not specified, then the outfile will
        be stdout

    Return: None
    """

    # Get the input name
    reader = None
    if not filename:
        reader = sys.stdin
    else:
        reader = open(filename)

    # Get the outfile name
    writer = None
    if not outfile:
        writer = sys.stdout
    else:
        writer = open(outfile, "w")

    # Include the errors
    for record in SeqIO.parse(reader, "fasta"):
        # Randomize the number of errors for the current read
        l = len(record.seq)
        nb_del = round(error_generators.dele.next_value()[0] * 0.01 * l)
        nb_sub = round(error_generators.sub.next_value()[0] * 0.01 * l)
        nb_ins = round(error_generators.ins.next_value()[0] * 0.01 * l)

        # Output the modified read
        writer.write(">{};errors:d{}_s{}_i{};\n{}\n".format(
            record.description, nb_del, nb_sub, nb_ins,
            modify_read2(str(record.seq), nb_del, nb_sub, nb_ins)))

    # Close the outfile if needed
    if not outfile:
        writer.close()
    if not filename:
        reader.close()


__sub_lists = {'A':'CGT', 'C':'AGT', 'G':'ACT', 'T':'ACG', 'a':'CGT', 'c':'AGT', 'g':'ACT', 't':'ACG'}

__ACGT = "ACGT"
__ACG = "ACG"
__ACT = "ACT"
__AGT = "AGT"
__CGT = "CGT"

def modify_read2(sequence, nb_del, nb_sub, nb_ins):
    """ Modify the input sequence adding deletions, substitutions and insertions. The errors are
    added simultaneously.

    Args:
        sequence: The DNA sequence to modify

        nb_del, nb_sub, nb_ins: The number of deteltions (resp. substitutions, insertions) to apply
        onto sequence in order to generate the errorneous read.

    Return: string
        Modified sequence with the inputed number of errors.
    """

    # Generate random operation sequence
    operations = ['d']*nb_del + ['s']*nb_sub + ['i']*nb_ins
    random.shuffle(operations)
    # Generate random idx selection
    idxs = list(random.sample(range(len(sequence)), nb_ins+nb_sub+nb_del))
    idxs.sort()

    # Init sequence and operations
    out_sequence = []

    prev_idx = 0
    for op, idx in zip(operations, idxs):
        # Copy all the sequence until the error
        out_sequence.append(sequence[prev_idx:idx])

        if op == 'd':
            prev_idx = idx+1
        elif op == 'i':
            out_sequence.append(__ACGT[random.randint(0,3)])
            prev_idx = idx
        else:
            n = sequence[idx]
            if n == 'A' or n == 'a':
                out_sequence.append(__CGT[random.randint(0,2)])
            elif n == 'C' or n == 'c':
                out_sequence.append(__AGT[random.randint(0,2)])
            elif n == 'G' or n == 'g':
                out_sequence.append(__ACT[random.randint(0,2)])
            else:
                out_sequence.append(__ACG[random.randint(0,2)])
            prev_idx = idx+1

    # Copy the end of the error
    out_sequence.append(sequence[prev_idx:])


    return "".join(out_sequence)


def modify_read(sequence, nb_del, nb_sub, nb_ins):
    """ Modify the input sequence adding deletions, substitutions and insertions. The errors are
    added deletions first, then substitutions and finally insertions. WARNING: deletion then
    insertion of the same letter at the same idx is not controlled (should be very rare).

    Args:
        sequence: The DNA sequence to modify

        nb_del, nb_sub, nb_ins: The number of deteltions (resp. substitutions, insertions) to apply
        onto sequence in order to generate the errorneous read.

    Return: string
        Modified sequence with the inputed number of errors.
    """

    # transform sequence to list
    seq_list = list(sequence)

    # Perform deletions
    del_idxs = random.sample(range(0, len(seq_list)), nb_del)
    del_idxs.sort(reverse=True)
    
    for idx in del_idxs:
        del(seq_list[idx])

    # Perform substitutions
    sub_idxs = random.sample(range(0, len(seq_list)), nb_sub)

    for idx in sub_idxs:
        seq_list[idx] = __sub_lists[seq_list[idx]][random.randint(0,2)]

    # Perform insertions
    ins_idxs = random.sample(range(0, len(seq_list)+1), nb_ins)
    ins_idxs.sort(reverse=True)

    for idx in ins_idxs:
        seq_list.insert(idx, "ACGT"[random.randint(0,3)])

    return "".join(seq_list)


def main():
    args = parse_arguments()
    rnd_gens = init_random_generators(
        args.insertions_distribution,
        args.deletions_distribution,
        args.substitutions_distribution,
        args.seed
    )
    introduce_errors(args.fasta, rnd_gens, args.outfile)


if __name__ == "__main__":
    main()

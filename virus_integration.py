#!/usr/bin/env python3

import sys
import argparse
import random
import time
from Bio import SeqIO


def parse_arguments():
    # Parse the command line
    parser = argparse.ArgumentParser(description="Include parts of the genome of a virus into an host genome. The sizes of inclusions are uniformaly distributed between min and max. The virus slices and the inclusion points in host are also uniformaly distributed into the genomes.")
    parser.add_argument('--min-size', '-min', type=int, default=0, help="Minimum size of the integrated parts [default=0].")
    parser.add_argument('--max-size', '-max', type=int, help="Maximum size of the integrated parts [default=len(virus)]")
    parser.add_argument('--random-seed', '-s', type=int, help="Random seed to use for reproducibility. If no specified, the random seed will be the current time.")

    parser.add_argument('--virus', '-vs', type=str, required=True, help="The virus reference (fasta file). Warning: The virus sequence is supposed to be formed by only one contig.")
    parser.add_argument('--host', '-ht', type=str, required=True, help="Host reference (fasta file). Can be a list of contigs if no full assembly exist.")
    parser.add_argument('--output', '-o', type=str, help="Fasta file that will contains the result. If not specified, redirect to stdout")
    parser.add_argument('--num-inclusions', '-n', type=int, required=True, help="Number of virus slices to integrate into the host")

    args = parser.parse_args()

    return args


def init(args):
    if args.random_seed:
        random.seed(args.random_seed)

    # Parse and return the virus sequence
    virus = SeqIO.read(args.virus, "fasta")
    if not args.max_size:
        args.max_size = len(virus.seq)

    return virus


def include_virus(host_filename, virus, slice_sizes, inclusion_positions, outfile=None):
    # open the out file
    fw = None
    if not outfile:
        fw = sys.stdout
    else:
        fw = open(outfile, "w")

    # Read the host and include
    for seq_idx, record in enumerate(SeqIO.parse(host_filename, "fasta")):
        local_inclusions = []

        # Get inclusions for this sequence
        while len(inclusion_positions) > 0 and inclusion_positions[0][0] == seq_idx:
            local_inclusions.append(inclusion_positions[0])
            inclusion_positions = inclusion_positions[1:]

        # Include the parts of the virus from the end of the sequence (allow to preserve inclusion idxs)
        sequence = record.seq
        for _, inclusion_idx in reversed(local_inclusions):
            # Get the virus slice
            slice = generate_virus_slice(virus, slice_sizes[0], slice_sizes[1])
            # Include the slice into the sequence
            sequence = "{}{}{}".format(sequence[:inclusion_idx], slice, sequence[inclusion_idx:])

        # Write the hybrid in the out fasta
        fw.write(">{};inclusions:{};\n{}\n".format(
            record.description.replace(" ", "_"),
            "_".join([str(j) for _,j in local_inclusions]),
            sequence))

    # Close the output file
    if not outfile:
        fw.close()


def generate_virus_slice(virus, min_size, max_size):
    # Size randomization
    size = random.randint(min_size, max_size)
    max_start_idx = len(virus.seq) - size

    # First nucleotide
    first_idx = random.randint(0, max_start_idx)
    # Last nucleotide
    last_idx = first_idx + size

    return virus.seq[first_idx:last_idx]


def generate_inclusion_positions(host_filename, nb_inclusions):
    host_sizes = []
    total_inclusion_points = 0

    # Read all the sequences composing the host to get all the possible inclusion points
    for record in SeqIO.parse(host_filename, "fasta"):
        host_sizes.append(len(record.seq))
        # +1 because there insertion points before AND after the sequence
        total_inclusion_points += host_sizes[-1] + 1

    # Randomize selection of insertion points
    selected_insertion_points = []
    for _ in range(nb_inclusions):
        selected_insertion_points.append(random.randint(0, total_inclusion_points))
    selected_insertion_points.sort()

    # Transform insertion points into tuples (seq_id, insertion_idx)
    selected_insertion_tuples = []
    cumulative_insert_idx = 0
    current_host_sequence = 0
    for idx in selected_insertion_points:
        # Continue to the next host sequence that contain the insertion point
        while idx > cumulative_insert_idx + (host_sizes[current_host_sequence] + 1):
            cumulative_insert_idx += host_sizes[current_host_sequence]
            current_host_sequence += 1

        # Add the tuple
        selected_insertion_tuples.append((current_host_sequence, idx - cumulative_insert_idx))

    return selected_insertion_tuples


def main():
    # Parse arguments
    args = parse_arguments()
    # Init the random generator and parse the virus sequence once for all
    virus = init(args)
    # Get inclusion positions on the host
    positions = generate_inclusion_positions(args.host, args.num_inclusions)
    # Outpur the sequences with inclusions
    include_virus(args.host, virus, (args.min_size, args.max_size), positions, args.output)


if __name__ == "__main__":
    main()

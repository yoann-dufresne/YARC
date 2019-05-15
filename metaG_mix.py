#!/usr/bin/env python3

import sys
import argparse
import os
import random
import shutil

import reads_generator


def parse_arguments():
    # Parse the command line
    parser = argparse.ArgumentParser(description="Use a bunch of fasta file to create a metagenomic synthetic sample.")
    parser.add_argument('--genome_directory', '-d', required=True, help="A directory containing fasta files that will be used for the sample creation.")
    parser.add_argument('--num_genomes', '-n', type=int, required=True, help="Number of genomes that must be includded")
    parser.add_argument('--coverage', '-c', type=int, required=True, help="Average coverage.")
    parser.add_argument('--outdir', '-o', default='.', help="Output directory")
    parser.add_argument('--seed', '-s', type=int, help="The random seed (use it only for reproducibility purpose).")
    parser.add_argument('--length_distribution', '-l', help="Read length distribution")
    parser.add_argument('--only_distribution', '-od', action='store_true', help="Stop the execution after creating a distribution file.")

    args = parser.parse_args()

    # Verification of directory correctness
    path = args.genome_directory
    if not os.path.isdir(path):
        print(f"{path} is not a directory", file=sys.stderr)
        exit(1)
    if path[-1] == "/":
        args.genome_directory = path[:-1]

    path = args.outdir
    if not os.path.isdir(path):
        os.mkdir(path)
    if path[-1] == "/":
        args.outdir = path[:-1]

    # Conflicting parameters
    if (not args.only_distribution) and (not args.length_distribution):
        print("A csv length distribution must be provided", file=sys.stderr)
        exit(1)

    # Fix random seed
    if args.seed:
        random.seed(args.seed)

    print("TODO: add threading support", file=sys.stderr)

    return args


def select_genomes(directory, n):
    """ Select n genomes from the directory.
        @args directory A directory containing fasta, fa or fna files.
        @args n The number of genome to select
        @except Raise an exception when not enought files (regarding n) are in directory.
        @return A list of n filenames.
    """
    # Select files from the directory, filetering out the wrong extentions
    supported_ententions = ["fa", "fasta", "fna"]
    files = [f for f in os.listdir(directory) if f.split('.')[-1] in supported_ententions]

    # Sub-select n files
    if len(files) < n:
        raise Exception(f"Impossible to find {n} distinct genomes in {directory}.")
        return []
    return random.sample(files, n)


def create_mix(files):
    """ Compute a relative abundance of each genome using a log-norm distribution.
        The distribution is centered on 1 with a standard deviation of 2 (from CAMISIM).
        A relative abundance of 1 significate that the genome associated must have
        a coverage equivalent to the average coverage.
        @args files All the files to compute abundances.
        @return A dictionnary associating each filename to a relative abundance.
    """
    relative_dist = {}
    for f in files:
        relative_dist[f] = random.lognormvariate(1, 1.5)

    # Normalise to center on 1
    vals = relative_dist.values()
    avg = sum(vals) / len(vals)
    relative_dist = {k:v/avg for k, v in relative_dist.items()}

    return relative_dist

# def generate(filename, len_generator, depth, paired=False, circular=False, subsample_ratio=1, outprefix=None)

from datetime import datetime
def create_out_tree(outdir):
    # Create output directory tree
    sample_dir = outdir + '/' + "_".join(str(datetime.now()).split())
    os.mkdir(sample_dir)
    shutil.move(outdir + '/abundances.csv', sample_dir + '/abundances.csv')
    out_gen_dir = sample_dir + '/genomes'
    os.mkdir(out_gen_dir)
    reads_dir = sample_dir + '/reads'
    os.mkdir(reads_dir)

    return sample_dir


def create_sequences(abundances, gen_dir, len_dist, sample_dir):
    """ Create one read file for each selected genome regarding the abundance in the sample.
        The used genomes will be copied into a genome folder into the sample directory.
        The reads will be copied into a reads folder.
        @args abundances A dictionary associating genomes to abundances
        @args gen_dir The directory containing the genome fasta files
        @args len_dist A read length distribution must be provided. See reads_generator for more details.
        @args sample_dir The directory where everything will be outputed.
    """
    # Sampling
    for genome, abundance in abundances.items():
        # Copy the complete genome
        shutil.copy(f"{gen_dir}/{genome}", f"{sample_dir}/genomes/{genome}")

        # Generate reads
        rnd_gens = reads_generator.init_random_generators(len_dist)
        reads_generator.generate(
            f"{sample_dir}/genomes/{genome}",
            rnd_gens,
            abundance,
            paired=False,
            subsample_ratio=0.5,
            outprefix=f"{sample_dir}/reads/{'.'.join(genome.split('.')[:-1])}"
        )


def pull_files(files, sample_dir):
    """ Take all the content of each file from files and pull them together.
        TODO: The reads must be shuffled.
        @args files the list of file to concatenate
        @args sample_dir the directory of the sample to pull
    """
    with open(f"{sample_dir}/pulled_reads.fa", "w") as pulled:
        for f in files:
            core_name = ".".join(f.split(".")[:-1])
            with open(f"{sample_dir}/reads/{core_name}.fasta") as fr:
                for line in fr:
                    pulled.write(line)


def main():
    args = parse_arguments()

    print("1 - Selecting genomes", file=sys.stderr)
    files = select_genomes(args.genome_directory, args.num_genomes)

    print("2 - Generate abundances", file=sys.stderr)
    # Compute relative distribution
    relative_abundances = create_mix(files)
    abundances = {k:v*args.coverage for k,v in relative_abundances.items()}

    # Save distribution
    with open(f"{args.outdir}/abundances.csv", "w") as ab_fw:
        ab_fw.write("genome,abundance\n")
        for gen, ab in abundances.items():
            ab_fw.write(f"{gen},{ab}\n")

    if args.only_distribution:
        return

    print("3 - Create reads", file=sys.stderr)
    sample_dir = create_out_tree(args.outdir)
    create_sequences(abundances, args.genome_directory, args.length_distribution, sample_dir)

    print("4 - Pull reads together", file=sys.stderr)
    pull_files(files, sample_dir)
    print("TODO: Shuffle the concatenate reads", file=sys.stderr)


if __name__ == "__main__":
    main()

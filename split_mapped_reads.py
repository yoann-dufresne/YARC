import sys, argparse
import paf_parser as paf
from Bio import SeqIO


def parse_arguments():
    # Parse the command line
    parser = argparse.ArgumentParser(description="Use a mapping file to extract mapped and flanking regions from the FASTA mapped file. (file mapped on the reference)")
    parser.add_argument('--alignment', '-a', type=str, required=True, help="Alignment file (PAF format)")
    parser.add_argument('--fasta', '-f', type=str, required=True, help="FASTA file with aligned reads")
    parser.add_argument('--reference', '-r', type=str, required=True, help="FASTA file with the reference where the reads are aligned")
    parser.add_argument('--size', '-s', type=int, help="Maximum size to extract (left and right)")
    parser.add_argument('--prefix', '-p', type=str, default="", help="Output FASTA file with flanking regions. Left and right are pasted together.")

    args = parser.parse_args()

    return args


def main(args):
    alignments = read_alignments(args.alignment)
    prefix = args.prefix
    if prefix == "":
        prefix = args.fasta[:args.fasta.rfind('.')]
    extract_flanking(args.fasta, args.reference, alignments, prefix)


def read_alignments(paf_filename):
    alignments = {}

    # Read alignment
    with open(paf_filename) as paf_reader:
        for idx, alignment in enumerate(paf.parse(paf_reader)):
            print(idx)
            if not alignment["query"]["name"] in alignments:
                alignments[alignment["query"]["name"]] = []
            alignments[alignment["query"]["name"]].append(alignment)
            
    return alignments


def extract_flanking(fasta, virus_infile, alignments, prefix):
    # Open the 3 files to outputs the parts of the sequences splitted
    with open(f"{prefix}_left.fasta", "w") as left, open(f"{prefix}_right.fasta", "w") as right, open(f"{prefix}_virus_ref.fasta", "w") as virus_ref, open(f"{prefix}_virus_reads.fasta", "w") as virus_reads:

        # load virus
        virus_sequences = {}
        for record in SeqIO.parse(virus_infile, "fasta"):
            virus_sequences[record.name] = str(record.seq)

        # For all sequences in the reads
        align_idx = 0
        for record in SeqIO.parse(fasta, "fasta"):
            # If the read is aligne
            if record.description in alignments:
                align_idx += 1
                # For each match on this read
                for sub_idx, alignment in enumerate(alignments[record.description]):
                    # Spit this read into the 3 files.
                    print(f"Process alignment {align_idx}[{sub_idx+1}]/{len(alignments)}:\n{alignment}", "\n")
                    # Left read extraction
                    seq = str(record.seq)
                    seq = seq[:alignment["query"]["start"]]
                    left.write(f">{alignment['query']['name']}\n{seq}\n")

                    # Right read extraction
                    seq = str(record.seq)
                    seq = seq[alignment["query"]["end"]+1:]
                    right.write(f">{alignment['query']['name']}\n{seq}\n")

                    # virus part from read
                    seq = str(record.seq)
                    seq = seq[alignment["query"]["start"]: alignment["query"]["end"]+1]
                    virus_reads.write(f">{alignment['query']['name']} read sequence\n{seq}\n")
                    
                    # virus part from reference
                    seq = virus_sequences[alignment["target"]["name"]]
                    seq = seq[alignment["target"]["start"]: alignment["target"]["end"]+1]
                    virus_ref.write(f">{alignment['target']['name']} virus sequence\n{seq}\n")
                

if __name__ == "__main__":
    args = parse_arguments()
    main(args)


# Yet Another Random Code
Repository of unrelated useful bioinformatics scripts

* `reads_generator.py` 

is an ad-hoc NGS reads simulator (sort of like wgsim), based on Metropolis-Hastings sampling for read lengths, which generates simulated reads *without any sequencing error*

* `introduce_errors`

takes a FASTA file as input, typically reads produced by `reads_generator.py`, and produces an output FASTA file with sequencing errors introduced randomly in the input, according to a customized error profile. Typically used to simulate PacBio reads.

* `metaG_mix.py` 

simulates metagenomics NGS reads (sort of like CAMISIM) with a lognorm abundance distribution over the genomes

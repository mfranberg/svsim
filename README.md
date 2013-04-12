# svsim

The svsim package provides a collection of tools for simulating structural variations. Specifically it automates the process of generating structural variations in given reference contigs, simulating reads from the donor contigs and mapping those reads back to the reference genome.

The toolbox consists of the following scripts:
* *simulate_genome.py*, generates a synthetic genome from given base pair frequencies.
* *create_donor_contigs.py*, takes a a set of contigs and inserts structural variations in them.
* *simulate_reads.py*, simulates reads with metasim or dwgas from given contigs.
* *map_reads.py*, maps reads with bwa to a given reference genome.
* *svsim_pipeline.py*, combines the previous 3 scripts into one.

## simulate_genome.py

Generates a genome of the given length and base pair frequencies. The generated genome will be outputted in fasta format. It is used as follows:

    > python simulate_genome.py 10000 genome.fa

If you want to change the base pair frequencies to be uniform, you can change them with the -p flag, and supply each base pair frequency for A, C, G and T separated by space:

    > python simulate_genome.py -p 0.2 0.3 0.3 0.2 10000 genome.fa

## create_donor_contigs.py

This script generates contigs with the given structural variations from a set of reference contigs. The structural variations are specified in a custom file format that support 4 types of structural variations: insertion, deletion, duplication and translocation.

Structural variations are defined by lines:

    contig_name insertion start length
    contig_name deletion start length
    contig_name duplication start length to
    contig_name translocation start length to
    
where *contig_name* is the name of a contig in the contig fasta file, *start* is the start position of the variation in the contig, *length* is the length of the variation, and to is the location where the copied or removed segment will end up. It is not allowed to have overlapping structural variations. Start positions are always defined as one position before the actual variation.

I will illustrate this with an example. Let genome.fa be:

    >contig1
    AAACCCGGGTTT
    >contig2
    AACCGGTT
    
and variations.txt be:

    contig1 deletion 3 3
    contig1 duplication 6 3 9
    
If we now run:

    > python create_donor_contigs.py genome.fa variations.txt indel_genome.fa
    > cat indel_genome.fa
    >contig1-donor
    AAAGGGGGGTTT
    >contig2-donor
    AACCGGTT
    
As we can see the three C's has been deleted and the three G's has been duplicated in contig1, the other contig is kept intact.

## simulate_reads.py

This program uses [MetaSim](http://ab.inf.uni-tuebingen.de/software/metasim/) or [dwgsim](https://github.com/nh13/DWGSIM) to simulate reads from a given set of contigs. To use this script you need to have MetaSim or dwgsim in your PATH.

    usage: simulate_reads.py [-h] [-c C] [-m M] [-s S] -t {metasim,dwgsim}
                             genome_file output_prefix
    
    Simulates illumina reads with metasim.
    
    usage: simulate_reads.py [-h] [-c C] [-m M] [-s S] -t {metasim,dwgsim}
                             genome_file output_prefix
    
    Simulates illumina reads with metasim.
    
    positional arguments:
      genome_file          Path to the genome
      output_prefix        Output prefix for the paired end files, will apped
                           _pe1.fa and pe2.fa to this.

    optional arguments:
      -h, --help           show this help message and exit
      -c C                 Coverage.
      -m M                 Mean of the library distribution.
      -s S                 Standard deviation of the library distribution.
      -t {metasim,dwgsim}  Type of simulator 'metasim' or 'dwgsim'.

To simulate reads with metasim with default parameters do:

    > python simulate_reads.py -t metasim genome.fa reads
    
This will produce reads_pe1.fa and reads_pe2.fa which are two fasta files.

To simulate reads with 30 coverage, 400 mean insert size and 50 in standard deviation run:

    > python simulate_reads.py -t metasim -c 30 -m 400 -s 50 genome.fa reads
    
This will again produce reads_pe1.fa and reads_pe2.fa.

Note: MetaSim can be complicated to get working in your PATH. One way is to create a bash-script that forwards the command to the real MetaSim command. The script is called MetaSim and looks like:

    #!/bin/sh
    
    /Applications/metasim/MetaSim $*
    exit $?


## map_reads.py

This script maps paired reads with bwa. You need to have bwa in your path in order to run this script. To map the reads *reads_pe1.fa* and *reads_pe2.fa* to the reference genome *genome.fa* do:

    > python map_reads.py reads_pe1.fa reads_pe2.fa genome.fa mapped_reads
    
This will produce the file mapped_reads.bam which is a sorted .bam file.

## svsim_pipeline.py

All of the above scripts are combined in a single script called svsim_pipeline.py. To generate mapped reads for normal (*genome.fa*) and mutated contigs, where the mutated contigs has structural variations defined in *variations.txt*, with 50 coverage and the default library distribution using *dwgsim* do:

    > python svsim_pipeline.py -t dwgsim -c 50 genome.fa variations.txt output_dir/

The ouput directory now contains the mutated contigs, the unmapped read pairs and the mapped read pairs:

    > ls output_dir/
    donor_contigs.fa
    mapped_donor.bam
    mapped_normal.bam
    reads_donor_pe1.fa
    reads_donor_pe2.fa
    reads_normal_pe1.fa
    reads_normal_pe2.fa

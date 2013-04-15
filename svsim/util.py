import pyfasta

##
# Calculate the number of reads required to get a certain coverage.
#
# @param coverage Desired mean genome coverage.
# @param read_length Length of each read in a pair.
# @param genome_length Estimated genome length.
#
# @return Number of reads need to get the given coverage.
#
def calculate_num_reads(coverage, read_length, genome_length):
    return int( ( genome_length * coverage ) / ( read_length * 2.0 ) )

##
# Opens the given fasta file and calculates the total length
# of all contigs.
#
# @param genome_path Path to a fasta file.
#
# @return The length of the genome.
#
def get_genome_length(genome_path):
    genome_fasta = pyfasta.Fasta( genome_path )
    return sum( len( genome_fasta[ contig ] ) for contig in genome_fasta.iterkeys( ) )

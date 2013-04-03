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
# Opens the given fasta file and calculates the length of the
# genome.
#
# @param genome_path Path to a fasta file containing a single genome.
#
# @return The length of the genome.
#
def get_genome_length(genome_path):
    with open( genome_path, "r" ) as genome_file:
        next( genome_file ) # Skip header
        return sum( len( line.strip( ) ) for line in genome_file )

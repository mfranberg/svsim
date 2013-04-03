import argparse
import sys
import random

##
# Generates base pairs according to the given
# distribution.
#
# @param probs A list of probabilities for each base pair.
# @param bp A list of base pairs.
#
def generate_bp(probs, bp):
    rand = random.random( )
    total = 0.0
    for i in range( len( probs ) ):
        total += probs[ i ]
        if rand <= total:
            return bp[ i ]

##
# Generates the genome a returns it.
#
# @param probs Base pair probabilities.
# @param length Length of the genome.
#
# @return The sequence of the generated genome.
#
def generate_genome(probs, length):
    bp = "ACGT"
    return ''.join( generate_bp( probs, bp ) for i in range( length ) )

##
# Writes the genome to the given output file in
# fasta format.
#
# @param genome Sequence of the genome to write.
# @param output_file Output file to write the genome to.
#
def write_genome(genome, output_file):
    output_file.write( ">normal-genome\n" )
    output_file.write( genome )
    output_file.write( '\n' )

USAGE = """Usage: simulate_genome length output_file"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser( description="Simulates a genome from the given base probabilities." )
    parser.add_argument( 'length', type=int, help='Path to the normal genome file.' )
    parser.add_argument( '-p', type=float, nargs=4, help="List of frequencies for A, C, G and T respectively.", default=[0.25,0.25,0.25,0.25] )
    parser.add_argument( 'output_file', type=argparse.FileType( 'w' ), help='Output path of the mutated genome.' )
    args = parser.parse_args( )
    
    assert( abs( sum( args.p ) - 1.0 ) <= 0.001 )
    genome = generate_genome( args.p, args.length )
    write_genome( genome, args.output_file )

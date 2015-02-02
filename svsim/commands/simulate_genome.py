
from __future__ import unicode_literals, print_function

import sys
import random

import click

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

def genome_generator(probabilities, length):
    """
    Generate bases with the given probabilities
    
    Args:
        probabilities (list): A list with the probabilities for each nucleotide
        length (int): The length of the genome
    
    Returns:
        A generator that produces letters of [A,C,G,T]
    """
    bp = "ACGT"
    for i in range( length ):
        yield generate_bp( probabilities, bp )

##
# Writes the genome to the given output file in
# fasta format.
#
# @param genome Sequence of the genome to write.
# @param output_file Output file to write the genome to.
#

def write_genome(genome_generator, output_file):
    output_file.write( ">normal-genome\n" )
    for base in genome_generator:
        output_file.write( base )
    output_file.write( '\n' )


@click.command()
@click.argument('length',
                    type=int
)
@click.argument('output',
                    type=click.File('wb'),
)
@click.option('-f', '--frequencies',
                    nargs=4,
                    default=[0.25,0.25,0.25,0.25],
                    help="List of frequencies for A, C, G and T respectively."
)
def simulate_genome(length, output, frequencies):
    """
    Simulates a genome from the given base probabilities..
    """
    frequencies = [float(frequency) for frequency in frequencies]
    try:
        assert( abs( sum( frequencies ) - 1.0 ) <= 0.001 )
    except AssertionError:
        print('The sum of the frequencies must equal 1', file=sys.stderr)
        sys.exit()
        
    write_genome( genome_generator( frequencies, length ), output )
    

if __name__ == '__main__':
    simulate_genome()

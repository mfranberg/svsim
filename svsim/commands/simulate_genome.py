
from __future__ import unicode_literals, print_function

import sys
import random

import click

def generate_bp(probs, bp):
    """
    Generates base pairs according to the given distribution.
    
    Arguments:
        probs (list): Probabilities for each base pair.
        bp (list): A list of base pairs.
    
    Returns:
        base_pair (str): The generated basepair
    """
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

def write_genome(output_file, contigs, genome_name, probabilities):
    """
    Writes the genome to the given output file in fasta format.
    
    Arguments:
        output_file (file_handle): Output file to write the genome to.
        contigs (dict): A dictionary with {<contig_name>: <contig_length>}
        genome_name (str): A string that represents the genome name
        probabilities (list): A list with the probabilities for each position. 
    """
    for contig in contigs:
        contig_length = contigs[contig]
        output_file.write(">{0}|dna:chromosome|chromosome:{1}:{0}:1:{2}:1|REF\n".format(
                                contig,
                                genome_name,
                                contig_length
                                )
                            )
        for base in genome_generator(probabilities, contig_length):
            output_file.write( base )
        output_file.write( '\n' )
    return


@click.command()
@click.argument('output',
                    type=click.File('wb'),
)
@click.option('-p', '--probabilities',
                    nargs=4,
                    type=float,
                    default=[0.25,0.25,0.25,0.25],
                    help="List of probabilities for A, C, G and T respectively."
)
@click.option('-f', '--contig_file',
                    type=click.File('r'),
                    help="Specify contigs and their length in a tab separated file."
)
@click.option('-n', '--genome_name',
                    default='normal-genome',
                    help="Add a contig that should be included."
)
@click.option('-c', '--contig',
                    multiple=True,
                    default=['1'],
                    help="Add a contig that should be included."
)
@click.option('-l','--contig_length',
                    multiple=True,
                    default=[1000],
                    help=("Specify the length of the contig"
                        "(in same order as the contigs are given).")
                    
)
def simulate_genome(output, probabilities, genome_name, contig, contig_length, contig_file):
    """
    Simulates a fasta genome from the given base probabilities.
    
    User can specify multiple contigs by using -c/--contig.
    Each contig needs to have a length specified (use -l/--contig_length).
    genome_name is the name of the genome (default normal-genome).
    
    Each contig will get a header line in the fasta file that looks like:
    ><contig>|dna:chromosome|chromosome:<genome_name>:<contig>:<start>:<stop>:1|REF
    """
    try:
        assert( abs( sum( probabilities ) - 1.0 ) <= 0.001 )
    except AssertionError:
        print('The sum of the probabilities must equal 1', file=sys.stderr)
        sys.exit()
    try:
        assert( len( contig ) == len( contig_length) )
    except AssertionError:
        print('Each contig must have a length', file=sys.stderr)
        sys.exit()
    contigs = {}
    if contig_file:
        try:
            for line in contig_file:
                line = line.rstrip().split()
                if len(line) > 1:
                    contigs[line[0]] = int(line[1])
        except ValueError:
            print("Contig length must be a integer.", file=sys.err)
            sys.exit()
    for i, contig in enumerate(contig):
        contigs[contig] = int(contig_length[i])
    
    write_genome( output, contigs, genome_name, probabilities)
    

if __name__ == '__main__':
    simulate_genome()


from __future__ import unicode_literals, print_function

import sys
import os
import click
import random

from StringIO import StringIO
from collections import defaultdict
from pprint import pprint as pp

import pyfasta

from svsim import variation, vcf

##
# Writes a donor contig to the given output file.
#
# @param donor_contig The contig to write.
# @param contig_name The id of the contig, prefix will be appeneded.
# @param output_file The file to write to.
# @param postfix Will be appended to the contig_name.
#
def write_donor_contig(donor_contig, contig_name, output_file, postfix="-donor"):
    output_file.write( ">{0}{1}\n".format( contig_name, postfix ) )
    
    for i in range( 0, len( donor_contig ), 80 ):
        output_file.write( donor_contig[ i:i+80 ] )
        output_file.write( '\n' )

##
# Create a variation object from the given variation_type
# and parameters.
#
# @param variation_type The type of variation, e.g. "insertions".
# @param params Parameters of the particular variation.
#
# @return A list of variations created from the given variation_type,
#         translocations for example create both deletions and insertions.
#
def create_sv(variation_type, params):
    if variation_type == "insertion":
        return [ variation.Insertion( params[ 0 ], params[ 1 ], -1 ) ]
    elif variation_type == "deletion":
        return [ variation.Deletion( params[ 0 ], params[ 1 ] ) ]
    elif variation_type == "duplication":
        return [ variation.Insertion( params[ 2 ], params[ 1 ], params[ 0 ] ) ]
    elif variation_type == "translocation":
        return [ variation.Deletion( params[ 0 ], params[ 1 ] ),
                 variation.Insertion( params[ 2 ], params[ 1 ], params[ 0 ] ) ]
    else:
        return None

##
# 
# 
#
# @param 
#
# @return 
#
def read_variations(variation_file):
    """
    Reads the variations defined in a file and returns
    them as a list of variation objects.
    
    Args:
        variation_file (file): A file that defines variations.
    
    Returns:
        variations (dict): A dict containing the variations for each contig.
    """
    variations = defaultdict( list )
    for line_number, line in enumerate( variation_file ):
        column = line.split( )
        
        contig_name = column[ 0 ]
        variation_type = column[ 1 ]
        params = [ int( v ) for v in column[ 2: ] ]
        
        sv = create_sv( variation_type, params )
        if sv:
            variations[ contig_name ].extend( sv )
        else:
            print( "warning: Ignored bad variation on line: {0}".format( line_number ) )

    return variations



##
# Creates donor contigs that contains the given
# variations, and writes them to an output file.
#
# Note: If a contig does not have any variations it will
#       still be written to the output file.
#
@click.command()
@click.argument('normal_contig_file',
                    type=click.Path(exists=True),
)
@click.argument('variation_file',
                    type=click.File('r'),
)
@click.argument('output',
                    type=click.File('w'),
)
@click.option('-d', '--delimiter',
                    default='|',
                    help='The fasta identifier separator, default is |.'
)
@click.option('-i', '--field_index',
                    nargs=1,
                    default=0,
                    help='The 0-based index of field separated by field_sep that contains the relevant contig name.'
)
@click.option('-c', '--chrom',
                    nargs=1,
                    type=str,
                    help='Sets chromosome that will be written in the vcf file, name of the contig by default.'
)
@click.option('-v', '--vcf_file',
                    type=click.Path(exists=False),
                    help='Output file of the structural variations in VCF format.'
)
def create_donor_contigs(normal_contig_file, variation_file, output, delimiter,
    field_index, chrom, vcf_file):
    """
    Creates the donor contigs with structural variations.
    
    The format of variation_file are lines consisting of the following: 
    contig_name insertion start length, contig_name deletion start length, 
    contig_name duplication start length to or contig_name translocation start length to.
    """
    
    if vcf_file and chrom:
        vcf_file = vcf.open_vcf_file
        vcf_file.set_chrom( args.chrom )
    
    normal_contigs = pyfasta.Fasta(
                                normal_contig_file,
                                key_fn = lambda x: x.split( delimiter )[ field_index ].strip( )
                                )
    
    variations = read_variations( variation_file )
    
    for contig_name, sequence in normal_contigs.iteritems( ):
        contig_variations = variations.get( contig_name, None )
        if contig_variations:
            donor_contig = variation.create_indel_genome( sequence, contig_variations )
            write_donor_contig( donor_contig, contig_name, output )

            if vcf_file:
                vcf_file.write( contig_name, sequence, contig_variations )
        else:
            write_donor_contig( sequence, contig_name, output )


if __name__ == '__main__':
    create_donor_contigs()

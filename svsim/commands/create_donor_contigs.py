import argparse
import sys
import random
import os.path
from StringIO import StringIO
from collections import defaultdict

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
# Reads the variations defined in a file and returns
# a them as a list of variation objects.
#
# @param variation_file A file that defines variations.
#
# @return A dict containing the variations for each contig.
#
def read_variations(variation_file):
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
# @param normal_contig_path Path of the fasta file that contains normal contigs.
# @param variation_file File that contains variations.
# @param donor_contig_file Output file where the donor contigs will be written.
#
# @param field_index 0-based index of contig_name in the fasta identifier.
# @param field_sep Field separator in the fasta identifier.
# @param vcf_file If not None vcf data will be written to this file object.
#
def create_donor_contigs(normal_contig_path, variation_file, donor_contig_file, field_index = 0, field_sep = "|", vcf_file = None):
    normal_contigs = pyfasta.Fasta( normal_contig_path, key_fn = lambda x: x.split( field_sep )[ field_index ].strip( ) )
    variations = read_variations( variation_file )

    for contig_name, sequence in normal_contigs.iteritems( ):
        contig_variations = variations.get( contig_name, None )
        if contig_variations:
            donor_contig = variation.create_indel_genome( sequence, contig_variations )
            write_donor_contig( donor_contig, contig_name, donor_contig_file )

            if vcf_file:
                vcf_file.write( contig_name, sequence, contig_variations )
        else:
            write_donor_contig( sequence, contig_name, donor_contig_file )

##
# Adds the arguments relevant for generation of the
# donor contigs.
#
# @param group An argument group to add options to.
#
def add_sv_arguments(group):
    group.add_argument( '-d', dest="field_sep", type=str, help='The fasta identifier separator, default is |.', default="|" ) 
    group.add_argument( '-i', dest="field_index", type=int, help='The 0-based index of field separated by field_sep that contains the relevant contig name.', default=0 ) 
    group.add_argument( '-v', dest="vcf_file", type=vcf.open_vcf_file, help='Output file of the structural variations in VCF format.' )
    group.add_argument( '-c', dest="chrom", type=str, help='Sets chromosome that will be written in the vcf file, name of the contig by default.' )


USAGE = """Usage: create_donor_contigs normal_contig_file variation_file donor_contig_file"""

VARIATION_USAGE = """Path to the variation file. The format of variation_file are lines consisting of the following: 
contig_name insertion start length, contig_name deletion start length, contig_name duplication start length to or contig_name translocation start length to."""

if __name__ == '__main__':
    parser = argparse.ArgumentParser( description=USAGE )
    parser.add_argument( 'normal_contig_file', type=str, help='Path to the normal contigs.' )
    parser.add_argument( 'variation_file', type=argparse.FileType( "r" ), help=VARIATION_USAGE )
    parser.add_argument( 'donor_contig_file', type=argparse.FileType( "w" ), help='Output file, the donor contigs will be written here.' )

    add_sv_arguments( parser )

    args = parser.parse_args( )

    if args.vcf_file and args.chrom:
        args.vcf_file.set_chrom( args.chrom )

    create_donor_contigs( args.normal_contig_file,
                          args.variation_file,
                          args.donor_contig_file,
                          field_index = args.field_index,
                          field_sep = args.field_sep,
                          vcf_file = args.vcf_file )

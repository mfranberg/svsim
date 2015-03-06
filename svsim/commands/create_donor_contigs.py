
from __future__ import unicode_literals, print_function

import sys
import os
import click
import random

from StringIO import StringIO
from collections import defaultdict
from pprint import pprint as pp

from pyfasta import Fasta

from svsim import variation, vcf

VARIATIONS = ("insertion", "deletion", "duplication", "translocation", 
                "transversion")

def write_donor_contig(donor_contig, contig_name, output_file, postfix="donor"):
    """
    Writes a donor contig to the given output file.
    
    Arguments:
        donor_contig The contig to write.
        contig_name The id of the contig, prefix will be appeneded.
        output_file The file to write to.
        postfix Will be appended to the contig_name.
    """
    output_file.write( ">{0}{1}\n".format( contig_name, postfix ) )
    
    for i in range( 0, len( donor_contig ), 80 ):
        output_file.write( donor_contig[ i:i+80 ] )
        output_file.write( '\n' )

def create_sv(variation_type, contig, pos, length, from_loc = -1, 
                from_contig=None, nr_duplications = 1):
    """
     Create a variation object from the given variation_type and parameters.
     
     Arguments:
         variation_type (str): The type of variation, e.g. "insertion".
         contig (str): The contig where the sv is located
         pos (int): The 0-based position of the base pair before
                   the sv.
         length (int): The length of the event.
         from_loc (int): The 0-based position from where the event sequence is
                         taken. Can be -1 then a random sequence is generated.
         from_contig (str): The contig where the sv originates. 
                            If None => same contig as sv is at.
         nr_duplications (int): If duplication this is the number of times that
                                 the sequence is duplicated.
         
    Returns:
        variations (list): A list ov variation objects
    
    """
    if variation_type == "insertion":
        return variation.Insertion( contig, pos, length, from_loc, from_contig)
    if variation_type == "transversion":
        return variation.Transversion( contig, pos, length)
    elif variation_type == "deletion":
        return variation.Deletion( contig, pos, length )
    elif variation_type == "duplication":
        return variation.Duplication( contig, pos, length, nr_duplications )
    elif variation_type == "translocation":
        return variation.Translocation( contig, pos, length, from_loc, from_contig )
    else:
        return None

def read_variations(variation_file, contigs):
    """
    Reads the variations defined in a file and returns
    them as a list of variation objects.
    
    Args:
        variation_file (file): A file that defines variations.
        contigs (list): A list with the existing contigs from the
                        reference genome
    
    Returns:
        variations (dict): A dict containing the variations for each contig.
    """
    variations = defaultdict( list )
    variations = []
    
    for line_number, line in enumerate( variation_file ):
        use_sv = True
        column = line.split( )
        contig_name = column[ 0 ]
        contig_from = contig_name
        from_loc = -1
        # Test if the contig exists
        if contig_name not in contigs:
            ##TODO Raise proper Error
            print( "Error: Contig does not exist in reference: {0}".format(
                                                                    line_number 
                                                                    ) 
                                                                )
            print("Skipping sv")
            use_sv = False
        
        variation_type = column[ 1 ]
        
        if variation_type not in VARIATIONS:
            print("Error: unknown variation: {0} on line number: {1}".format(
                                                            variation_type,
                                                            line_number
                                                            )
                                                           )
            print("Skipping sv")
            use_sv = False
        
        position = int(column[ 2 ])
        length = int(column[3])
        
        nr_duplications = 1
        if variation_type == 'duplication':
            if len(column) > 4:
                nr_duplications = int(column[4])
        
        elif variation_type == 'insertion':
            if len(column) > 5:
                contig_from = column[4]
                from_loc = int(column[5])
        
        elif variation_type == 'translocation':
            if len(column) < 5:
                print("Translocation must have contig_from and from_loc")
                print("contig_from and or from_loc is missing in"\
                        " line {0}".format(line_number))
                print("Skipping sv")
                use_sv = False
            else:
                contig_from = column[4]
                from_loc = int(column[5])
        
        print(line_number, variation_type)
        if use_sv:
            variations.append(
                        create_sv( 
                                    variation_type, 
                                    contig_name, 
                                    position, 
                                    length, 
                                    from_loc, 
                                    contig_from, 
                                    nr_duplications
                                )
                            )
        # for variation in sv:
        #     print(variation)
        # if sv:
        #     variations[ contig_name ].extend( sv )
        # else:
        #     print( "warning: Ignored bad variation on line: {0}".format( line_number ) )
    for variation in variations:
        print(variation)
    return variations



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
                    help='The fasta identifier separator, default is "|".'
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
    
    normal_contigs = Fasta(
                            normal_contig_file,
                            key_fn = lambda x: x.split( delimiter )[ field_index ].strip( )
                            )
    
    contigs = normal_contigs.keys()
    
    variations = read_variations( variation_file, contigs )
    # print(variations)
    sys.exit()
    
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

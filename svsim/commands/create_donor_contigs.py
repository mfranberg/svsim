
from __future__ import unicode_literals, print_function

import sys
import os
import click

from pprint import pprint as pp

from pyfasta import Fasta

from svsim.utils import (read_variations, write_donor_contigs)
from svsim.variations import (check_variations)
from svsim.warnings import OverlappingFeaturesError
from svsim import (init_log, vcf)

import logging

@click.command()
@click.argument('normal_contig_file',
                    type=click.Path(exists=True),
)
@click.argument('variation_file',
                    type=click.File('r'),
)
@click.argument('outfile',
                    type=click.File('w'),
)
@click.option('-d', '--delimiter',
                    default='|',
                    help='The fasta identifier separator, default is "|".'
)
@click.option('-n', '--genome_name',
                    default='donor-genome',
                    help="The name of the genome."
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
@click.option('-l', '--logfile',
                    type=click.Path(exists=False),
                    help="Path to log file. If none logging is "\
                          "printed to stderr."
)
@click.option('--loglevel',
                    type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 
                                        'CRITICAL']),
                    default='INFO',
                    help="Set the level of log output."
)
def create_donor_contigs(normal_contig_file, variation_file, outfile, delimiter,
    genome_name, field_index, chrom, vcf_file, logfile, loglevel):
    """
    Creates the donor contigs with structural variations.
    
    The format of variation_file are lines consisting of the following: 
    contig_name insertion start length, contig_name deletion start length, 
    contig_name duplication start length to or contig_name translocation start length to.
    """
    
    logger = logging.getLogger("svsim.create_donor_contigs")
    
    init_log(logger, logfile, loglevel)
    
    if vcf_file and chrom:
        vcf_file = vcf.open_vcf_file
        vcf_file.set_chrom(args.chrom)
    
    logger.debug("Getting original sort order of contigs from contig file")
    sorted_contigs = []
    with open(normal_contig_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                sorted_contigs.append(
                    line.rstrip().split(delimiter)[0][1:]
                )
    
    logger.info("Reading reference genome")
    normal_contigs = Fasta(
                            normal_contig_file,
                            key_fn = lambda x: x.split(delimiter)[field_index].strip( )
                            )
    logger.debug("Reading reference genome done.")
    
    logger.info("Reading variations")
    variations = read_variations(variation_file, sorted_contigs, logger)
    
    logger.info("Check if variations overlap")    
    for contig in variations:
        try:
            variations[contig] = check_variations(variations[contig], logger)
        except OverlappingFeaturesError as e:
            logger.critical("Interval at contig {0}, position {1} is overlapping "\
                    "another interval. Pleast check your variations file."\
                    " Exiting.".format(e.contig, e.position)
                    )
            sys.exit(1)
            
    logger.info("Write donor contigs")
    write_donor_contigs(normal_contigs, variations, sorted_contigs, 
                        genome_name, outfile)
    

if __name__ == '__main__':
    create_donor_contigs()


from __future__ import print_function

import os
import sys

import click
import logging
import subprocess
import tempfile

import pysam

import svsim.log as log

##
# Converts a sam file to a bam file using pysam.
#
# @param sam_path Path of the .sam file.
# @param bam_path Path of the resulting .bam file.
#
def sam_to_bam(sam_path, bam_path):
    sam_file = pysam.AlignmentFile( sam_path, "r" )
    bam_file = pysam.AlignmentFile( bam_path, "wb", template = sam_file )
    for alignment in sam_file:
        bam_file.write( alignment )

##
# Maps the given paired end reads using bwa, and writes a
# sorted .bam file in the given output file.
#
# @param pe1_path Path of the first reads.
# @param pe2_path Path of the second reads.
# @param genome_path Path to the reference genome.
# @param output_path Path of the output file without extension ".bam".
#
@click.command()
@click.argument('genome',
                    type=click.Path(exists=True),
)
@click.argument('output',
                    type=click.Path(exists=False),
)
@click.option('-pe1', '--pe1_path',
                    type=click.Path(exists=True),
                    help='Path to the first pairs.'
)
@click.option('-pe2', '--pe2_path',
                    type=click.Path(exists=True),
                    help='Path to the second pairs.'
)
@click.option('-l', '--logfile',
                    type=click.Path(exists=False),
                    default='./svsim_log.txt',
                    help='Path to log file.'
)
def map_reads(pe1_path, pe2_path, genome, output, logfile):
    """
    Map reads to a given genome.
    """
    log.init_log( logfile )
    if not (pe1_path and pe1_path):
        print('Please provide paths to the read files with -pe1 and -pe2', file=sys.stderr)
        sys.exit()
    
    
    work_dir = tempfile.mkdtemp( )
    genome_db = os.path.join( work_dir, "genome" )
    pe1_output = os.path.join( work_dir, "pe1.sai" )
    pe2_output = os.path.join( work_dir, "pe2.sai" )
    bwa_output = os.path.join( work_dir, "output.sam" )
    
    log_file = log.get_log_stream( )
    
    logging.info( "Running bwa index, aln and sampe in {0}".format( work_dir ) )
    
    subprocess.check_call( [ "bwa", "index", "-p", genome_db, genome ], stderr = log_file )
    
    with open( pe1_output, "w" ) as pe1_file:
        subprocess.check_call( [ "bwa", "aln", genome_db, pe1_path ], stdout = pe1_file, stderr = log_file )
    
    with open( pe2_output, "w" ) as pe2_file:
        subprocess.check_call( [ "bwa", "aln", genome_db, pe2_path ], stdout = pe2_file, stderr = log_file )
    
    with open( bwa_output, "w" ) as bwa_file:
        subprocess.check_call( [ "bwa", "sampe",
                           "-r", "@RG\tID:ILLUMINA\tSM:48_2\tPL:ILLUMINA\tLB:LIB1",
                           genome_db,
                           pe1_output, pe2_output,
                           pe1_path, pe2_path ], stdout = bwa_file, stderr = log_file )
    
    logging.info( "Sorting .bam file" )
    
    sam_to_bam( bwa_output, bwa_output + ".bam" )
    pysam.sort( bwa_output + ".bam", output )

if __name__ == '__main__':
    map_reads()

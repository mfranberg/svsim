
from __future__ import print_function

import os
import sys

import click
import logging
import subprocess
import tempfile

import pysam

from logbook import Logger

from svsim import init_log

logger = Logger('map_reads logger')

def sam_to_bam(sam_path, bam_path):
    """
    Converts a sam file to a bam file using pysam.
    
    Arguments:
        sam_path (str): Path of the .sam file.
        bam_path (str): Path of the resulting .bam file.
    """
    sam_file = pysam.AlignmentFile(sam_path, "r")
    bam_file = pysam.AlignmentFile(bam_path, "wb", template = sam_file)
    for alignment in sam_file:
        bam_file.write(alignment)

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
                    help="Path to log file. If none logging is "\
                          "printed to stderr."
)
def map_reads(pe1_path, pe2_path, genome, output, logfile):
    """
    Maps the given paired end reads using bwa, and writes a sorted .bam file 
    in the given output file.
    
    genome: Path to the reference genome
    output: Path of the output file without extension ".bam".
    """
    
    log = init_log(logfile)
    log.push_application()
    
    if not (pe1_path and pe1_path):
        logger.error('Please provide paths to the read files with -pe1 and -pe2')
        sys.exit(1)
    
    
    work_dir = tempfile.mkdtemp()
    genome_db = os.path.join(work_dir, "genome")
    pe1_output = os.path.join(work_dir, "pe1.sai")
    pe2_output = os.path.join(work_dir, "pe2.sai")
    bwa_output = os.path.join(work_dir, "output.sam")
        
    logger.info("Running bwa index, aln and sampe in {0}".format(work_dir))
    
    subprocess.check_call(["bwa", "index", "-p", genome_db, genome], stderr = log_file)
    
    with open(pe1_output, "w") as pe1_file:
        subprocess.check_call(["bwa", "aln", genome_db, pe1_path], stdout = pe1_file, stderr = log_file)
    
    with open(pe2_output, "w") as pe2_file:
        subprocess.check_call(["bwa", "aln", genome_db, pe2_path], stdout = pe2_file, stderr = log_file)
    
    with open(bwa_output, "w") as bwa_file:
        subprocess.check_call(["bwa", "sampe",
                           "-r", "@RG\tID:ILLUMINA\tSM:48_2\tPL:ILLUMINA\tLB:LIB1",
                           genome_db,
                           pe1_output, pe2_output,
                           pe1_path, pe2_path], stdout = bwa_file, stderr = log_file)
    
    logger.info("Sorting .bam file")
    
    sam_to_bam(bwa_output, bwa_output + ".bam")
    pysam.sort(bwa_output + ".bam", output)

if __name__ == '__main__':
    map_reads()

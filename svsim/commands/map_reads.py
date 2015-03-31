
from __future__ import print_function

import os
import sys

import click
import logging
import subprocess
import tempfile
from pprint import pprint as pp

import pysam

from svsim import init_log, get_log_stream

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
                    help='Path to the first pairs.',
                    required=True
)
@click.option('-pe2', '--pe2_path',
                    type=click.Path(exists=True),
                    help='Path to the second pairs.',
                    required=True
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
def map_reads(pe1_path, pe2_path, genome, output, logfile, loglevel):
    """
    Maps the given paired end reads using bwa, and writes a sorted .bam file 
    in the given output file.
    
    genome: Path to the reference genome
    output: Path of the output file without extension ".bam".
    """
    
    logger = logging.getLogger("svsim.map_reads")
    
    init_log(logger, logfile, loglevel)
    
    log_stream = get_log_stream(logger)
     
    logger.debug("Creating temporary directory.")
    work_dir = tempfile.mkdtemp()
    logger.debug("Temp dir is: {0}.".format(work_dir))
    genome_db = os.path.join(work_dir, "genome")
    bwa_output = os.path.join(work_dir, "output.sam")
    logger.debug("BWA output: {0}.".format(bwa_output))
   
    if not os.path.exists( genome_db + ".bwt" ):
        logger.info("Running bwa index and mem in {0}".format(work_dir))
        
        bwa_index_call = [
                            "bwa", 
                            "index", 
                            "-p", 
                            genome_db, 
                            genome
                        ]
        
        logger.info("Running bwa index with command {0}".format(' '.join(bwa_index_call)))
        
        subprocess.check_call(bwa_index_call, stderr = log_stream)
        
        logger.debug("BWA index done")
    
    bwa_mem_call =["bwa", "mem", genome_db, pe1_path, pe2_path]
    
    logger.debug("Running BWA align read pairs with command:{0}.".format(
                    ' '.join(bwa_mem_call)
                    )
                )
    
    with open(bwa_output, "w") as bwa_file:
        subprocess.check_call(
                                bwa_mem_call, 
                                stdout = bwa_file, 
                                stderr = log_stream
                            )
    logger.debug("Bwa mem done.") 
    
    logger.info("Converting .sam to .bam.")
    sam_to_bam(bwa_output, bwa_output + ".bam")
    logger.debug("Converting .sam to .bam done.")
    logger.info("Sorting .bam file")
    
    pysam.sort(bwa_output + ".bam", output)
    
    logger.debug("Sorting .bam file done.")

if __name__ == '__main__':
    map_reads()

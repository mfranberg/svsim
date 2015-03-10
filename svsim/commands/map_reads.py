
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
    
    if not (pe1_path != None and pe1_path != None):
        logger.error('Please provide paths to the read files with -pe1 and -pe2')
        sys.exit(1)
        
    
    logger.debug("Creating temporary directory.")
    work_dir = tempfile.mkdtemp()
    logger.debug("Temp dir is: {0}.".format(work_dir))
    genome_db = os.path.join(work_dir, "genome")
    pe1_output = os.path.join(work_dir, "pe1.sai")
    pe2_output = os.path.join(work_dir, "pe2.sai")
    bwa_output = os.path.join(work_dir, "output.sam")
    logger.debug("BWA output: {0}.".format(bwa_output))
    
    logger.info("Running bwa index, aln and sampe in {0}".format(work_dir))
    
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
    
    bwa_aln_1_call =["bwa", "aln", genome_db, pe1_path]
    
    logger.debug("Running BWA align on first read pairs with command:{0}.".format(
                    ' '.join(bwa_aln_1_call)
                    )
                )
    
    with open(pe1_output, "w") as pe1_file:
        subprocess.check_call(
                                bwa_aln_1_call, 
                                stdout = pe1_file, 
                                stderr = log_stream
                            )
    logger.debug("First read pair done.")
    
    bwa_aln_2_call =["bwa", "aln", genome_db, pe2_path]
    
    logger.debug("Running BWA align on second read pairs with command:{0}.".format(
        ' '.join(bwa_aln_2_call)
        )
    )
    
    with open(pe2_output, "w") as pe2_file:
        subprocess.check_call(
                                bwa_aln_2_call, 
                                stdout = pe2_file, 
                                stderr = log_stream
                            )
    
    logger.debug("Second read pair done.")
    
    bwa_sampe_call = [
        "bwa", 
        "sampe",
        "-r", 
        "@RG\tID:ILLUMINA\tSM:48_2\tPL:ILLUMINA\tLB:LIB1",
        genome_db,
        pe1_output, 
        pe2_output,
        pe1_path, 
        pe2_path
    ]
    
    logger.debug("Merging read pairs with command:{0}.".format(
        ' '.join(bwa_sampe_call)
        )
    )
    
    with open(bwa_output, "w") as bwa_file:
        subprocess.check_call(
                                bwa_sampe_call, 
                                stdout = bwa_file, 
                                stderr = log_stream
                            )
    
    logger.debug("Merging read pairs done.")
    
    logger.info("Converting .sam to .bam.")
    sam_to_bam(bwa_output, bwa_output + ".bam")
    logger.debug("Converting .sam to .bam done.")
    logger.info("Sorting .bam file")
    
    pysam.sort(bwa_output + ".bam", output)
    
    logger.debug("Sorting .bam file done.")

if __name__ == '__main__':
    map_reads()

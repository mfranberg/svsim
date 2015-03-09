import os
import subprocess
import tempfile
import fnmatch
import logging

from pkg_resources import resource_filename

import svsim.log as log
from svsim.util import calculate_num_reads, get_genome_length

class IReadSimulator(object):
    """
    Abstract base class for read simulators.
    """
    def __init__(self):
        super(IReadSimulator, self).__init__()
        ##
        # Genome coverage.
        #
        self.coverage = 10
        ##
        # Read length.
        #
        self.read_length = 100
        
        ##
        # Mean distance between read pairs.
        #
        
        self.mean = 500
        ##
        # Standard deviation for distance between pairs.
        # 
        self.std = 50
    
    def simulate(genome_path, output_path):
        """
        Simulates reads from the given genome and writes them in the 
        given output file.
        
        Arguments:
            genome_path (str): Path to file that contains a genome to simulate 
                                reads from.
            output_path (str): Path to file to write simulated reads to.
        """
        raise NotImplementedError( "ReadSimulator: Unimplemented abstract method." )


class DwgsimSimulator(IReadSimulator):
    """Simulates reads using dwgsim with Illumina-like properties."""
    def __init__(self):
        super(DwgsimSimulator, self).__init__()
    
    def simulate(self, genome_path, output_file):
        output_prefix = os.path.join(tempfile.mkdtemp(), "dwgsim")

        read_error_rate = "{0}-{0}".format(self.read_error)
        subprocess.check_call(
            [ 
                "dwgsim",
                "-d", str(int(self.mean)),
                "-s", str(int( self.std)),
                "-C", str(self.coverage),
                "-1", str(int(self.read_length)),
                "-2", str(int(self.read_length)),
                "-e", read_error_rate,
                "-E", read_error_rate,
                "-c", "0",
                "-r", "0",
                "-R", "0",
                "-y", "0",
                genome_path,
                output_prefix 
            ], 
            stdout = open("/dev/null"), 
            stderr = open("/dev/null"))
        
        subprocess.check_call(
            [
                "cp", 
                output_prefix + ".bwa.read1.fastq", 
                output_file + "_pe1.fa"
            ]
        )
        subprocess.check_call(
            [
                "cp", 
                output_prefix + ".bwa.read2.fastq", 
                output_file + "_pe2.fa"
            ]
        )

class MetaSimSimulator(IReadSimulator):
    """Simulates reads using metasim."""
    def __init__(self):
        super(MetaSimSimulator, self).__init__()
        self.error_model = resource_filename( "svsim", "data/errormodel-100bp.mconf" )
    
    def simulate(self, genome_path, output_file):
        if not self.read_length == 100:
            raise ValueError("Read length must be 100 for metasim.")
        
        genome_length = get_genome_length(genome_path)
        num_reads = calculate_num_reads(self.coverage, self.read_length, genome_length)
        output_dir = tempfile.mkdtemp()
        
        logging.info( "Starting metasim:" )
        subprocess.check_call(
            [
                 "MetaSim", 
                 "cmd",
                 "-r", 
                 str(num_reads), 
                 "-m",
                 "-g", 
                 self.error_model,
                 "-2", 
                 self.error_model,
                 "--empirical-pe-probability", 
                 "100",
                 "--clones-mean", 
                 str(self.mean),
                "--clones-param2", 
                str(self.std),
                "-d", 
                output_dir,
                genome_path 
            ], 
            stdout = log.get_log_stream()
        )
        
        # Metasim outputs reads for each contig, gather them into one file
        found_match = False
        for i, metasim_output_file in enumerate(os.listdir(output_dir)):
            metasim_output_name = os.path.basename(genome_path).split(".")[0] + "-Empirical.*fna"
            if fnmatch.fnmatch(metasim_output_file, metasim_output_name):
                metasim_output_file_absolute = os.path.join(output_dir, metasim_output_file)
                convert_to_pe(
                                metasim_output_file_absolute, 
                                output_file + "_pe1.fa", 
                                output_file + "_pe2.fa", 
                                append = i > 0
                            )
                found_match = True
        
        if not found_match:
            raise RuntimeError( "MetaSim could not simulate reads, is the genome too short?" )

def convert_to_pe(metasim_output_path, pe1_path, pe2_path, append = False ):
    """
    Converts a metasim output fasta file in which all reads are
    in one file, to two separate files containg the first and second
    pair respectively.
    
    Arguments:
        metasim_output_path (str): An output file from metasim.
        pe1_path (str): The path to fasta file that will contain the first read.
        pe2_path (str): The path to fasta file that will contain the second read.
        append (bool): Determines whether to append to the output files instead.
    """
    open_type = "w"
    if append:
        open_type = "a"
    
    output_pe1 = open(pe1_path, open_type)
    output_pe2 = open(pe2_path, open_type)
    output_file = None
    
    with open(metasim_output_path, "r") as metasim_output_file:
        for line in metasim_output_file:
            columns = line.split()
            
            if columns[0].startswith(">"):
                if output_file:
                    output_file.write("\n")
                
                if columns[0].endswith(".1"):
                    output_file = output_pe1
                else:
                    output_file = output_pe2
                
                output_file.write(line)
            else:
                output_file.write(line.strip())
    
    output_pe1.close()
    output_pe2.close()

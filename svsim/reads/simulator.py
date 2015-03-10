import os
import subprocess
import tempfile
import fnmatch
import shutil

from pkg_resources import resource_filename

from svsim.util import (calculate_num_reads, get_genome_length)
from svsim.log import get_log_stream

class IReadSimulator(object):
    """
    Abstract base class for read simulators.
    """
    def __init__(self, logger):
        super(IReadSimulator, self).__init__()
        
        self.logger = logger
        self.log_stream = get_log_stream(logger)
        
        self.logger.info("Initiating read simulator")
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
    
    def simulate(genome_path, output_path, second_genome=None):
        """
        Simulates reads from the given genome and writes them in the 
        given output file.
        
        Arguments:
            genome_path (str): Path to file that contains a genome to simulate 
                                reads from.
            output_path (str): Path to file to write simulated reads to.
        """
        raise NotImplementedError( "ReadSimulator: Unimplemented abstract method." )
    
    def __repr__(self):
        return "IReadSimulator(coverage={0}, read_length={1}, mean={2}, "\
               "std={3})".format(
                   self.coverage,
                   self.read_length,
                   self.mean,
                   self.std
               )

class DwgsimSimulator(IReadSimulator):
    """Simulates reads using dwgsim with Illumina-like properties."""
    def __init__(self, logger):
        super(DwgsimSimulator, self).__init__(logger)
    
    def simulate(self, genome_path, output_file, second_genome_path=None):
        output_prefix = os.path.join(tempfile.mkdtemp(), "dwgsim")
        
        read_error_rate = "{0}-{0}".format(self.read_error)
        self.logger.info("Starting dwgsim")
        dwgsim_call = [
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
        ]
        self.logger.debug("Running dwgsim with call {0}".format(
            ' '.join(dwgsim_call)
            )
        )
        
        subprocess.check_call(
            dwgsim_call, 
            stdout = self.log_stream, 
            stderr = self.log_stream
        )
        
        if second_genome_path:
            second_output_prefix = os.path.join(tempfile.mkdtemp(), "dwgsim2")
            self.logger.info("Starting dwgsim for second genome")
            dwgsim_call2 = [
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
                    second_genome_path,
                    second_output_prefix
            ]
            self.logger.debug("Running dwgsim with call {0}".format(
                ' '.join(dwgsim_call2)
                )
            )
            subprocess.check_call(
                dwgsim_call2, 
                stdout = self.log_stream, 
                stderr = self.log_stream
            )
            
            first_pair_files = [
                output_prefix + ".bwa.read1.fastq",
                second_output_prefix + ".bwa.read1.fastq"
            ]
            second_pair_files = [
                output_prefix + ".bwa.read2.fastq",
                second_output_prefix + ".bwa.read2.fastq"
            ]
        else:
            first_pair_files = [
                output_prefix + ".bwa.read1.fastq"
            ]
            second_pair_files = [
                output_prefix + ".bwa.read2.fastq"
            ]
            
        with open(output_file + "_pe1.fa", 'w') as outfile_1:
            for infile in first_pair_files:
                shutil.copyfileobj(open(infile), outfile_1)
        
        with open(output_file + "_pe2.fa", 'w') as outfile_2:
            for infile in first_pair_files:
                shutil.copyfileobj(open(infile), outfile_2)
            
    
    def __repr__(self):
        return "DwgsimSimulator(coverage={0}, read_length={1}, mean={2}, "\
               "std={3})".format(
                   self.coverage,
                   self.read_length,
                   self.mean,
                   self.std
               )

class MetaSimSimulator(IReadSimulator):
    """Simulates reads using metasim."""
    def __init__(self, logger):
        super(MetaSimSimulator, self).__init__(logger)
        self.error_model = resource_filename( "svsim", "data/errormodel-100bp.mconf" )
    
    def simulate(self, genome_path, output_file, second_genome=None):
        if not self.read_length == 100:
            self.logger.critical("Read length must be 100 for metasim.")
            raise ValueError()
        
        genome_length = get_genome_length(genome_path)
        num_reads = calculate_num_reads(self.coverage, self.read_length, genome_length)
        output_dir = tempfile.mkdtemp()
        
        self.logger.info( "Starting metasim:" )
        meta_sim_call = [
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
        ]
        self.logger.debug("Running metasim with call {0}".format(
            ' '.join(meta_sim_call)
            )
        )
        subprocess.check_call(
            meta_sim_call, 
            stdout = self.log_stream,
            stderr = self.log_stream
        )
        
        # Metasim outputs reads for each contig, gather them into one file
        found_match = False
        for i, metasim_output_file in enumerate(os.listdir(output_dir)):
            metasim_output_name = os.path.basename(genome_path).split(".")[0] + "-Empirical.*fna"
            if fnmatch.fnmatch(metasim_output_file, metasim_output_name):
                metasim_output_file_absolute = os.path.join(output_dir, metasim_output_file)
                self.convert_to_pe(
                                metasim_output_file_absolute, 
                                output_file + "_pe1.fa", 
                                output_file + "_pe2.fa", 
                                append = i > 0
                            )
                found_match = True
        
        if not found_match:
            self.logger.critical("MetaSim could not simulate reads, "\
                                 "is the genome too short?")
            raise RuntimeError()
    
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
    
    def __repr__(self):
        return "MetaSimSimulator(coverage={0}, read_length={1}, mean={2}, "\
               "std={3})".format(
                   self.coverage,
                   self.read_length,
                   self.mean,
                   self.std
               )
    


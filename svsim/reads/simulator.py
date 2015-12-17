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
        
        self.logger.debug("Merging files for first pair")
            
        with open(output_file + "_pe1.fa", 'w') as outfile_1:
            for infile in first_pair_files:
                shutil.copyfileobj(open(infile), outfile_1)
        
        self.logger.debug("Merging files for second pair")
        
        with open(output_file + "_pe2.fa", 'w') as outfile_2:
            for infile in second_pair_files:
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
    


##
# Simulates reads using a lognormal distribution.
#
class LogNSimulator( IReadSimulator ):
    ##
    # Constructor.
    #
    def __init__(self, logger):
        super(LogNSimulator, self).__init__(logger)
        if self.mean < 10 and  1 < self.std <  2: # Well defined distribution
            pass
        else: # user probably forgot that mu and std needs to be specified in log base
            self.logger.critical("mu and std needs to be specified in log base (usually mu < 10 and 1 < sigma < 2).\
                    You specified mean: {0}, sigma: {1}".format(self.mean, self.std))
            raise ValueError()


    ##
    # @see IReadSimulator.simulate
    #
    def simulate(self, genome_path, output_prefix):
        if not self.read_length == 100:
            self.logger.critical("Read length must be 100 for metasim.")
            raise ValueError()

        genome_sequence = get_genome_sequence( genome_path )
        genome_length = len(genome_sequence)
        dna_library = DNAseq(self.read_length, self.coverage, mean=self.mean, stddev=self.std, distribution='lognormal')
        dna_library.simulate_pe_reads(genome_sequence)

        reads1 = open(os.path.join(output_prefix, "_pe1.fa"),'w')
        reads2 = open(os.path.join(output_prefix, "_pe2.fa"),'w')
        i=0
        for read in dna_library.fasta_format():
            if i%2==0:
                reads1.write(read)
            else:
                reads2.write(read)
            i+=1

        # output_dir = tempfile.mkdtemp( )

        # subprocess.call( [ "MetaSim", "cmd",
        #                    "-r", str( num_reads ), 
        #                    "-m",
        #                    "-g", self.error_model,
        #                    "-2", self.error_model,
        #                    "--empirical-pe-probability", "100",
        #                    "--clones-mean", str( self.mean ),
        #                    "--clones-param2", str( self.std ),
        #                    "-d", output_dir,
        #                    genome_path ], stdout = open( "/dev/null", "w" ) )


class PairedEndRead(object):
    """docstring for PairedEndRead"""
    def __init__(self,distribution = 'normal',mean=None,sigma=None,read_length= None, min_size=None, max_size=None):
        super(PairedEndRead, self).__init__()
        self.distribution = distribution
        self.mean = mean
        self.sigma = sigma
        self.read_length = read_length
        self.min_size = min_size
        self.max_size = max_size

    def generate(self, reference_sequence, read_index):
        if self.distribution == 'normal':
            self.fragment_length = int(random.gauss(self.mean,self.sigma))
        elif self.distribution == 'uniform':
            self.fragment_length = int(random.uniform(self.min_size,self.max_size))
        elif self.distribution == 'lognormal':
            self.fragment_length = int(lognormal(self.mean, self.sigma)[0]) # one sample at a time to conform with the implementation...

        if self.fragment_length >= len(reference_sequence): 
            raise Exception("To short reference sequence length for \
                simulated read. \nRead fragment: {0}\nTranscript \
                length:{1}".format(self.fragment_length,len(reference_sequence)))
        
        self.start_pos = random.randrange(len(reference_sequence) - self.fragment_length)
        self.read1 = reference_sequence[self.start_pos : self.start_pos + self.read_length]
        self.read2 = reverse_complement(reference_sequence[self.start_pos + self.fragment_length - self.read_length : self.start_pos+self.fragment_length])
        self.reference_accession = "reference_genome"
        self.read_index = read_index

    def fastq_format(self):
        r1= '@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\
        \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
            self.reference_accession,self.read1,'J'*self.read_length)
        r2= '@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\
        \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
            self.reference_accession,self.read2,'J'*self.read_length)        
        yield r1
        yield r2

    def fasta_format(self):
        r1= '>HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\n'.format(self.read_index,
            self.reference_accession,self.read1)
        r2= '>HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\n'.format(self.read_index,
            self.reference_accession,self.read2)        
        yield r1
        yield r2


class DNAseq(object):
    """docstring for DNAseq"""
    def __init__(self,read_length, coverage, mean=None,stddev=None, min_size=None, max_size = None, distribution='normal'):
        super(DNAseq, self).__init__()
        self.distribution = distribution
        
        if self.distribution == 'normal' or self.distribution == "lognormal":
            self.mean = mean
            self.stddev = stddev
        elif self.distribution == 'uniform':
            self.min_size = min_size
            self.max_size = max_size

        self.read_length = read_length
        self.coverage = coverage

    def simulate_pe_reads(self, genome_sequence):
        """
        Arguments:
        """
        genome_length = len(genome_sequence)
        number_of_reads = calculate_num_reads(self.coverage, self.read_length, genome_length)  # Specifies the number of simulated read pairs (related to insertion size length of genome and coverage
    
        self.reads = []
        
        
        for i in range(number_of_reads):
            if self.distribution == 'normal':
                read_pair = PairedEndRead(distribution=self.distribution, mean=self.mean, sigma=self.stddev, read_length=self.read_length)
            elif self.distribution == 'lognormal':
                read_pair = PairedEndRead(distribution=self.distribution, mean=self.mean, sigma=self.stddev, read_length=self.read_length)
            elif self.distribution == 'uniform':
                read_pair = PairedEndRead(distribution=self.distribution, min_size=self.min_size,max_size=self.max_size,read_length=self.read_length)
     
            read_pair.generate(genome_sequence, i)
            self.reads.append(read_pair)


    def fastq_format(self):
        for pe_read in self.reads:
            for mate in pe_read.fastq_format():
                yield mate

    def fasta_format(self):
        for pe_read in self.reads:
            for mate in pe_read.fasta_format():
                yield mate

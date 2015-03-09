import os
import subprocess
import tempfile

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
        output_prefix = os.path.join( tempfile.mkdtemp( ), "dwgsim" )

        read_error_rate = "{0}-{0}".format( self.read_error )
        subprocess.check_call(
            [ 
                "dwgsim",
                "-d", str( int( self.mean ) ),
                "-s", str( int( self.std ) ),
                "-C", str( self.coverage ),
                "-1", str( int( self.read_length ) ),
                "-2", str( int( self.read_length ) ),
                "-e", read_error_rate,
                "-E", read_error_rate,
                "-c", "0",
                "-r", "0",
                "-R", "0",
                "-y", "0",
                genome_path,
                output_prefix 
            ], 
            stdout = open( "/dev/null" ), 
            stderr = open( "/dev/null" ) )
        
        subprocess.check_call( [ "cp", output_prefix + ".bwa.read1.fastq", output_file + "_pe1.fa" ] )
        subprocess.check_call( [ "cp", output_prefix + ".bwa.read2.fastq", output_file + "_pe2.fa" ] )

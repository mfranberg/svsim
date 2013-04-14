import fnmatch
import os
import subprocess
import tempfile

from pkg_resources import resource_filename

from svsim.util import calculate_num_reads, get_genome_length
from svsim.reads.isim import IReadSimulator

##
# Simulates reads using metasim.
#
class MetaSimSimulator( IReadSimulator ):
    ##
    # Constructor.
    #
    def __init__(self):
        self.error_model = resource_filename( "svsim", "data/errormodel-100bp.mconf" )

    ##
    # @see IReadSimulator.simulate
    #
    def simulate(self, genome_path, output_file):
        if not self.read_length == 100:
            raise ValueError( "Read length must be 100 for metasim." )

        genome_length = get_genome_length( genome_path )
        num_reads = calculate_num_reads( self.coverage, self.read_length, genome_length )
        output_dir = tempfile.mkdtemp( )

        subprocess.check_call( [ "MetaSim", "cmd",
                           "-r", str( num_reads ), 
                           "-m",
                           "-g", self.error_model,
                           "-2", self.error_model,
                           "--empirical-pe-probability", "100",
                           "--clones-mean", str( self.mean ),
                           "--clones-param2", str( self.std ),
                           "-d", output_dir,
                           genome_path ], stdout = open( "/dev/null", "w" ) )

        # Metasim outputs reads for each contig, gather them into one file
        for i, metasim_output_file in enumerate( os.listdir( output_dir ) ):
            metasim_output_name = os.path.basename( genome_path ).split( "." )[ 0 ] + "-Empirical.*fna"
            if fnmatch.fnmatch( metasim_output_file, metasim_output_name ):
                metasim_output_file_absolute = os.path.join( output_dir, metasim_output_file )
                convert_to_pe( metasim_output_file_absolute, output_file + "_pe1.fa", output_file + "_pe2.fa", append = i > 0 )

##
# Converts a metasim output fasta file in which all reads are
# in one file, to two separate files containg the first and second
# pair respectively.
#
# @param metasim_output_path An output file from metasim.
# @param pe1_path The path to fasta file that will contain the first read.
# @param pe2_path The path to fasta file that will contain the second read.
# @param append Determines whether to append to the output files instead.
#
def convert_to_pe(metasim_output_path, pe1_path, pe2_path, append = False ):
    open_type = "w"
    if append:
        open_type = "a"

    output_pe1 = open( pe1_path, open_type )
    output_pe2 = open( pe2_path, open_type )
    output_file = None

    with open( metasim_output_path, "r" ) as metasim_output_file:
        for line in metasim_output_file:
            columns = line.split( )

            if columns[ 0 ].startswith( ">" ):
                if output_file:
                    output_file.write( "\n" )

                if columns[ 0 ].endswith( ".1" ):
                    output_file = output_pe1
                else:
                    output_file = output_pe2

                output_file.write( line )
            else:
                output_file.write( line.strip( ) )

    output_pe1.close( )
    output_pe2.close( )


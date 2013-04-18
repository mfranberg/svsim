import os
import subprocess
import tempfile

from svsim.reads.isim import IReadSimulator

##
# Simulates reads using dwgsim with Illumina-like
# properties.
#
class DwgsimSimulator( IReadSimulator ):
    def simulate(self, genome_path, output_file):
        output_prefix = os.path.join( tempfile.mkdtemp( ), "dwgsim" )

        read_error_rate = "{0}-{0}".format( self.read_error )
        subprocess.check_call( [ "dwgsim",
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
                           output_prefix ], stdout = open( "/dev/null" ), stderr = open( "/dev/null" ) )

        subprocess.check_call( [ "cp", output_prefix + ".bwa.read1.fastq", output_file + "_pe1.fa" ] )
        subprocess.check_call( [ "cp", output_prefix + ".bwa.read2.fastq", output_file + "_pe2.fa" ] )

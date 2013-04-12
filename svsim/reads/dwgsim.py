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

        subprocess.call( [ "dwgsim",
                           "-d", str( int( self.mean ) ),
                           "-s", str( int( self.std ) ),
                           "-C", str( self.coverage ),
                           "-1", str( int( self.read_length ) ),
                           "-2", str( int( self.read_length ) ),
                           "-c", "0",
                           genome_path,
                           output_prefix ], stdout = open( "/dev/null" ), stderr = open( "/dev/null" ) )

        subprocess.call( [ "cp", output_prefix + ".bwa.read1.fastq", output_file + "_pe1.fa" ] )
        subprocess.call( [ "cp", output_prefix + ".bwa.read2.fastq", output_file + "_pe2.fa" ] )

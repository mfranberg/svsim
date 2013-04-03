##
# Abstract base class for read simulators.
#
#
class IReadSimulator:
    def __init__(self):
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

    ##
    # Simulates reads from the given genome and writes
    # them in the given output file.
    #
    # @param genome_path File that contains a genome to simulate reads from.
    # @param output_path File to write simulated reads to.
    #
    def simulate(genome_path, output_path):
        raise Exception( "ReadSimulator: Unimplemented abstract method." )

import random
from StringIO import StringIO

##
# Represents a insertion of a sequence in the reference
# genome.
#
class Insertion:
    ##
    # Constructor.
    #
    # @param pos The 0-based position of the base pair before
    #            the insertion in the refrence genome.
    # @param length The length of the inserted sequence.
    # @param from_loc The 0-based position from where the
    #                 inserted sequence is taken. Can be -1
    #                 then a random sequence is inserted.
    #
    def __init__(self, pos, length, from_loc):
        self.pos = pos
        self.length = length
        self.from_loc = from_loc
        if from_loc < 0:
            self.sequence =  ''.join( random.choice( "ACGT" ) for i in range( self.length ) )

    def get_sequence(self, normal_genome):
        if self.from_loc >= 0:
            return normal_genome[ (self.from_loc + 1):(self.from_loc + 1 + self.length) ] 
        else:
            return self.sequence

    def get_delta(self):
        return 0

##
# Represents a deletion of a sequence in the
# genome starting at pos and ending at pos + length.
#
class Deletion:
    ##
    # Constructor.
    #
    # @param pos The 0-based position of the base pair before
    #            the deletion in the reference genome.
    # @param length Length of the deletion.
    #
    def __init__(self, pos, length):
        self.pos = pos
        self.length = length

    def get_sequence(self, normal_genome):
        return [ ]

    def get_delta(self):
        return self.length

##
# Represents a chunk of the reference genome.
#
class NullVariation:
    ##
    # Constructor.
    #
    # @param pos The 0-based start position of the chunk (inclusive).
    # @param length The length of the chunk.
    #
    def __init__(self, pos, length):
        self.pos = pos
        self.length = length

    def get_sequence(self, normal_genome):
        return normal_genome[ self.pos:(self.pos + self.length) ]

    def get_delta(self):
        return self.length

##
# Chunks the genome according the variations. So that each
# part of the genome is either has no varation (Null) or has
# a variation (Insertion or Deletion). These chunks are then
# sorted, and the sequence for the new genome can easily be found
# by calling get_sequence on each chunk in order.
#
# @param variations List of variations.
# @param genome_length Length of the genome.
#
# @return A list of chunks.
#
def create_chunks(variations, genome_length):
    variations = sorted( variations, key = lambda v: v.pos )
    chunks = [ ]

    pos = 0
    for v in variations:
        # Variations always point to the base pair before, and
        # we want to include that in the null chunk.
        chunks.append( NullVariation( pos, v.pos - pos + 1 ) )
        chunks.append( v )
       
        pos += ( v.pos - pos + 1 ) + v.get_delta( )

    chunks.append( NullVariation( pos, genome_length - pos ) )

    return chunks

##
# Generates the genome with the given variations and returns it.
#
# @param normal_genome The sequence for the normal genome.
# @param variations A list of variations.
#
# @return The sequence for the new genome.
#
def create_indel_genome(normal_genome, variations):
    mutated_genome = StringIO( )
    chunks = create_chunks( variations, len( normal_genome ) )
    for chunk in chunks:
        sequence = chunk.get_sequence( normal_genome )
        
        if len( sequence ) > 0:
            mutated_genome.write( sequence )

    return mutated_genome.getvalue( )

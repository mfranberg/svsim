import random
from StringIO import StringIO

class StructuralVariation(object):
    """Base class for structural variations"""
    def __init__(self, contig, pos, length):
        """
        Arguments:
            contig (str): The contig where the SV is located.
            pos (int): The 0-based position of the base pair before
                       the sv.
            length (int): The length of the event.
        """
        super(StructuralVariation, self).__init__()
        self.contig = contig
        self.pos = pos
        self.length = length
    
    def get_sequence(self, normal_genome):
        return normal_genome[ self.pos:(self.pos + self.length) ]
    
    def get_delta(self):
        return self.length
    
    def __repr__(self):
        return "StructuralVariation(contig={0}, pos={1}, length={2})".format(
            self.contig,
            self.pos,
            self.length
        )
    
    # def __str__(self):
    #     return "Contig:{0}\nPos:{1}\nLength:{2}\n".format(
    #                                                     self.contig,
    #                                                     self.pos,
    #                                                     self.length
    #                                                     )
        
class Insertion(StructuralVariation):
    """
    Represents a insertion of a sequence in the reference genome.
    If from_loc == -1 a random sequence will be inserted.
    Otherwise a sequence from the specified location is inserted
    """
    def __init__(self, contig, pos, length, from_loc = -1, from_contig=None):
        """
        Arguments:
            from_loc (int): The 0-based position from where the
                            inserted sequence is taken. Can be -1
                            then a random sequence is inserted.
        """
        super(Insertion, self).__init__(contig, pos, length)
        self.from_loc = from_loc
        if from_contig:
            self.from_contig = from_contig
        else:
            self.from_contig = contig
        
        if from_loc < 0:
            self.sequence =  ''.join( random.choice( "ACGT" ) for i in range( self.length ) )
        
    def get_sequence(self, normal_genome):
        if self.from_loc >= 0:
            return normal_genome[ (self.from_loc):(self.from_loc + self.length) ] 
        else:
            return self.sequence
        
    def get_delta(self):
        return 0
    
    def __repr__(self):
        return "Insertion(contig={0}, pos={1}, length={2}), "\
                "from_loc={3}, from_contig={4}".format(
                                                    self.contig,
                                                    self.pos,
                                                    self.length,
                                                    self.from_loc,
                                                    self.from_contig
                                                    )
    

class Deletion(StructuralVariation):
    """
    Represents a deletion of a sequence in the
    genome starting at pos and ending at pos + length.
    """
    def __init__(self, contig, pos, length):
        super(Deletion, self).__init__(contig, pos, length)
    
    def get_sequence(self, normal_genome):
        return []
    
    def __repr__(self):
        return "Deletion(contig={0}, pos={1}, length={2})".format(
            self.contig,
            self.pos,
            self.length
        )
    
class Transversion(StructuralVariation):
    """
    Represents a transversion of a sequence in the
    genome starting at pos and ending at pos + length.
    """
    def __init__(self, contig, pos, length):
        super(Transversion, self).__init__(contig, pos, length)
    
    def get_sequence(self, normal_genome):
        return []
    
    def __repr__(self):
        return "Transversion(contig={0}, pos={1}, length={2})".format(
            self.contig,
            self.pos,
            self.length
        )
    

class Duplication(StructuralVariation):
    """
    Represents a copy number variation of a sequence in the
    genome starting at pos and ending at pos + length.
    """
    def __init__(self, contig, pos, length, nr_copies=1):
        super(Duplication, self).__init__(contig, pos, length)
        self.nr_copies = nr_copies
    
    def get_sequence(self, normal_genome):
        return []
    
    def __repr__(self):
        return "Duplication(contig={0}, pos={1}, length={2}, nr_copies={3})".format(
            self.contig,
            self.pos,
            self.length,
            self.nr_copies
        )
    

class Translocation(StructuralVariation):
    """
    Represents a translocation of a sequence in the reference genome.
    This means that a part of the genome has been moved.
    """
    def __init__(self, contig, pos, length, from_contig, from_loc):
        """
        Arguments:
            from_contig (str): The contig from where the sequnce is taken.
            
            from_loc (int): The 0-based position from where the
                            inserted sequence is taken.
        """
        super(Translocation, self).__init__(contig, pos, length)
        self.from_loc = from_loc
        self.from_contig = from_contig
        
    def get_sequence(self, normal_genome):
        if self.from_loc >= 0:
            return normal_genome[ (self.from_loc):(self.from_loc + self.length) ] 
        else:
            return self.sequence
        
    def get_delta(self):
        return 0
    
    def __repr__(self):
        return "Translocation(contig={0}, pos={1}, length={2}), "\
                "from_loc={3}, from_contig={4}".format(
                                                    self.contig,
                                                    self.pos,
                                                    self.length,
                                                    self.from_loc,
                                                    self.from_contig
                                                    )
    
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

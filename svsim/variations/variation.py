import random

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
        """
        Takes the normal genome and return the appropriate segment.
        """
        return normal_genome[self.contig][self.pos:(self.pos + self.length)]

    def get_delta(self):
        """
        Returns how much we move forward on the reference genome when this variant occurs.
        """
        return self.length
    
    def __repr__(self):
        return "StructuralVariation(contig={0}, pos={1}, length={2})".format(
            self.contig,
            self.pos,
            self.length
        )

class Insertion(StructuralVariation):
    """
    Represents a insertion of a sequence in the reference genome.
    If from_loc == -1 a random sequence will be inserted.
    Otherwise a sequence from the specified location is inserted
    """
    def __init__(self, contig, pos, length, from_contig=None, from_loc = -1):
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
        
        
    def get_sequence(self, normal_genome):
        """
        The sequence of an insertion consists of the insertion followed by the
        reference sequence of the length.
        """
        if self.from_loc >= 0:
            insertion = normal_genome[self.from_contig][(self.from_loc):(self.from_loc + self.length)]
        else:
            insertion = self.sequence =  ''.join(random.choice( "ACGT" ) for i in range(self.length))
        
        return insertion
        
    def get_delta(self):
        return 0
    
    def __repr__(self):
        return "Insertion(contig={0}, pos={1}, length={2}, "\
                "from_contig={3}, from_loc={4})".format(
                                                    self.contig,
                                                    self.pos,
                                                    self.length,
                                                    self.from_contig,
                                                    self.from_loc
                                                    )
    
class Deletion(StructuralVariation):
    """
    Represents a deletion of a sequence in the
    genome starting at pos and ending at pos + length.
    """
    def __init__(self, contig, pos, length):
        super(Deletion, self).__init__(contig, pos, length)
    
    def get_sequence(self, normal_genome):
        return ''

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
        "Return the reversed sequence"
        return normal_genome[self.contig][self.pos:(self.pos + self.length)][::-1]
    
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
    def __init__(self, contig, pos, length, nr_copies=2):
        super(Duplication, self).__init__(contig, pos, length)
        self.nr_copies = nr_copies
    
    def get_sequence(self, normal_genome):
        return normal_genome[self.contig][self.pos:(self.pos + self.length)]*(self.nr_copies)
    
    def __repr__(self):
        return "Duplication(contig={0}, pos={1}, length={2}, nr_copies={3})".format(
            self.contig,
            self.pos,
            self.length,
            self.nr_copies
        )



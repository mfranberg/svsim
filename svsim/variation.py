import random
from collections import defaultdict
from interval_tree import IntervalTree

VARIATIONS = ("insertion", "deletion", "duplication", "translocation", 
                "transversion")

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
        ref_seq = normal_genome[self.contig][(self.pos):(self.pos + self.length)]
        
        if self.from_loc >= 0:
            insertion = normal_genome[self.from_contig][(self.from_loc):(self.from_loc + self.length)]
        else:
            insertion = self.sequence =  ''.join(random.choice( "ACGT" ) for i in range(self.length))
        
        return insertion + ref_seq
        
    def get_delta(self):
        return 0
    
    def __repr__(self):
        return "Insertion(contig={0}, pos={1}, length={2}), "\
                "from_contig={3}, from_loc={4}".format(
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
    

def read_variations(variation_file, contigs, logger):
    """
    Reads the variations defined in a file and returns
    them as a list of variation objects.
    
    Args:
        variation_file (file): A file that defines variations.
        contigs (list): A list with the existing contigs from the
                        reference genome
        logger (logger): An logbook object to print messages.
    
    Returns:
        variations (dict): A dictionary with contigs as keys and a list of 
                          variation objects as values.
    """
    variations = defaultdict(list)
    
    for line_number, line in enumerate(variation_file):
        use_sv = True
        column = line.split()
        contig_name = column[0]
        contig_from = contig_name
        from_loc = -1
        # Test if the contig exists
        if contig_name not in contigs:
            logger.warning("Contig does not exist in "\
                    "reference.\n Skipping sv on line {0}".format(line_number))
            use_sv = False
        
        variation_type = column[1]
        
        if variation_type not in VARIATIONS:
            logger.warning("Unknown variation: {0} at "\
                    "line number: {1}.\n Skipping sv.".format(
                                                        variation_type, 
                                                        line_number
                                                        )
                                                    )
            use_sv = False
        
        position = int(column[2])
        length = int(column[3])
        
        nr_duplications = 2
        if variation_type == 'duplication':
            if len(column) > 4:
                nr_duplications = int(column[4])
        
        elif variation_type == 'insertion':
            if len(column) > 5:
                contig_from = column[4]
                from_loc = int(column[5])
        
        elif variation_type == 'translocation':
            # We need to treat translocations as one deletion
            # and one insertion
            use_sv = False
            if len(column) < 5:
                logger.warning("Translocation must have contig_from and"\
                " from_loc. contig_from and/or from_loc is missing at "\
                "line number: {0}.\n Skipping sv.".format(
                                                        line_number
                                                        )
                                                    )
            else:
                contig_from = column[4]
                from_loc = int(column[5])
                variations[contig_from].append(
                    create_sv(
                                'deletion',
                                contig_from,
                                from_loc,
                                length, 
                    )
                )
                variations[contig_name].append(
                    create_sv(
                                'insertion',
                                contig_name,
                                position,
                                length,
                                contig_from,
                                from_loc
                    )
                )
        
        if use_sv:
            variations[contig_name].append(
                        create_sv( 
                                    variation_type, 
                                    contig_name, 
                                    position, 
                                    length, 
                                    contig_from, 
                                    from_loc, 
                                    nr_duplications
                                )
                            )
    
    return variations

def create_sv(variation_type, contig, pos, length, from_contig=None, 
                from_loc = -1, nr_duplications = 2):
    """
     Create a variation object from the given variation_type and parameters.
     
     Arguments:
         variation_type (str): The type of variation, e.g. "insertion".
         contig (str): The contig where the sv is located
         pos (int): The 0-based position of the base pair before
                   the sv.
         length (int): The length of the event.
         from_loc (int): The 0-based position from where the event sequence is
                         taken. Can be -1 then a random sequence is generated.
         from_contig (str): The contig where the sv originates. 
                            If None => same contig as sv is at.
         nr_duplications (int): If duplication this is the number of times that
                                 the sequence is duplicated.
         
    Returns:
        variations (list): A list ov variation objects
    
    """
    if variation_type == "insertion":
        return Insertion(contig, pos, length, from_contig, from_loc)
    
    elif variation_type == "transversion":
        return Transversion(contig, pos, length)
    
    elif variation_type == "deletion":
        return Deletion(contig, pos, length)
    
    elif variation_type == "duplication":
        return Duplication(contig, pos, length, nr_duplications)
    
    else:
        return None

def check_variations(variations, logger):
    """
    Sort the variations and check if any of the variations overlap.
    Variations can not overlap, if any variations ovelap the program
    exits.
    
    Arguments:
        variations (list): List of variation objects
        logger (logger): An logbook object to print messages.
    
    Returns:
        sorted_variations (list): A list of sorted variations
    """
    # We first check if any of the intervals are overlapping
    intervals = []
    start = 0
    stop = 0
    i = 0
    for variation in variations:
        i += 1
        interval_start = variation.pos
        interval_stop = interval_start + variation.length
        interval_id = '_'.join(['id', str(i)])
        intervals.append(
            [
                interval_start,
                interval_stop,
                interval_id
            ]
        )
        if interval_stop > stop:
            stop = interval_stop
    
    interval_tree = IntervalTree(intervals, start, stop)
    
    for variation in variations:
        # If features are overlapping we abbort
        interval_start = variation.pos
        interval_stop = interval_start + variation.length
        interval = [interval_start, interval_stop]
        if len(interval_tree.find_range(interval)) > 1:
            logger.critical("Interval at contig {0}, position {1} is overlapping "\
                    "another interval. Pleast check your variations file."\
                    "Variants are not allowed to overlap.\n Exiting.".format(
                        variation.contig,
                        variation.pos
                    ))
            sys.exit(1)
    
    return sorted(variations, key=lambda v: v.pos)

def write_donor_contigs(normal_contigs, variations, sorted_contigs, 
                        genome_name, outfile):
    """
    docstring for print_donor_contigs
    """
    
    for contig in sorted_contigs:
        start_position = 0
        contig_variants = variations[contig]
        normal_sequence = normal_contigs[contig]
        contig_length = len(normal_sequence)
        
        outfile.write(">{0}|dna:chromosome|chromosome:{1}:{0}:1:{2}:1|DONOR\n".format(
                                contig,
                                genome_name,
                                contig_length
                                )
                            )
        for variant in contig_variants:
            variant_start = variant.pos
            variant_end = variant.pos + variant.length
            #print the part of genome without svs
            outfile.write(normal_sequence[start_position:variant_start])
            #then print the sv sequence
            outfile.write(variant.get_sequence(normal_contigs))
            start_position = variant_end
        
        outfile.write(normal_sequence[start_position:contig_length] + '\n')



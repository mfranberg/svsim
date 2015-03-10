from collections import defaultdict
from svsim import create_sv

VARIATIONS = ("insertion", "deletion", "duplication", "translocation", 
                "transversion")


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

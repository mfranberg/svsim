from interval_tree import IntervalTree

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

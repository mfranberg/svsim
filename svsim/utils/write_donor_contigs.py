
def write_donor_contigs(normal_contigs, variations, sorted_contigs, 
                        genome_name, outfile):
    """
    Write the donor contigs to a file.
    
    Arguments:
        normal_contigs (pyfasta): A pyfasta object with the normal sequences
        variations (dict): A dictionary with contig names as keys and list of
                            variations as values
        sorted_contigs (list): A list with the contig names in original order
        genome_name (str): then genome name
        outfile (file_stream): The output file stream
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
            variant_end = variant.pos + variant.get_delta()
            # Print the part of genome without svs
            outfile.write(normal_sequence[start_position:variant_start])
            # Then print the sv sequence
            outfile.write(variant.get_sequence(normal_contigs))
            start_position = variant_end
        
        outfile.write(normal_sequence[start_position:contig_length] + '\n')

import pyfasta

def calculate_num_reads(coverage, read_length, genome_length):
    """
    Calculate the number of reads required to get a certain coverage.
    
    Arguments:
        coverage (int): Desired mean genome coverage.
        read_length (int): Length of each read in a pair.
        genome_length (int): Estimated genome length.
    
    Returns:
        Number of reads need to get the given coverage.
    """
    return int((genome_length * coverage) / (read_length * 2.0))

def get_genome_length(genome_path):
    """
    Opens the given fasta file and calculates the total length 
    of all contigs.
    
    Arguments:
        genome_path (str): Path to a fasta file.
        
    Returns:
        The length of the genome.
    """
    genome_fasta = pyfasta.Fasta(genome_path)
    return sum(len(genome_fasta[contig]) for contig in genome_fasta.iterkeys())


def get_genome_sequence(genome_path):
    genome_fasta = pyfasta.Fasta(genome_path)
    return [genome_fasta[contig] for contig in genome_fasta.iterkeys()]


def reverse_complement(string):
    """
        Reverse complements a DNA string
        Arguments:
        string  - A DNA string
        Returns:
            A python string object that represents the
            reverse complement of the input (DNA) string.
    """
    rev_nuc={'A':'T','C':'G','G':'C','T':'A','N':'N','X':'X', 
    'a':'t', 't':'a','g':'c','c':'g'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)
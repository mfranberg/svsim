import click
from svsim.variations import (Insertion, Transversion, Deletion, Duplication)

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


@click.command()
def cli():
    """docstring for cli"""
    small_insertion = create_sv('insertion', '1', 1, 1)
    print(small_insertion)
    

if __name__ == '__main__':
    cli()
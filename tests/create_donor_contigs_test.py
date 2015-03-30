
import logging

from svsim.vcf import *
from svsim.variations import *
from svsim import init_log

from svsim.utils import read_variations

from StringIO import StringIO


def test_read_insertion():
    """
    Tests that insertions are parsed properly.    
    """
    variation_file = StringIO()
    logger = logging.getLogger("svsim.create_donor_contigs.test")
    init_log(logger)
    contigs = ['contig1']
    
    variation_file.write("contig1 insertion 100 50\n")
    variation_file.seek(0)

    variations = read_variations(variation_file, contigs, logger)
    
    assert len(variations) == 1
    assert len(variations["contig1"]) == 1
    v = variations["contig1"][0]
    assert isinstance(v,variation.Insertion)
    assert v.pos == 100
    assert v.get_delta() == 0
    assert len(v.get_sequence("")) == 50
     
def test_read_deletion():
    """
    Tests that deletions are parsed properly.
    """
    variation_file = StringIO()
    logger = logging.getLogger("svsim.create_donor_contigs.test")
    init_log(logger)
    contigs = ['contig1']
    
    variation_file.write("contig1 deletion 100 50\n")
    variation_file.seek(0)
    
    variations = read_variations(variation_file, contigs, logger)
    
    assert len(variations) == 1
    assert len(variations["contig1"]) == 1
    v = variations["contig1"][0]
    assert isinstance(v, variation.Deletion)
    assert v.pos == 100
    assert v.get_delta() == 50
    assert len(v.get_sequence("")) == 0


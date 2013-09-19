from svsim.vcf import *
from svsim.variation import *

from scripts.create_donor_contigs import read_variations

from StringIO import StringIO

##
# Tests that insertions are parsed properly.
#
def test_read_insertion():
    variation_file = StringIO( )
    variation_file.write( "contig1 insertion 100 50\n" )
    variation_file.seek( 0 )

    variations = read_variations( variation_file )
    
    assert len( variations ) == 1
    assert len( variations[ "contig1" ] ) == 1
    v = variations[ "contig1" ][ 0 ]
    assert isinstance( v, variation.Insertion )
    assert v.pos == 100
    assert v.get_delta( ) == 0
    assert len( v.get_sequence( "" ) ) == 50
     
##
# Tests that deletions are parsed properly.
#
def test_read_deletion():
    variation_file = StringIO( )
    variation_file.write( "contig1 deletion 100 50\n" )
    variation_file.seek( 0 )

    variations = read_variations( variation_file )
    
    assert len( variations ) == 1
    assert len( variations[ "contig1" ] ) == 1
    v = variations[ "contig1" ][ 0 ]
    assert isinstance( v, variation.Deletion )
    assert v.pos == 100
    assert v.get_delta( ) == 50
    assert len( v.get_sequence( "" ) ) == 0


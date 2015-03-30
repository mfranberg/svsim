from svsim.vcf import *
from svsim.variations import *

from StringIO import StringIO

##
# Converts a line in the vcf file to a list of
# columns.
#
# @return A list of columns found in the line.
#
def get_fields(line):
    return line.strip( ).split( )

##
# Returns the info column as a dict, e.g. mapping "SVTYPE" to "2".
#
# @param fields The list of columns of one line in a vcf file.
#
# @return Return the info column as a dict.
#
def get_info(fields):
    return dict( tuple( item.split( "=" ) for item in fields[ 7 ].split( ";" ) ) )

##
# Tests that deletions are printed properly.
#
def test_vcf():
    string_file = StringIO( )
    vcf_file = VCFFile( string_file, False )

    genome = "123456789"
    deletion = Deletion( 1, 2 )
    vcf_file.write_variation( "test", genome, deletion )

    output = string_file.getvalue( )

    fields = get_fields( output )

    assert fields[ 0 ] == "test"
    assert fields[ 1 ] == "2"
    assert fields[ 3 ] == "234"
    assert fields[ 4 ] == "2"
    
    info = get_info( fields )
    assert info[ "SVTYPE" ] == "DEL"
    assert info[ "SVLEN" ] == "-2"
    assert info[ "END" ] == "4"

##
# Tests that insertions are printed properly.
#
def test_vcf2():
    string_file = StringIO( )
    vcf_file = VCFFile( string_file, False )

    genome = "123456789"
    insertion = Insertion( 1, 2, 7 )
    vcf_file.write_variation( "test", genome, insertion )

    output = string_file.getvalue( )

    fields = get_fields( output )

    assert fields[ 0 ] == "test"
    assert fields[ 1 ] == "2"
    assert fields[ 3 ] == "2"
    assert fields[ 4 ] == "289"
    
    info = get_info( fields )
    assert info[ "SVTYPE" ] == "INS"
    assert info[ "SVLEN" ] == "2"
    assert info[ "END" ] == "4"


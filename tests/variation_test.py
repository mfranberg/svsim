from svsim.variations import *


def test_null():
    """
    Tests basic properties of the SV base class object.
    """
    null = StructuralVariation('1', 2, 2)
    assert isinstance(null,variation.StructuralVariation)
    
    assert null.get_sequence({'1':"ABCDEFGHIJ"}) == "CD"
    assert null.get_delta() == 2

def test_insertion():
    """
    Tests basic properties of the Insertion object.
    """
    # Create a random insertion on contig '1', position 2
    # with length 5
    insertion = Insertion('1', 2, 5)
    assert isinstance(insertion,variation.Insertion)
    
    assert insertion.contig == '1'
    assert insertion.pos == 2
    assert insertion.length == 5
    assert insertion.from_contig == '1'
    assert insertion.from_loc == -1

    assert len(insertion.get_sequence("")) == 5
    assert insertion.get_delta() == 0
    
    insertion = Insertion('2', 2, 2, '2', 4)
    assert insertion.get_sequence({'2':"ABCDEFGH"}) == "EF"
    
def test_deletion():
    """
    Tests basic properties of the Deletion object.
    """
    deletion = Deletion('1', 2, 5)
    assert isinstance(deletion,variation.Deletion)
    
    assert deletion.get_sequence("ABCDEFGH") == ""
    assert deletion.get_delta() == 5

def test_transversion():
    """
    Tests basic properties of the Transversion object.
    """
    transversion = Transversion('1', 2, 5)
    assert isinstance(transversion,variation.Transversion)
    
    assert transversion.get_sequence({'1':"ABCDEFGH"}) == "GFEDC"
    assert transversion.get_delta() == 5

#
# ##
# # Tests that the reference genome is chunked properly.
# #
# # Insertion
# # Reference: 123456789AB
# # Donor:     12AB3456789AB
# #
# # Deletion
# # Reference: 123456789AB
# # Donor:     1236789
# #
# def test_create_chunks():
#     genome = "123456789AB"
#     variations = [ Insertion( 1, 2, 9 ) ]
#     chunks = create_chunks( variations, len( genome ) )
#
#     assert len( chunks ) == 3
#
#     assert chunks[ 0 ].get_sequence( genome ) == "12"
#     assert chunks[ 1 ].get_sequence( genome ) == "AB"
#     assert chunks[ 2 ].get_sequence( genome ) == "3456789AB"
#
#     variations = [ Deletion( 2, 2 ) ]
#
#     chunks = create_chunks( variations, len( genome ) )
#
#     assert len( chunks ) == 3
#     assert chunks[ 0 ].get_sequence( genome ) == "123"
#     assert chunks[ 1 ].get_sequence( genome ) == []
#     assert chunks[ 2 ].get_sequence( genome ) == "6789AB"

from svsim import variation

# Header in the VCF file
VCF_HEADER = """#fileformat=VCFv4.1
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

##
# Opens a vcf file for writing.
#
# @param Path to the file to open.
#
def open_vcf_file(path):
    vcf_file = open( path, "w" )
    return VCFFile( vcf_file )

##
# Represents a line in the vcf file. The idea is to
# set the parameters and the convert it to a string.
#
class VCFLine:
    def __init__(self):
        ##
        # Chromosome or contig name.
        #
        self.chrom = "."

        ##
        # 1-based position.
        #
        self.pos = "."

        ##
        # Identifier.
        #
        self.id = "."

        ##
        # Reference sequence, always includes the base pair before.
        #
        self.ref = "."

        ##
        # Alternative sequence, always includes the base pair before.
        #
        self.alt = "."

        ##
        # Info string.
        #
        self.info = "."

    ##
    # Creates and sets an id based on the contig name, sv type
    # and a integer sv_number.
    #
    # @param contig_name Name of the contig.
    # @param sv_type Type of sv.
    # @param sv_number An integer.
    #
    def set_id(self, contig_name, sv_type, sv_number):
        self.id = "{0}_{1}_{2}".format( contig_name, sv_type, sv_number )

    ##
    # Sets the info string for a variation.
    #
    # @param sv_type Type of sv, INS or DEL.
    # @param sv_len Length of the sv.
    # @param end The end position of the sv.
    #
    def set_info(self, sv_type, sv_len, end):
        self.info = "SVTYPE={0};SVLEN={1};END={2}".format( sv_type, sv_len, end )

    ##
    # Converts the object to a string that can be written as a
    # valid vcf line.
    # 
    def __str__(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t.\t.\t{5}".format( self.chrom, self.pos, self.id, self.ref, self.alt, self.info )

##
# Represents a vcf file.
# 
class VCFFile:
    ##
    # Constructor.
    #
    # @param vcf_file A writeable file.
    #
    def __init__(self, vcf_file):
        ##
        # The underlying file object.
        #
        self.vcf_file = vcf_file
        
        ##
        # An integer representing the number of vcf lines written
        # so far.
        #
        self.sv_number = 0

        ##
        # Chromosome name, by default the name of the contig.
        #
        self.chrom = -1

        self.vcf_file.write( VCF_HEADER )

    ##
    # Sets the name of the current chromosome.
    #
    # @param chrom Name of the chromosome.
    #
    def set_chrom(self, chrom):
        self.chrom = chrom

    ##
    # Writes a list of variation objects to the vcf file.
    #
    # @param contig_name Name of the contig.
    # @param sequence The sequence of the contig.
    # @param variations List of variations.
    #
    def write(self, contig_name, sequence, variations):
        for sv in variations:
            self.write_variation( contig_name, sequence, sv )

    ##
    # Writes a single variation to the vcf file.
    #
    # @param contig_name Name of the contig.
    # @param sequence The sequence of the contig.
    # @param variation The variation to write.
    #
    def write_variation(self, contig_name, sequence, sv):
        chrom = self.chrom
        if chrom == -1:
            chrom = contig_name

        vcf_line = VCFLine( )
        if isinstance( sv, variation.Deletion ):
            vcf_line.chrom = chrom
            vcf_line.pos = sv.pos + 1
            vcf_line.set_id( contig_name, "del", self.sv_number )
            
            vcf_line.ref = sequence[ sv.pos : sv.pos + sv.length + 1 ]
            vcf_line.alt = sequence[ sv.pos ]

            vcf_line.set_info( "DEL", -sv.length, sv.pos + sv.length + 1 )
        elif isinstance( sv, variation.Insertion ):
            vcf_line.chrom = chrom
            vcf_line.pos = sv.pos + 1
            vcf_line.set_id( contig_name, "ins", self.sv_number )

            vcf_line.ref = sequence[ sv.pos ]
            vcf_line.alt = sequence[ sv.pos ] + sv.get_sequence( sequence )

            vcf_line.set_info( "INS", sv.length, sv.pos + sv.length + 1 )
        else:
            raise ValueError( "write_variation: Bad variation {0}.".format( type( sv ) ) )
        
        self.vcf_file.write( str( vcf_line ) )
        self.vcf_file.write( "\n" )
        self.sv_number += 1

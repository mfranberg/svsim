import argparse
import os

from create_donor_contigs import create_donor_contigs
from simulate_reads import get_simulator
from map_reads import map_reads

from svsim import vcf

if __name__ == '__main__':
    parser = argparse.ArgumentParser( description="Simulates contigs with indels, simulates reads and map them." )
    parser.add_argument( 'normal_contig_file', type=argparse.FileType('r'), help='Path to the reference genome.' )
    parser.add_argument( 'variation_file', type=argparse.FileType('r'), help='A file containg a list of variations to simulate.' )
    parser.add_argument( 'output_dir', type=str, help='The output directory, where final and transitory files will be stored.' )
    parser.add_argument( '-t', type=str, choices=[ "metasim", "dwgsim" ], help='The simulator to use.', required=True )
    parser.add_argument( '-m', type=float, help='Mean of the library.', default=500 )
    parser.add_argument( '-s', type=float, help='Standard deviation of the library.', default=50 )
    parser.add_argument( '-C', type=float, help='Coverage.', default=10.0 )
    parser.add_argument( '-d', dest="field_sep", type=str, help='The field delimiter for the fasta identifier, default is |.', default="|" ) 
    parser.add_argument( '-i', dest="field_index", type=int, help='The 0-based index of field separated by field_sep that contains the relevant contig name.', default=0 ) 
    parser.add_argument( '-v', dest="vcf_file", type=vcf.open_vcf_file, help='Output file of the structural variations in VCF format.' )
    parser.add_argument( '-c', dest="chrom", type=str, help='Sets chromosome that will be written in the vcf file, name of the contig by default.' )

    args = parser.parse_args( )

    if not os.path.exists( args.output_dir ):
        os.makedirs( args.output_dir )

    donor_contig_path = os.path.join( args.output_dir, "donor_contigs.fa" )
    normal_reads_path = os.path.join( args.output_dir, "reads_normal" )
    donor_reads_path = os.path.join( args.output_dir, "reads_donor" )
    normal_bwa_path = os.path.join( args.output_dir, "mapped_normal" )
    donor_bwa_path = os.path.join( args.output_dir, "mapped_donor" )

    # Create indel genome
    if args.vcf_file and args.chrom:
        args.vcf_file.set_chrom( args.chrom )

    with open( donor_contig_path, "w" ) as donor_contig_file:
        create_donor_contigs( args.normal_contig_file.name,
                              args.variation_file,
                              donor_contig_file,
                              field_index = args.field_index,
                              field_sep = args.field_sep,
                              vcf_file = args.vcf_file )

    # Simulate reads
    simulator = get_simulator( args.t )
    simulator.coverage = args.C
    simulator.mean = args.m
    simulator.std = args.s
    simulator.read_length = 100

    simulator.simulate( args.normal_contig_file.name, normal_reads_path )
    simulator.simulate( donor_contig_path, donor_reads_path )

    # Map reads
    map_reads( normal_reads_path + "_pe1.fa", normal_reads_path + "_pe2.fa",
            args.normal_contig_file.name, normal_bwa_path )
    map_reads( donor_reads_path + "_pe1.fa", donor_reads_path + "_pe2.fa",
            args.normal_contig_file.name, donor_bwa_path )

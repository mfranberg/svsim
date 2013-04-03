import argparse
import os

from create_indel_genome import read_genome, read_variations, write_donor_genome
from simulate_reads import get_simulator
from map_reads import map_reads

if __name__ == '__main__':
    parser = argparse.ArgumentParser( description="Simulates a genome with indels, simulates reads and map them." )
    parser.add_argument( 'genome_file', type=argparse.FileType('r'), help='Path to the reference genome.' )
    parser.add_argument( 'variation_file', type=argparse.FileType('r'), help='A file containg a list of variations to simulate.' )
    parser.add_argument( 'output_dir', type=str, help='The output directory, where final and transitory files will be stored.' )
    parser.add_argument( '-t', type=str, choices=[ "metasim", "dwgsim" ], help='The simulator to use.', required=True )
    parser.add_argument( '-m', type=str, help='Mean of the library.', default=500 )
    parser.add_argument( '-s', type=str, help='Standard deviation of the library.', default=50 )
    parser.add_argument( '-c', type=str, help='Coverage.', default=10.0 )

    args = parser.parse_args( )

    if not os.path.exists( args.output_dir ):
        os.makedirs( args.output_dir )

    indel_genome_path = os.path.join( args.output_dir, "indel_genome.fa" )
    normal_reads_path = os.path.join( args.output_dir, "reads_normal" )
    indel_reads_path = os.path.join( args.output_dir, "reads_indel" )
    normal_bwa_path = os.path.join( args.output_dir, "mapped_normal" )
    indel_bwa_path = os.path.join( args.output_dir, "mapped_indel" )

    # Create indel genome
    normal_genome = read_genome( args.genome_file )
    variations = read_variations( args.variation_file )
    write_donor_genome( normal_genome, variations, indel_genome_path )
    del normal_genome # Make sure to clear memory if this is big

    # Simulate reads
    simulator = get_simulator( args.t )
    simulator.coverage = args.c
    simulator.mean = args.m
    simulator.std = args.s
    simulator.read_length = 100

    simulator.simulate( args.genome_file.name, normal_reads_path )
    simulator.simulate( indel_genome_path, indel_reads_path )

    # Map reads
    map_reads( normal_reads_path + "_pe1.fa", normal_reads_path + "_pe2.fa",
            args.genome_file.name, normal_bwa_path )
    map_reads( indel_reads_path + "_pe1.fa", indel_reads_path + "_pe2.fa",
            args.genome_file.name, indel_bwa_path )

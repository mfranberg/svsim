import argparse
import logging
import os

from create_donor_contigs import create_donor_contigs, add_sv_arguments
from simulate_reads import get_simulator, add_simulator_arguments
from map_reads import map_reads

from svsim import vcf
import svsim.log as log

if __name__ == '__main__':
    parser = argparse.ArgumentParser( description="Simulates contigs with indels, simulates reads and map them." )
    parser.add_argument( 'normal_contig_file', type=argparse.FileType('r'), help='Path to the reference genome.' )
    parser.add_argument( 'variation_file', type=argparse.FileType('r'), help='A file containg a list of variations to simulate.' )
    parser.add_argument( 'output_dir', type=str, help='The output directory, where final and transitory files will be stored.' )
    parser.add_argument( '-l', type=str, help='Log file path' )
    
    sim_group = parser.add_argument_group( title = "Simulation parameters", description = "Parameters that control the read simulation." )
    add_simulator_arguments( sim_group )

    sv_group = parser.add_argument_group( title = "Donor contig generation", description = "Parameters that control how donor contigs are generated" )
    add_sv_arguments( sv_group )

    args = parser.parse_args( )

    # Set up logging
    log.init_log( args.l )

    if not os.path.exists( args.output_dir ):
        os.makedirs( args.output_dir )

    donor_contig_path = os.path.join( args.output_dir, "donor_contigs.fa" )
    normal_reads_path = os.path.join( args.output_dir, "reads_normal" )
    donor_reads_path = os.path.join( args.output_dir, "reads_donor" )
    normal_bwa_path = os.path.join( args.output_dir, "mapped_normal" )
    donor_bwa_path = os.path.join( args.output_dir, "mapped_donor" )

    # Create indel genome
    logging.info( "Pipeline: Creating donor genome." )
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
    logging.info( "Pipeline: Simulating reads." )
    simulator = get_simulator( args.t )
    simulator.coverage = args.C
    simulator.mean = args.m
    simulator.std = args.s
    simulator.read_length = 100
    simulator.read_error = args.r

    simulator.simulate( args.normal_contig_file.name, normal_reads_path )
    simulator.simulate( donor_contig_path, donor_reads_path )

    # Map reads
    logging.info( "Pipeline: Mapping reads." )
    map_reads( normal_reads_path + "_pe1.fa", normal_reads_path + "_pe2.fa",
            args.normal_contig_file.name, normal_bwa_path )
    map_reads( donor_reads_path + "_pe1.fa", donor_reads_path + "_pe2.fa",
            args.normal_contig_file.name, donor_bwa_path )

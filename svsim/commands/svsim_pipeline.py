import click
import logging
import os

# from svsim.commands import create_donor_contigs, add_sv_arguments
# from simulate_reads import get_simulator, add_simulator_arguments
# from map_reads import map_reads

from svsim import vcf
from svsim.commands import (create_donor_contigs, simulate_reads, map_reads)
import svsim.log as log


@click.command()
@click.argument('genome',
                    type=click.Path(exists=True)
)
@click.argument('variations',
                    type=click.File('r')
)
@click.argument('output_dir',
                    type=click.Path(),
)
@click.option('-l', '--log_file',
                    type=click.Path(exists=False),
                    help="Path to log file."
)
@click.option('-vcf', '--vcf_file',
                    type=click.Path(exists=False),
                    help="Path to vcf file."
)
@click.option('-c', '--chrom',
                    nargs=1,
                    type=str,
                    help='Sets chromosome that will be written in the vcf file, name of the contig by default.'
)
@click.pass_context
def pipeline(ctx, genome, variations, output_dir, log_file, vcf_file, chrom):
    """
    Run the svsim pipeline with default arguments.
    """
    
    # Set up logging
    log.init_log( log_file )
    if not os.path.exists( output_dir ):
        os.makedirs( output_dir )
    
    donor_contig_path = os.path.join( output_dir, "donor_contigs.fa" )
    normal_reads_path = os.path.join( output_dir, "reads_normal" )
    donor_reads_path = os.path.join( output_dir, "reads_donor" )
    normal_bwa_path = os.path.join( output_dir, "mapped_normal" )
    donor_bwa_path = os.path.join( output_dir, "mapped_donor" )
    # Create indel genome
    logging.info( "Pipeline: Creating donor genome." )
    with open( donor_contig_path, "w" ) as donor_contig_file:
        ctx.invoke(create_donor_contigs,
                    normal_contig_file = genome,
                    variation_file = variations,
                    output = donor_contig_file,
                    delimiter='|',
                    field_index=0,
                    chrom=chrom,
                    vcf_file=vcf_file
                )

    # Simulate reads
    logging.info( "Pipeline: Simulating reads." )

    ctx.invoke(simulate_reads,
                    genome_file = genome,
                    output_prefix=donor_reads_path,
                    coverage=10.0,
                    mean=550.0,
                    standard_deviation=50.0,
                    simulator='dwgsim',
                    read_error_rate=0.02
                    )

    # Map reads
    logging.info( "Pipeline: Mapping reads." )
    
    ctx.invoke(map_reads,
                 genome = genome,
                 output = donor_bwa_path,
                 pe1_path = donor_reads_path + "_pe1.fa", 
                 pe2_path = donor_reads_path + "_pe2.fa",
                 logfile = log_file
                 )


if __name__ == '__main__':
    pipeline()




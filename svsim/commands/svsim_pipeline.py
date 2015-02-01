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
@click.option('-d', '--delimiter',
                    default='|',
                    help='The fasta identifier separator, default is |.'
)
@click.option('-i', '--field_index',
                    nargs=1,
                    default=0,
                    help='The 0-based index of field separated by field_sep that contains the relevant contig name.'
)
@click.option('-c', '--coverage',
                    nargs=1,
                    default=10.0,
                    help="The medium coverage of the simulated reads. Default 10."
)
@click.option('-m', '--mean',
                    nargs=1,
                    default=550.0,
                    help="Mean insert size of the library distribution for the simulated reads. Default 550 bp."
)
@click.option('-s', '--standard_deviation',
                    nargs=1,
                    default=50.0,
                    help="Standard deviation of the library distribution for the simulated reads. Default 50 bp."
)
@click.option('-t', '--simulator',
                    nargs=1,
                    type=click.Choice(["metasim", "dwgsim"]),
                    default='dwgsim',
                    help="Type of simulator 'metasim' or 'dwgsim', default 'dwgsim'."
)
@click.option('-r', '--read_error_rate',
                    nargs=1,
                    default=0.02,
                    help="Probability of a read error (not used in metasim)."
)
@click.option('-l', '--log_file',
                    type=click.Path(exists=False),
                    help="Path to log file."
)
@click.option('-vcf', '--vcf_file',
                    type=click.Path(exists=False),
                    help="Path to vcf file."
)
@click.option('-chr', '--chrom',
                    nargs=1,
                    type=str,
                    help='Sets chromosome that will be written in the vcf file, name of the contig by default.'
)
@click.pass_context
def pipeline(ctx, genome, variations, output_dir, delimiter, field_index, 
            coverage, standard_deviation, mean, simulator, read_error_rate, 
            log_file, vcf_file, chrom):
    """
    Run the svsim pipeline.
    
    That is, a genome should already exist then this command will run:\n
    svsim create_donor_contigs\n
    svsim simulate_reads\n
    svsim map_reads\n
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
                    delimiter=delimiter,
                    field_index=field_index,
                    chrom=chrom,
                    vcf_file=vcf_file
                )

    # Simulate reads
    logging.info( "Pipeline: Simulating reads." )

    ctx.invoke(simulate_reads,
                    genome_file = genome,
                    output_prefix=donor_reads_path,
                    coverage=coverage,
                    mean=mean,
                    standard_deviation=standard_deviation,
                    simulator=simulator,
                    read_error_rate=read_error_rate
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




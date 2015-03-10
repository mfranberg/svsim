

from __future__ import (division)
import sys
import os
import tempfile
import click

import logging

from svsim.reads import (MetaSimSimulator, DwgsimSimulator)
from svsim.warnings import SimulatorNotFoundError
from svsim import init_log

def get_simulator(simulator, logger):
    """
    Returns a simulator that can be used to simulate reads.
    
    Arguments:
        simulator (str): A string representing the simulator to use,
                        either 'metasim' or 'dwgsim'.
    
    Returns:
        (simulator): A matching simulator, or raises ValueError if no 
                     simulator was found.
    
    """
    if simulator == "metasim":
        return MetaSimSimulator(logger)
    elif simulator == "dwgsim":
        return DwgsimSimulator(logger)
    else:
        raise SimulatorNotFoundError(simulator)

@click.command()
@click.argument('genome_file',
                    nargs=-1,
                    type=click.Path(exists=True),
)
@click.argument('output_prefix',
                    nargs=1,
                    type=click.Path(exists=False),
)
@click.option('-c', '--coverage',
                    nargs=1,
                    default=20.0,
                    help="The medium coverage."
)
@click.option('-m', '--mean',
                    nargs=1,
                    default=550.0,
                    help="Mean insert size of the library distribution."
)
@click.option('-s', '--standard_deviation',
                    nargs=1,
                    default=50.0,
                    help="Standard deviation of the library distribution."
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
@click.option('-l', '--logfile',
                    type=click.Path(exists=False),
                    help="Path to log file. If none logging is "\
                          "printed to stderr."
)
@click.option('-h', '--heterozygous',
                    is_flag=True,
                    help="Simulate heterozygous variations. This is done "\
                         "by simulating reads from both donor contigs and "\
                         "reference contigs."
)
@click.option('--loglevel',
                    type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 
                                        'CRITICAL']),
                    default='INFO',
                    help="Set the level of log output."
)
def simulate_reads(genome_file, output_prefix, coverage, mean, standard_deviation,
                    simulator, read_error_rate, logfile, heterozygous, loglevel):
    """
    Simulate reads from a given genome.
    
    If heterozygous is used please provide two genomes.
    """
    logger = logging.getLogger("svsim.simulate_reads")
    
    init_log(logger, logfile, loglevel)
    
    genome = genome_file[0]
    second_genome = None
    
    if heterozygous:
        if len(genome_file) != 2:
            logger.critical("Provide two genomes when heterozygous is used. "\
                            "Please see documentation")
            sys.exit(1)
        
        logger.info('Simulating reads from two genomes to get heterozygosity')
        second_genome = genome_file[1]
        # We will try to simulate half of the reads from each genome
        coverage = int(coverage/2)
    
    try:
        simulator = get_simulator(simulator, logger)
    except SimulatorNotFoundError as e:
        logger.critical("No such simulator found: {0}".format(e.name))
        sys.exit(1)
    
    logger.debug('Set simulator coverage to {0}'.format(coverage))
    simulator.coverage = coverage
    logger.debug('Set simulator mean to {0}'.format(mean))
    simulator.mean = mean
    logger.debug('Set simulator std deviation to {0}'.format(standard_deviation))
    simulator.std = standard_deviation
    logger.debug('Set simulator read length to {0}'.format(100))
    simulator.read_length = 100
    logger.debug('Set simulator read error rate to {0}'.format(read_error_rate))
    simulator.read_error = read_error_rate
    
    simulator.simulate(genome, output_prefix, second_genome)
    
    
    

if __name__ == '__main__':
    simulate_reads()

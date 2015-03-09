
import sys
import os

import click

from svsim.reads.metasim import MetaSimSimulator
from svsim.reads.dwgsim import DwgsimSimulator

def get_simulator(simulator):
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
        return MetaSimSimulator()
    elif simulator == "dwgsim":
        return DwgsimSimulator()
    else:
        raise ValueError("No such simulator found")

@click.command()
@click.argument('genome_file',
                    type=click.Path(exists=True),
)
@click.argument('output_prefix',
                    type=click.Path(exists=False),
)
@click.option('-c', '--coverage',
                    nargs=1,
                    default=10.0,
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
def simulate_reads(genome_file, output_prefix, coverage, mean, standard_deviation,
                    simulator, read_error_rate):
    """
    Simulate reads from a given genome.
    
    Adds arguments relevant to simulation to the parser.
    """
    simulator = get_simulator(simulator)
    simulator.coverage = coverage
    simulator.mean = mean
    simulator.std = standard_deviation
    simulator.read_length = 100
    simulator.read_error = read_error_rate
    
    simulator.simulate(genome_file, output_prefix)
    

if __name__ == '__main__':
    simulate_reads()

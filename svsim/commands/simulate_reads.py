import argparse
import os
import subprocess
import sys
import tempfile

from svsim.reads.metasim import MetaSimSimulator
from svsim.reads.dwgsim import DwgsimSimulator

##
# Returns a simulator that can be used to simulate reads.
#
# @param simulator A string representing the simulator to use,
#                  either 'metasim' or 'dwgsim'.
#
# @return A matching simulator, or raises ValueError if no simulator
#         was found.
#
def get_simulator(simulator):
    if simulator == "metasim":
        return MetaSimSimulator( )
    elif simulator == "dwgsim":
        return DwgsimSimulator( )
    else:
        raise ValueError( "No such simulator found" )

##
# Adds arguments relevant to simulation to the parser.
#
# @param group A parser group to add arguments to.
#
def add_simulator_arguments(group):
    group.add_argument( '-C', type=float, help='Coverage.', default=10.0 )
    group.add_argument( '-m', type=float, help='Mean of the library distribution.', default=550.0 )
    group.add_argument( '-s', type=float, help='Standard deviation of the library distribution.', default=50.0 )
    group.add_argument( '-t', type=str, choices=[ "metasim", "dwgsim" ], help="Type of simulator 'metasim' or 'dwgsim', default 'metasim'.", default="metasim" )
    group.add_argument( '-r', type=float, help='Probability of a read error (not used in metasim).', default=0.02 )

if __name__ == '__main__':
    parser = argparse.ArgumentParser( description="Simulates illumina reads with metasim." )
    parser.add_argument( 'genome_file', type=str, help='Path to the genome' )
    parser.add_argument( 'output_prefix', type=str, help='Output prefix for the paired end files, will apped _pe1.fa and pe2.fa to this.' )
    add_simulator_arguments( parser )

    args = parser.parse_args( )

    simulator = get_simulator( args.t )
    simulator.coverage = args.C
    simulator.mean = args.m
    simulator.std = args.s
    simulator.read_length = 100
    simulator.read_error = args.r

    simulator.simulate( args.genome_file, args.output_prefix )

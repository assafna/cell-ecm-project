import os

from libs.paths_lib import working_directory

SIMULATIONS = os.path.join(working_directory(), 'Simulations')

# Main
RAW = os.path.join(SIMULATIONS, 'Raw')
MANIPULATIONS = os.path.join(SIMULATIONS, 'Manipulations')
OUTPUTS = os.path.join(SIMULATIONS, 'Outputs')

# Manipulations
STRUCTURED = os.path.join(MANIPULATIONS, 'Structured')
FIBER_LENGTHS = os.path.join(MANIPULATIONS, 'Fibers Lengths')
NORMALIZATION = os.path.join(MANIPULATIONS, 'Normalization')
PLOTS = os.path.join(OUTPUTS, 'Plots')


def raw(_simulation):
    return os.path.join(RAW, _simulation)


def structured(_simulation):
    return os.path.join(STRUCTURED, _simulation)


def fiber_lengths(_simulation):
    return os.path.join(FIBER_LENGTHS, _simulation)


def normalization(_simulation):
    return os.path.join(NORMALIZATION, _simulation)

import os

from libs.paths_lib import working_directory

SIMULATIONS = os.path.join(working_directory(), 'Simulations')

# Main
RAW = os.path.join(SIMULATIONS, 'Raw')
MANIPULATIONS = os.path.join(SIMULATIONS, 'Manipulations')
OUTPUTS = os.path.join(SIMULATIONS, 'Outputs')

# Manipulations
STRUCTURED = os.path.join(MANIPULATIONS, 'Structured')
FIBERS_LENGTHS = os.path.join(MANIPULATIONS, 'Fibers Lengths')
FIBERS_DENSITIES = os.path.join(MANIPULATIONS, 'Fibers Densities')
NORMALIZATION = os.path.join(MANIPULATIONS, 'Normalization')
PLOTS = os.path.join(OUTPUTS, 'Plots')
PLOTS_HISTORY = os.path.join('Plots History')
RESULTS = os.path.join(OUTPUTS, 'Results')


def raw(_simulation):
    return os.path.join(RAW, _simulation)


def structured(_simulation):
    return os.path.join(STRUCTURED, _simulation)


def fibers_lengths(_simulation):
    return os.path.join(FIBERS_LENGTHS, _simulation)


def fibers_densities(_simulation):
    return os.path.join(FIBERS_DENSITIES, _simulation)


def normalization(_simulation):
    return os.path.join(NORMALIZATION, _simulation)


def plots(_simulation):
    return os.path.join(PLOTS, _simulation)


def plots_history(_simulation):
    return os.path.join(PLOTS_HISTORY, _simulation)


def results(_simulation):
    return os.path.join(RESULTS, _simulation)

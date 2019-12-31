import os

EXPERIMENTS = 'G:\\My Drive\\BGU\\Thesis\\Cell-ECM & Cell-ECM-Cell Project\\Data\\Experiments'
# EXPERIMENTS = os.path.join(working_directory(), 'Experiments')

# Main
RAW = os.path.join(EXPERIMENTS, 'Raw')
MANIPULATIONS = os.path.join(EXPERIMENTS, 'Manipulations')
OUTPUTS = os.path.join(EXPERIMENTS, 'Outputs')

# Outputs
OBJECTS = os.path.join(OUTPUTS, 'Objects')
OBJECTS_Z = os.path.join(OUTPUTS, 'Objects Z')
CELL_COORDINATES_TRACKED = os.path.join(OUTPUTS, 'Cell Coordinates Tracked')
CELL_COORDINATES_TRACKED_Z = os.path.join(OUTPUTS, 'Cell Coordinates Tracked Z')
PAIRS = os.path.join(OUTPUTS, 'Pairs')
FIBERS_DENSITY = os.path.join(OUTPUTS, 'Fibers Density')
NORMALIZATION = os.path.join(OUTPUTS, 'Normalization')
PLOTS = os.path.join(OUTPUTS, 'Plots')


def objects(_experiment, _is_z):
    return os.path.join(OBJECTS, _experiment) if not _is_z else os.path.join(OBJECTS_Z, _experiment)


def cell_coordinates_tracked(_experiment, _is_z):
    return os.path.join(CELL_COORDINATES_TRACKED, _experiment) if not _is_z else os.path.join(
        CELL_COORDINATES_TRACKED_Z, _experiment)


def pairs(_experiment):
    return os.path.join(PAIRS, _experiment)


def fibers_density(_experiment):
    return os.path.join(FIBERS_DENSITY, _experiment)


def normalization(_experiment):
    return os.path.join(NORMALIZATION, _experiment)


def plots(_experiment):
    return os.path.join(PLOTS, _experiment)

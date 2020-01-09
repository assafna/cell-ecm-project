import os

EXPERIMENTS = 'G:\\My Drive\\BGU\\Thesis\\Cell-ECM & Cell-ECM-Cell Project\\Data\\Experiments'
# EXPERIMENTS = os.path.join(working_directory(), 'Experiments')

# Main
RAW = os.path.join(EXPERIMENTS, 'Raw')
MANIPULATIONS = os.path.join(EXPERIMENTS, 'Manipulations')
OUTPUTS = os.path.join(EXPERIMENTS, 'Outputs')

# Outputs
INFORMATION = os.path.join(OUTPUTS, 'Information')
OBJECTS = os.path.join(OUTPUTS, 'Objects')
OBJECTS_Z = os.path.join(OUTPUTS, 'Objects Z')
CELL_COORDINATES_TRACKED = os.path.join(OUTPUTS, 'Cell Coordinates Tracked')
CELL_COORDINATES_TRACKED_Z = os.path.join(OUTPUTS, 'Cell Coordinates Tracked Z')
PAIRS = os.path.join(OUTPUTS, 'Pairs')
FIBERS_DENSITY = os.path.join(OUTPUTS, 'Fibers Density')
NORMALIZATION_LINES = os.path.join(OUTPUTS, 'Normalization Lines')
NORMALIZATION = os.path.join(OUTPUTS, 'Normalization')
PLOTS = os.path.join(OUTPUTS, 'Plots')


def get_series_text_file_name(_series):
    return 'series_' + str(_series.split()[1]) + '.txt'


def folders(_path):
    return [_folder for _folder in os.listdir(_path) if os.path.isdir(os.path.join(_path, _folder))]


def text_files(_path):
    return [_file for _file in os.listdir(_path) if _file.endswith('.txt')]


def information(_experiment, _series=None):
    _experiment_path = os.path.join(INFORMATION, _experiment)

    if _series is None:
        return _experiment_path

    return os.path.join(_experiment_path, _series)


def objects(_experiment, _series=None, _time_point=None):
    _experiment_path = os.path.join(OBJECTS, _experiment)

    if _series is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, _series)
    if _time_point is None:
        return _series_path

    return os.path.join(_series_path, _time_point)


def objects_z(_experiment, _series=None, _group=None, _time_point=None):
    _experiment_path = os.path.join(OBJECTS_Z, _experiment)

    if _series is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, _series)
    if _group is None:
        return _series_path

    _group_path = os.path.join(_series_path, _group)
    if _time_point is None:
        return _group_path

    return os.path.join(_group_path, _time_point)


def cell_coordinates_tracked(_experiment, _series=None):
    _experiment_path = os.path.join(CELL_COORDINATES_TRACKED, _experiment)

    if _series is None:
        return _experiment_path

    return os.path.join(_experiment_path, _series)


def cell_coordinates_tracked_z(_experiment, _series=None, _group=None):
    _experiment_path = os.path.join(CELL_COORDINATES_TRACKED_Z, _experiment)

    if _series is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, _series)
    if _group is None:
        return _series_path

    return os.path.join(_series_path, _group)


def pairs(_experiment, _series=None, _group=None):
    _experiment_path = os.path.join(PAIRS, _experiment)

    if _series is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, _series)
    if _group is None:
        return _series_path

    return os.path.join(_series_path, _group)


def fibers_density(_experiment, _series=None, _group=None, _z_group=None, _time_point=None):
    _experiment_path = os.path.join(FIBERS_DENSITY, _experiment)

    if _series is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, _series)
    if _group is None:
        return _series_path

    _group_path = os.path.join(_series_path, _group)
    if _z_group is None:
        return _group_path

    _z_group_path = os.path.join(_group_path, _z_group)
    if _time_point is None:
        return _z_group_path

    return os.path.join(_z_group_path, _time_point)


def normalization_lines(_experiment, _series=None):
    _experiment_path = os.path.join(NORMALIZATION_LINES, _experiment)

    if _series is None:
        return _experiment_path

    return os.path.join(_experiment_path, _series)


def normalization(_experiment, _series=None):
    _experiment_path = os.path.join(NORMALIZATION, _experiment)

    if _series is None:
        return _experiment_path

    return os.path.join(_experiment_path, get_series_text_file_name(_series))


def plots(_experiment, _series=None, _group=None, _z_group=None):
    _experiment_path = os.path.join(PLOTS, _experiment)

    if _series is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, _series)
    if _group is None:
        return _series_path

    _group_path = os.path.join(_series_path, _group)
    if _z_group is None:
        return _group_path

    return os.path.join(_group_path, _z_group)

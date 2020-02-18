import os

from libs.paths_lib import working_directory

# EXPERIMENTS = 'G:\\My Drive\\BGU\\Thesis\\Cell-ECM & Cell-ECM-Cell Project\\Data\\Experiments'
EXPERIMENTS = os.path.join(working_directory(), 'Experiments')

# Main
RAW = os.path.join(EXPERIMENTS, 'Raw')
MANIPULATIONS = os.path.join(EXPERIMENTS, 'Manipulations')
OUTPUTS = os.path.join(EXPERIMENTS, 'Outputs')

# Manipulations
STRUCTURED = os.path.join(MANIPULATIONS, 'Structured')
SERIESES = os.path.join(MANIPULATIONS, 'Serieses')
IMAGE_PROPERTIES = os.path.join(MANIPULATIONS, 'Image Properties')
OBJECTS = os.path.join(MANIPULATIONS, 'Objects')
CELL_COORDINATES_TRACKED = os.path.join(MANIPULATIONS, 'Cell Coordinates Tracked')
FIBERS_DENSITIES = os.path.join(MANIPULATIONS, 'Fibers Densities')
NORMALIZATION_LINES = os.path.join(MANIPULATIONS, 'Normalization Lines')
NORMALIZATION = os.path.join(MANIPULATIONS, 'Normalization')
BLACKLIST = os.path.join(MANIPULATIONS, 'Blacklist')

# Outputs
PAIRS = os.path.join(OUTPUTS, 'Pairs')
CELLS = os.path.join(OUTPUTS, 'Cells')
PLOTS = os.path.join(OUTPUTS, 'Plots')
IMAGES = os.path.join(OUTPUTS, 'Images')
ANIMATIONS = os.path.join(OUTPUTS, 'Animations')


def get_series_text_file_name(_series):
    return 'series_' + str(_series.split()[1]) + '.json'


def folders(_path):
    return [_folder for _folder in os.listdir(_path) if os.path.isdir(os.path.join(_path, _folder))]


def text_files(_path):
    return [_file for _file in os.listdir(_path) if _file.endswith('.txt')]


def image_files(_path):
    return [_file for _file in os.listdir(_path) if _file.endswith('.tif')]


def tif(_experiment, _series=None):
    _experiment_path = os.path.join(MANIPULATIONS, _experiment)

    if _series is None:
        return _experiment_path

    return os.path.join(_experiment_path, _series)


def structured(_experiment, _series=None, _group=None, _time_point=None):
    _experiment_path = os.path.join(STRUCTURED, _experiment)

    if _series is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, _series)
    if _group is None:
        return _series_path

    _group_path = os.path.join(_series_path, _group)
    if _time_point is None:
        return _group_path

    return os.path.join(_group_path, _time_point)


def serieses(_experiment, _series=None):
    _experiment_path = os.path.join(SERIESES, _experiment)

    if _series is None:
        return _experiment_path

    return os.path.join(_experiment_path, _series)


def image_properties(_experiment, _series=None):
    _experiment_path = os.path.join(IMAGE_PROPERTIES, _experiment)

    if _series is None:
        return _experiment_path

    return os.path.join(_experiment_path, _series)


def group_properties(_experiment, _series, _group):
    return os.path.join(structured(_experiment, _series, _group), 'properties.json')


def objects(_experiment, _series=None, _time_point=None):
    _experiment_path = os.path.join(OBJECTS, _experiment)

    if _series is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, _series)
    if _time_point is None:
        return _series_path

    return os.path.join(_series_path, _time_point)


def cell_coordinates_tracked(_experiment, _series=None):
    _experiment_path = os.path.join(CELL_COORDINATES_TRACKED, _experiment)

    if _series is None:
        return _experiment_path

    return os.path.join(_experiment_path, _series)


def pairs(_experiment, _series=None, _group=None):
    _experiment_path = os.path.join(PAIRS, _experiment)

    if _series is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, _series)
    if _group is None:
        return _series_path

    return os.path.join(_series_path, _group)


def fibers_densities(_experiment, _series=None, _group=None, _time_point=None):
    _experiment_path = os.path.join(FIBERS_DENSITIES, _experiment)

    if _series is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, _series)
    if _group is None:
        return _series_path

    _group_path = os.path.join(_series_path, _group)
    if _time_point is None:
        return _group_path

    return os.path.join(_group_path, _time_point)


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


def plots(_experiment, _series=None, _group=None):
    _experiment_path = os.path.join(PLOTS, _experiment)

    if _series is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, _series)
    if _group is None:
        return _series_path

    return os.path.join(_series_path, _group)


def images(_experiment, _series=None, _group=None):
    _experiment_path = os.path.join(IMAGES, _experiment)

    if _series is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, _series)
    if _group is None:
        return _series_path

    return os.path.join(_series_path, _group)


def blacklist(_experiment, _series_id=None):
    if _series_id is None:
        return os.path.join(BLACKLIST, _experiment + '.pkl')

    _experiment_path = os.path.join(BLACKLIST, _experiment)

    return os.path.join(_experiment_path, 'series_' + str(_series_id) + '.pkl')

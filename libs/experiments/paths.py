import os

from libs.paths_lib import working_directory

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
NORMALIZATION = os.path.join(MANIPULATIONS, 'Normalization')
BLACKLIST = os.path.join(MANIPULATIONS, 'Blacklist')

# Outputs
PLOTS = os.path.join(OUTPUTS, 'Plots')
IMAGES = os.path.join(OUTPUTS, 'Images')


def folders(_path):
    return [_folder for _folder in os.listdir(_path) if os.path.isdir(os.path.join(_path, _folder))]


def text_files(_path):
    return [_file for _file in os.listdir(_path) if _file.endswith('.txt')]


def image_files(_path):
    return [_file for _file in os.listdir(_path) if _file.endswith('.tif')]


def structured(_experiment, _series_id=None, _group=None, _time_frame=None):
    _experiment_path = os.path.join(STRUCTURED, _experiment)

    if _series_id is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, 'Series ' + str(_series_id))
    if _group is None:
        return _series_path

    _group_path = os.path.join(_series_path, _group)
    if _time_frame is None:
        return _group_path

    return os.path.join(_group_path, str(_time_frame) + '.pkl')


def serieses(_experiment, _series_id=None):
    _experiment_path = os.path.join(SERIESES, _experiment)

    if _series_id is None:
        return _experiment_path

    return os.path.join(_experiment_path, 'series_' + str(_series_id) + '_bc.tif')


def image_properties(_experiment, _series_id=None):
    _experiment_path = os.path.join(IMAGE_PROPERTIES, _experiment)

    if _series_id is None:
        return _experiment_path

    return os.path.join(_experiment_path, 'series_' + str(_series_id) + '.json')


def group_properties(_experiment, _series_id, _group):
    return os.path.join(structured(_experiment, _series_id, _group), 'properties.json')


def objects(_experiment, _series_id=None, _time_frame=None):
    _experiment_path = os.path.join(OBJECTS, _experiment)

    if _series_id is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, 'Series ' + str(_series_id))
    if _time_frame is None:
        return _series_path

    return os.path.join(_series_path, 'tp_' + str(_time_frame) + '.txt')


def cell_coordinates_tracked(_experiment, _series_id=None):
    _experiment_path = os.path.join(CELL_COORDINATES_TRACKED, _experiment)

    if _series_id is None:
        return _experiment_path

    return os.path.join(_experiment_path, 'series_' + str(_series_id) + '.txt')


def normalization(_experiment, _series_id=None):
    _experiment_path = os.path.join(NORMALIZATION, _experiment)

    if _series_id is None:
        return _experiment_path

    return os.path.join(_experiment_path, 'series_' + str(_series_id) + '.json')


def images(_experiment, _series_id=None, _group=None):
    _experiment_path = os.path.join(IMAGES, _experiment)

    if _series_id is None:
        return _experiment_path

    _series_path = os.path.join(_experiment_path, 'Series ' + str(_series_id))
    if _group is None:
        return _series_path

    return os.path.join(_series_path, _group)


def blacklist(_experiment, _series_id=None):
    if _series_id is None:
        return os.path.join(BLACKLIST, _experiment + '.pkl')

    _experiment_path = os.path.join(BLACKLIST, _experiment)

    return os.path.join(_experiment_path, 'series_' + str(_series_id) + '.pkl')

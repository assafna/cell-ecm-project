import os

from PIL import Image
from PIL.TiffTags import TAGS

from libs.experiments import paths, save, config


def process_series(_experiment, _series_id, _overwrite=False):
    _experiment_path = paths.image_properties(_experiment)
    _properties_path = os.path.join(_experiment_path, 'series_' + str(_series_id) + '.json')
    if not _overwrite and os.path.isfile(_properties_path):
        return

    print('Creating image properties for:', _experiment, 'Series ' + str(_series_id), sep='\t')
    _path = paths.serieses(_experiment, _series_id)
    with Image.open(_path) as _img:
        _meta_dict = {TAGS[_key]: _img.tag[_key] for _key in _img.tag}

        _properties = {
            'experiment': _experiment,
            'series': int(_series_id),
            'dimensions': {
                'width': _meta_dict['ImageWidth'][0],
                'height': _meta_dict['ImageLength'][0]
            },
            'resolutions': {
                'x': _meta_dict['XResolution'][0][1] / _meta_dict['XResolution'][0][0],
                'y': _meta_dict['YResolution'][0][1] / _meta_dict['YResolution'][0][0],
                'z': float(str(_meta_dict['ImageDescription'][0].split('spacing=')[1]).split()[0])
            },
            'slices': int(str(_meta_dict['ImageDescription'][0].split('slices=')[1]).split()[0]),
            'frames': int(str(_meta_dict['ImageDescription'][0].split('frames=')[1]).split()[0]),
            'frames_interval': float(str(_meta_dict['ImageDescription'][0].split('finterval=')[1]).split()[0])
            # TODO: add location position X, Y & Z
        }

    save.image_properties(_experiment, _series_id, _properties)


def process_experiment(_experiment, _overwrite=False):
    for _series in paths.image_files(paths.serieses(_experiment)):
        process_series(_experiment, int(_series.split('_')[1]), _overwrite)


def process_all_experiments(_overwrite=False):
    for _experiment in config.experiments():
        process_experiment(_experiment, _overwrite)


if __name__ == '__main__':
    process_all_experiments()

import os

from PIL import Image
from PIL.TiffTags import TAGS

from libs.experiments import paths, save, config


def process_series(_experiment, _series):
    _series_id = str(_series.split()[1])
    _path = os.path.join(paths.tif(_experiment, _series), 'series_' + _series_id + '_bc.tif')
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
            'frames': int(str(_meta_dict['ImageDescription'][0].split('frames=')[1]).split()[0])
        }

    save.image_properties(_experiment, _series_id, _properties)


def process_experiment(_experiment):
    for _series in paths.folders(paths.tif(_experiment)):
        process_series(_experiment, _series)


def process_all_experiments():
    for _experiment in config.experiments():
        process_experiment(_experiment)


if __name__ == '__main__':
    process_all_experiments()

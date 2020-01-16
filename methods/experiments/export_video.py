import numpy as np
from PIL import Image

from libs.experiments import load


def process_group(_experiment, _series_id, _group):
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    for _time_point in range(len(_group_properties['time_points'])):
        print(_time_point)
        _time_point_image = load.structured_image(_experiment, _series_id, _group, _time_point)
        _z_slice = _group_properties['time_points'][_time_point]['left_cell']['coordinates']['z']
        _z_image = _time_point_image[_z_slice - 20:_z_slice + 20]
        _z_image_stretched = np.repeat(_z_image, repeats=int(round(_z_image.shape[2] / _z_image.shape[1])), axis=1)
        _average_across_z = np.rint(np.mean(_z_image_stretched, axis=0))
        _img = Image.fromarray(_average_across_z).convert('L')
        _img.save(str(_time_point) + '.png')


if __name__ == '__main__':
    process_group('SN41', 3, 'cells_0_1')

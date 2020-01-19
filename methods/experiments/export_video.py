import os
from multiprocessing.pool import Pool

import numpy as np
from PIL import Image

from libs.experiments import load, paths
from libs.experiments.config import CELL_DIAMETER_IN_MICRONS


def process_group(_experiment, _series_id, _group):
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    _group_path = paths.images(_experiment, 'Series ' + str(_series_id), _group)
    os.makedirs(_group_path, exist_ok=True)
    for _time_point in range(0, len(_group_properties['time_points']), 5):
        _time_point_image = load.structured_image(_experiment, _series_id, _group, _time_point)
        _z_slice = _group_properties['time_points'][_time_point]['left_cell']['coordinates']['z']
        _z_cell_diameter = CELL_DIAMETER_IN_MICRONS / _group_properties['time_points'][_time_point]['resolutions']['z']
        _z_image = _time_point_image[_z_slice - int(round(_z_cell_diameter / 2)):_z_slice + int(round(_z_cell_diameter / 2))]
        # _z_image_stretched = np.repeat(_z_image, repeats=int(round(_z_image.shape[2] / _z_image.shape[1])), axis=1)
        _average_across_z = np.rint(np.mean(_z_image, axis=0))
        _img = Image.fromarray(_average_across_z).convert('L')
        _img_path = os.path.join(_group_path, str(_time_point) + '.png')
        _img.save(_img_path)


def process_experiment(_experiment):
    _p = Pool()
    _answers = _p.starmap(process_group, load.experiment_groups_as_tuples(_experiment))
    _p.close()


if __name__ == '__main__':
    process_experiment('SN41')

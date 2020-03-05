import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import seaborn as sns
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH
from plotting import heatmap, save

EXPERIMENTS = ['SN16']
TIME_POINT = 18
EXPERIMENTS_STR = '_'.join(EXPERIMENTS)
REAL_CELLS = True
STATIC = False
BAND = True
DOMAIN_ABSOLUTE_VALUE = 5
VALUES_BY_CELL_DIAMETER = [
    round(_value, 10) for _value in
    np.arange(start=0, stop=DOMAIN_ABSOLUTE_VALUE * 2 + 0.1, step=0.1) - DOMAIN_ABSOLUTE_VALUE]
OFFSET_X = 0
CELLS_DISTANCES = [6, 7, 8, 9]
DIRECTION = 'inside'

# globals
_experiments = []
_experiments_fibers_densities = None
_z_array = None


def compute_data(_arguments):
    _offset_y_index, _offset_y, _offset_z_index, _offset_z = _arguments
    _fibers_densities_array = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        _normalization = load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))
        for _cell_id in ['left_cell', 'right_cell']:
            _fibers_density = \
                _experiments_fibers_densities[(_experiment, _series_id, _group, _offset_y, _offset_z, _cell_id)]

            # remove out of boundaries
            if _fibers_density[1]:
                continue

            _normalized_fibers_density = compute_lib.z_score(
                _x=_fibers_density[0],
                _average=_normalization['average'],
                _std=_normalization['std']
            )
            _fibers_densities_array.append(_normalized_fibers_density)

    if len(_fibers_densities_array) > 0:
        return _offset_z_index, _offset_y_index, np.mean(_fibers_densities_array)
    else:
        return _offset_z_index, _offset_y_index, None


def main():
    global _experiments, _experiments_fibers_densities, _z_array

    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_time_points_amount(_experiments, TIME_POINT)
    _experiments = filtering.by_distances(_experiments, CELLS_DISTANCES)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    if BAND:
        _experiments = filtering.by_band(_experiments)

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _offset_y, _offset_z, _cell_id in \
                product(VALUES_BY_CELL_DIAMETER, VALUES_BY_CELL_DIAMETER, ['left_cell', 'right_cell']):
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': ROI_LENGTH,
                'length_y': ROI_HEIGHT,
                'length_z': ROI_WIDTH,
                'offset_x': OFFSET_X,
                'offset_y': _offset_y,
                'offset_z': _offset_z,
                'cell_id': _cell_id,
                'direction': DIRECTION,
                'time_point': TIME_POINT - 1
            })

    _rois_dictionary, _rois_to_compute =\
        compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'offset_y', 'offset_z', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {}
    for _key in tqdm(_rois_dictionary, desc='Organizing Fibers Densities'):
        _experiments_fibers_densities[_key] = _fibers_densities[_rois_dictionary[_key][0]]

    # clean
    _fibers_densities = None
    _rois_dictionary = None
    _rois_to_compute = None

    _arguments = []
    for (_offset_y_index, _offset_y), (_offset_z_index, _offset_z) in \
            product(enumerate(VALUES_BY_CELL_DIAMETER), enumerate(VALUES_BY_CELL_DIAMETER)):
        _arguments.append(
            (_offset_y_index, _offset_y, _offset_z_index, _offset_z))

    _z_array = np.zeros(shape=(len(VALUES_BY_CELL_DIAMETER), len(VALUES_BY_CELL_DIAMETER)))
    with Pool(1) as _p:
        for _answer in tqdm(_p.imap_unordered(compute_data, _arguments), total=len(_arguments),
                            desc='Computing Heatmap'):
            _offset_z_index, _offset_y_index, _mean = _answer
            _z_array[_offset_z_index, _offset_y_index] = _mean
        _p.close()
        _p.join()

    # plot
    _fig = heatmap.create_plot(
        _x_labels=VALUES_BY_CELL_DIAMETER,
        _y_labels=VALUES_BY_CELL_DIAMETER,
        _z_array=_z_array,
        _x_axis_title='Offset in Z axis',
        _y_axis_title='Offset in XY axis',
        _color_scale=sns.color_palette('BrBG').as_hex(),
        _zmin=-6,
        _zmax=6
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_' + EXPERIMENTS_STR + '_real_' + str(REAL_CELLS) + '_static_' + str(STATIC) + '_band_' +
                  str(BAND) + '_heatmap'
    )


if __name__ == '__main__':
    main()

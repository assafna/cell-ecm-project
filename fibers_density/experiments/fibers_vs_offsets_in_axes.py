import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
import seaborn as sns
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH
from plotting import save

# based on time resolution
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41', 'SN44', 'SN45']
}
TIME_POINT = 18
REAL_CELLS = True
STATIC = False
BAND = True
OFFSET_X = 0
OFFSET_Y_START = -1.1
OFFSET_Y_END = 2.6
OFFSET_Y_STEP = 0.1
OFFSET_Z_START = -5
OFFSET_Z_END = 5
OFFSET_Z_STEP = 0.1
CELLS_DISTANCE_RANGE = [4, 10]
DIRECTION = 'inside'

# globals
_experiments = []
_experiments_fibers_densities = {}
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
        return _offset_y_index, _offset_z_index, np.mean(_fibers_densities_array)
    else:
        return _offset_y_index, _offset_z_index, None


def compute_z_array(_high_time_resolution):
    global _experiments, _experiments_fibers_densities, _z_array

    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_time_points_amount(_experiments, TIME_POINT)
    _experiments = filtering.by_distance_range(_experiments, CELLS_DISTANCE_RANGE)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    if BAND:
        _experiments = filtering.by_band(_experiments)
    print('Total experiments:', len(_experiments))

    _offsets_y = np.arange(start=OFFSET_Y_START, stop=OFFSET_Y_END + OFFSET_Y_STEP, step=OFFSET_Y_STEP)
    _offsets_z = np.arange(start=OFFSET_Z_START, stop=OFFSET_Z_END + OFFSET_Z_STEP, step=OFFSET_Z_STEP)
    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _offset_y, _offset_z, _cell_id in product(_offsets_y, _offsets_z, ['left_cell', 'right_cell']):
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

    _rois_dictionary, _rois_to_compute = \
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
            product(enumerate(_offsets_y), enumerate(_offsets_z)):
        _arguments.append((_offset_y_index, _offset_y, _offset_z_index, _offset_z))

    _z_array = np.zeros(shape=(len(_offsets_y), len(_offsets_z)))
    with Pool(CPUS_TO_USE) as _p:
        for _answer in tqdm(_p.imap_unordered(compute_data, _arguments), total=len(_arguments),
                            desc='Computing Heatmap'):
            _offset_y_index, _offset_z_index, _mean = _answer
            _z_array[_offset_y_index, _offset_z_index] = _mean
        _p.close()
        _p.join()

    return _z_array


def main(_high_time_resolution=False):
    global _experiments, _experiments_fibers_densities, _z_array

    compute_z_array(_high_time_resolution)

    # plot
    _offsets_y = np.arange(start=OFFSET_Y_START, stop=OFFSET_Y_END + OFFSET_Y_STEP, step=OFFSET_Y_STEP)
    _offsets_z = np.arange(start=OFFSET_Z_START, stop=OFFSET_Z_END + OFFSET_Z_STEP, step=OFFSET_Z_STEP)
    _colors_array = ['white', '#ea8500']
    _fig = go.Figure(
        data=go.Heatmap(
            x=_offsets_z,
            y=_offsets_y,
            z=_z_array,
            colorscale=sns.color_palette(_colors_array).as_hex(),
            colorbar={
                'tickmode': 'array',
                'tickvals': [-3, 2.5, 8],
                'ticktext': ['-3', 'Z-score', '8'],
                'tickangle': -90
            },
            showscale=True,
            zmin=-3,
            zmax=8
        ),
        layout={
            'xaxis': {
                'title': 'Offset in XY axis (cell diameter)',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [-4, -2, 0, 2, 4]
            },
            'yaxis': {
                'title': 'Offset in Z axis (cell diameter)',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [-1, 0, 1, 2]
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_high_time_' + str(_high_time_resolution) + '_real_' + str(REAL_CELLS) + '_static_' + str(STATIC) + '_band_' +
                  str(BAND)
    )


if __name__ == '__main__':
    main()

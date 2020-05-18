import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
import seaborn as sns
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import load, filtering, compute, paths, organize
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import save

# based on time resolution
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41', 'SN44', 'SN45']
}
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
DERIVATIVE = 1
CELLS_DISTANCE_RANGE = [4, 10]
DIRECTION = 'inside'
MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50
}

# globals
_experiments = []
_tuples_by_experiment = {}
_experiments_fibers_densities = {}
_z_array = None
_annotations_array = None


def compute_data(_arguments):
    _offset_y_index, _offset_y, _offset_z_index, _offset_z = _arguments
    _same_correlations_array = []
    _different_correlations_array = []
    for _experiment in _tuples_by_experiment:
        _experiment_tuples = _tuples_by_experiment[_experiment]

        for _same_index in range(len(_experiment_tuples)):
            _same_tuple = _experiment_tuples[_same_index]
            _same_experiment, _same_series, _same_group = _same_tuple

            _same_left_cell_fibers_densities = _experiments_fibers_densities[(
                _same_experiment,
                _same_series,
                _same_group,
                _offset_y,
                _offset_z,
                'left_cell'
            )]
            _same_right_cell_fibers_densities = _experiments_fibers_densities[(
                _same_experiment,
                _same_series,
                _same_group,
                _offset_y,
                _offset_z,
                'right_cell'
            )]

            _same_properties = \
                load.group_properties(_same_experiment, _same_series, _same_group)
            _same_left_cell_fibers_densities = compute.remove_blacklist(
                _same_experiment, _same_series, _same_properties['cells_ids']['left_cell'],
                _same_left_cell_fibers_densities)
            _same_right_cell_fibers_densities = compute.remove_blacklist(
                _same_experiment, _same_series, _same_properties['cells_ids']['right_cell'],
                _same_right_cell_fibers_densities)

            _same_left_cell_fibers_densities_filtered, _same_right_cell_fibers_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _same_left_cell_fibers_densities, _same_right_cell_fibers_densities
                )

            # ignore small arrays
            if len(_same_left_cell_fibers_densities_filtered) < \
                    MINIMUM_CORRELATION_TIME_POINTS[_same_experiment]:
                continue

            _same_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_same_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
            )
            for _different_index in range(len(_experiment_tuples)):
                if _same_index != _different_index:
                    _different_tuple = _experiment_tuples[_different_index]
                    _different_experiment, _different_series, _different_group = _different_tuple

                    for _same_cell_id, _different_cell_id in product(['left_cell', 'right_cell'],
                                                                     ['left_cell', 'right_cell']):
                        _same_fibers_densities = _experiments_fibers_densities[(
                            _same_experiment,
                            _same_series,
                            _same_group,
                            _offset_y,
                            _offset_z,
                            _same_cell_id
                        )]
                        _different_fibers_densities = _experiments_fibers_densities[(
                            _different_experiment,
                            _different_series,
                            _different_group,
                            _offset_y,
                            _offset_z,
                            _different_cell_id
                        )]

                        _different_properties = load.group_properties(
                            _different_experiment,
                            _different_series,
                            _different_group
                        )
                        _same_fibers_densities = compute.remove_blacklist(
                            _same_experiment,
                            _same_series,
                            _same_properties['cells_ids'][_same_cell_id],
                            _same_fibers_densities
                        )
                        _different_fibers_densities = compute.remove_blacklist(
                            _different_experiment,
                            _different_series,
                            _different_properties['cells_ids'][_different_cell_id],
                            _different_fibers_densities
                        )

                        _same_fibers_densities_filtered, _different_fibers_densities_filtered = \
                            compute.longest_same_indices_shared_in_borders_sub_array(
                                _same_fibers_densities, _different_fibers_densities
                            )

                        # ignore small arrays
                        if len(_same_fibers_densities_filtered) < \
                                MINIMUM_CORRELATION_TIME_POINTS[_different_experiment]:
                            continue

                        _different_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_fibers_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_different_fibers_densities_filtered, _n=DERIVATIVE)
                        )
                        _same_correlations_array.append(_same_correlation)
                        _different_correlations_array.append(_different_correlation)

    # compute fraction
    _annotation = None
    _same_minus_different = np.array(_same_correlations_array) - np.array(_different_correlations_array)
    _same_count = len(_same_minus_different[_same_minus_different > 0])
    if len(_same_minus_different) > 0:
        _same_fraction = round(_same_count / len(_same_minus_different), 10)
        _wilcoxon = wilcoxon(_same_minus_different)
        _p_value = _wilcoxon[1]
        if _p_value > 0.05:
            _annotation = {
                'text': 'x',
                'showarrow': False,
                'x': _offset_z,
                'y': _offset_y,
                'font': {
                    'size': 6,
                    'color': 'red'
                }
            }
    else:
        _same_fraction = None

    return _offset_y_index, _offset_z_index, _same_fraction, _annotation


def compute_z_array(_high_time_resolution):
    global _experiments, _tuples_by_experiment, _experiments_fibers_densities, _z_array, _annotations_array

    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_distance_range(_experiments, CELLS_DISTANCE_RANGE)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    _experiments = filtering.by_band(_experiments, _band=BAND)
    print('Total experiments:', len(_experiments))

    _offsets_y = np.arange(start=OFFSET_Y_START, stop=OFFSET_Y_END + OFFSET_Y_STEP, step=OFFSET_Y_STEP)
    _offsets_z = np.arange(start=OFFSET_Z_START, stop=OFFSET_Z_END + OFFSET_Z_STEP, step=OFFSET_Z_STEP)
    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple

        # stop when windows are overlapping
        _properties = load.group_properties(_experiment, _series_id, _group)
        _latest_time_point = len(_properties['time_points'])
        for _time_point in range(len(_properties['time_points'])):
            _cells_distance = \
                compute.cells_distance_in_cell_size_time_point(_experiment, _series_id, _group, _time_point)
            if _cells_distance - 1 - OFFSET_X * 2 < ROI_LENGTH * 2:
                _latest_time_point = _time_point - 1
                break

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
                'time_points': _latest_time_point
            })

    _rois_dictionary, _rois_to_compute = \
        compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'offset_y', 'offset_z', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {}
    for _key in tqdm(_rois_dictionary, desc='Organizing Fibers Densities'):
        _experiments_fibers_densities[_key] = [_fibers_densities[_tuple] for _tuple in _rois_dictionary[_key]]

    _tuples_by_experiment = organize.by_experiment(_experiments)

    # clean
    _fibers_densities = None
    _rois_dictionary = None
    _rois_to_compute = None

    _arguments = []
    for (_offset_y_index, _offset_y), (_offset_z_index, _offset_z) in \
            product(enumerate(_offsets_y), enumerate(_offsets_z)):
        _arguments.append(
            (_offset_y_index, _offset_y, _offset_z_index, _offset_z))

    _z_array = np.zeros(shape=(len(_offsets_y), len(_offsets_z)))
    _annotations_array = []
    with Pool(CPUS_TO_USE) as _p:
        for _answer in tqdm(_p.imap_unordered(compute_data, _arguments), total=len(_arguments),
                            desc='Computing Heatmap'):
            _offset_y_index, _offset_z_index, _same_fraction, _annotation = _answer
            _z_array[_offset_y_index, _offset_z_index] = _same_fraction
            if _annotation is not None:
                _annotations_array.append(_annotation)
        _p.close()
        _p.join()

    return _z_array


def main(_high_time_resolution=False):
    global _experiments, _tuples_by_experiment, _experiments_fibers_densities, _z_array, _annotations_array

    compute_z_array(_high_time_resolution)

    # plot
    _offsets_y = np.arange(start=OFFSET_Y_START, stop=OFFSET_Y_END + OFFSET_Y_STEP, step=OFFSET_Y_STEP)
    _offsets_z = np.arange(start=OFFSET_Z_START, stop=OFFSET_Z_END + OFFSET_Z_STEP, step=OFFSET_Z_STEP)
    _colors_array = ['black', 'white', '#ea8500']
    _fig = go.Figure(
        data=go.Heatmap(
            x=_offsets_z,
            y=_offsets_y,
            z=_z_array,
            colorscale=sns.color_palette(_colors_array).as_hex(),
            colorbar={
                'tickmode': 'array',
                'tickvals': [0, 0.5, 1],
                'ticktext': ['0', 'Fraction', '1'],
                'tickangle': -90
            },
            showscale=True,
            zmin=0,
            zmax=1
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
            },
            'annotations': _annotations_array
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_high_time_' + str(_high_time_resolution) + '_real_' + str(REAL_CELLS) + '_static_' +
                  str(STATIC) + '_band_' + str(BAND)
    )


if __name__ == '__main__':
    main()

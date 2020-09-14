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
from libs.experiments import load, filtering, compute, paths, organize, config
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments, \
    DERIVATIVE
from plotting import save

OFFSET_X = 0
OFFSET_Y_START = -1.1
OFFSET_Y_END = 2.6
OFFSET_Y_STEP = 0.1
OFFSET_Z_START = -5
OFFSET_Z_END = 5
OFFSET_Z_STEP = 0.1

PAIR_DISTANCE_RANGE = [4, 10]

# globals
_tuples = []
_tuples_by_experiment = {}
_experiments_fiber_densities = {}
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

            _same_left_cell_fiber_densities = _experiments_fiber_densities[(
                _same_experiment,
                _same_series,
                _same_group,
                _offset_y,
                _offset_z,
                'left_cell'
            )]
            _same_right_cell_fiber_densities = _experiments_fiber_densities[(
                _same_experiment,
                _same_series,
                _same_group,
                _offset_y,
                _offset_z,
                'right_cell'
            )]

            _same_properties = \
                load.group_properties(_same_experiment, _same_series, _same_group)
            _same_left_cell_fiber_densities = compute.remove_blacklist(
                _same_experiment, _same_series, _same_properties['cells_ids']['left_cell'],
                _same_left_cell_fiber_densities)
            _same_right_cell_fiber_densities = compute.remove_blacklist(
                _same_experiment, _same_series, _same_properties['cells_ids']['right_cell'],
                _same_right_cell_fiber_densities)

            _same_left_cell_fiber_densities_filtered, _same_right_cell_fiber_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _same_left_cell_fiber_densities, _same_right_cell_fiber_densities
                )

            # ignore small arrays
            if len(_same_left_cell_fiber_densities_filtered) < compute.minimum_time_frames_for_correlation(_same_experiment):
                continue

            _same_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_same_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
            )
            for _different_index in range(len(_experiment_tuples)):
                if _same_index != _different_index:
                    _different_tuple = _experiment_tuples[_different_index]
                    _different_experiment, _different_series, _different_group = _different_tuple

                    for _same_cell_id, _different_cell_id in product(['left_cell', 'right_cell'],
                                                                     ['left_cell', 'right_cell']):
                        _same_fiber_densities = _experiments_fiber_densities[(
                            _same_experiment,
                            _same_series,
                            _same_group,
                            _offset_y,
                            _offset_z,
                            _same_cell_id
                        )]
                        _different_fiber_densities = _experiments_fiber_densities[(
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
                        _same_fiber_densities = compute.remove_blacklist(
                            _same_experiment,
                            _same_series,
                            _same_properties['cells_ids'][_same_cell_id],
                            _same_fiber_densities
                        )
                        _different_fiber_densities = compute.remove_blacklist(
                            _different_experiment,
                            _different_series,
                            _different_properties['cells_ids'][_different_cell_id],
                            _different_fiber_densities
                        )

                        _same_fiber_densities_filtered, _different_fiber_densities_filtered = \
                            compute.longest_same_indices_shared_in_borders_sub_array(
                                _same_fiber_densities, _different_fiber_densities
                            )

                        # ignore small arrays
                        if len(_same_fiber_densities_filtered) < compute.minimum_time_frames_for_correlation(_different_experiment):
                            continue

                        _different_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_fiber_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_different_fiber_densities_filtered, _n=DERIVATIVE)
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


def compute_z_array(_band=True, _high_temporal_resolution=False, _offset_x=OFFSET_X, _offset_y_start=OFFSET_Y_START,
                    _offset_y_end=OFFSET_Y_END, _offset_y_step=OFFSET_Y_STEP, _offset_z_start=OFFSET_Z_START,
                    _offset_z_end=OFFSET_Z_END, _offset_z_step=OFFSET_Z_STEP):
    global _tuples, _tuples_by_experiment, _experiments_fiber_densities, _z_array, _annotations_array

    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=_high_temporal_resolution,
        _is_bleb=False,
        _is_bleb_from_start=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_pair_distance_range(_tuples, PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples, _band=_band)
    print('Total tuples:', len(_tuples))

    _offsets_y = np.arange(start=_offset_y_start, stop=_offset_y_end + _offset_y_step, step=_offset_y_step)
    _offsets_z = np.arange(start=_offset_z_start, stop=_offset_z_end + _offset_z_step, step=_offset_z_step)
    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        for _offset_y, _offset_z, _cell_id in product(_offsets_y, _offsets_z, ['left_cell', 'right_cell']):
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'length_z': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'offset_x': OFFSET_X,
                'offset_y': _offset_y,
                'offset_z': _offset_z,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, OFFSET_X)
            })

    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'offset_y', 'offset_z', 'cell_id'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _experiments_fiber_densities = {}
    for _key in tqdm(_windows_dictionary, desc='Organizing fiber Densities'):
        _experiments_fiber_densities[_key] = [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]

    _tuples_by_experiment = organize.by_experiment(_tuples)

    # clean
    _fiber_densities = None
    _windows_dictionary = None
    _windows_to_compute = None

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


def main(_band=True, _high_temporal_resolution=False, _offset_x=OFFSET_X, _offset_y_start=OFFSET_Y_START,
         _offset_y_end=OFFSET_Y_END, _offset_y_step=OFFSET_Y_STEP, _offset_z_start=OFFSET_Z_START,
         _offset_z_end=OFFSET_Z_END, _offset_z_step=OFFSET_Z_STEP):
    global _tuples, _tuples_by_experiment, _experiments_fiber_densities, _z_array, _annotations_array

    compute_z_array(_band=_band, _high_temporal_resolution=_high_temporal_resolution, _offset_x=OFFSET_X,
                    _offset_y_start=OFFSET_Y_START, _offset_y_end=OFFSET_Y_END, _offset_y_step=OFFSET_Y_STEP,
                    _offset_z_start=OFFSET_Z_START, _offset_z_end=OFFSET_Z_END, _offset_z_step=OFFSET_Z_STEP)

    # plot
    _offsets_y = np.arange(start=_offset_y_start, stop=_offset_y_end + _offset_y_step, step=_offset_y_step)
    _offsets_z = np.arange(start=_offset_z_start, stop=_offset_z_end + _offset_z_step, step=_offset_z_step)
    _colors_array = ['black', 'white', config.colors(1)]
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
        _filename='plot_high_time_' + str(_high_temporal_resolution) + '_band_' + str(_band)
    )


if __name__ == '__main__':
    main()

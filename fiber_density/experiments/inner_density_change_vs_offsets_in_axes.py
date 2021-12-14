import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
import seaborn as sns
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import load, filtering, compute, paths, config
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, all_experiments
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
_experiments_fiber_densities = {}
_z_array = None


def compute_data(_arguments):
    _offset_y_index, _offset_y, _offset_z_index, _offset_z = _arguments
    _fiber_density_changes_array = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _normalization = load.normalization_series_file_data(_experiment, _series_id)
        _properties = load.group_properties(_experiment, _series_id, _group)
        for _cell_id in ['left_cell', 'right_cell']:
            _cell_fiber_densities = \
                _experiments_fiber_densities[(_experiment, _series_id, _group, _offset_y, _offset_z, _cell_id)]

            _cell_fiber_densities = compute.remove_blacklist(
                _experiment,
                _series_id,
                _properties['cells_ids'][_cell_id],
                _cell_fiber_densities
            )

            _previous_cell_fiber_density_normalized = None
            _cell_values = []
            for _time_frame, _cell_fiber_density in enumerate(_cell_fiber_densities):

                # not out of border
                if _cell_fiber_density[1]:
                    _previous_cell_fiber_density_normalized = None
                    continue

                # normalize
                _cell_fiber_density_normalized = compute_lib.z_score(
                    _x=_cell_fiber_density[0],
                    _average=_normalization['average'],
                    _std=_normalization['std']
                )

                # no previous
                if _previous_cell_fiber_density_normalized is None:
                    _previous_cell_fiber_density_normalized = _cell_fiber_density_normalized
                    continue

                # change
                _cell_fiber_density_normalized_change = abs(_cell_fiber_density_normalized - _previous_cell_fiber_density_normalized)
                _previous_cell_fiber_density_normalized = _cell_fiber_density_normalized

                # save
                _cell_values.append(_cell_fiber_density_normalized_change)

            # save
            if len(_cell_values) > 0:
                _fiber_density_changes_array.append(np.mean(_cell_values))

    if len(_fiber_density_changes_array) > 0:
        return _offset_y_index, _offset_z_index, np.mean(_fiber_density_changes_array)
    else:
        return _offset_y_index, _offset_z_index, None


def compute_z_array(_band=True, _high_temporal_resolution=False, _offset_x=OFFSET_X, _offset_y_start=OFFSET_Y_START,
                    _offset_y_end=OFFSET_Y_END, _offset_y_step=OFFSET_Y_STEP, _offset_z_start=OFFSET_Z_START,
                    _offset_z_end=OFFSET_Z_END, _offset_z_step=OFFSET_Z_STEP):
    global _tuples, _experiments_fiber_densities, _z_array

    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=_high_temporal_resolution,
        _is_bleb=False,
        _is_dead_dead=False,
        _is_live_dead=False,
        _is_bead=False,
        _is_metastasis=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_time_frames_amount(_tuples, compute.minimum_time_frames_for_correlation(_experiments[0]))
    _tuples = filtering.by_pair_distance_range(_tuples, PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples, _band=_band)
    print('Total tuples:', len(_tuples))

    _offsets_y = np.arange(start=_offset_y_start, stop=_offset_y_end + _offset_y_step, step=_offset_y_step)
    _offsets_z = np.arange(start=_offset_z_start, stop=_offset_z_end + _offset_z_step, step=_offset_z_step)
    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _latest_time_frame = compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, OFFSET_X)
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
                'time_points': _latest_time_frame
            })

    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'offset_y', 'offset_z', 'cell_id'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    # clean
    _fiber_densities = None
    _windows_dictionary = None
    _windows_to_compute = None

    _arguments = []
    for (_offset_y_index, _offset_y), (_offset_z_index, _offset_z) in \
            product(enumerate(_offsets_y), enumerate(_offsets_z)):
        _arguments.append((_offset_y_index, _offset_y, _offset_z_index, _offset_z))

    _z_array = np.zeros(shape=(len(_offsets_y), len(_offsets_z)))
    with Pool(CPUS_TO_USE) as _p:
        for _answer in tqdm(_p.imap_unordered(compute_data, _arguments), total=len(_arguments),
                            desc='Computing heatmap'):
            _offset_y_index, _offset_z_index, _mean = _answer
            _z_array[_offset_y_index, _offset_z_index] = _mean
        _p.close()
        _p.join()

    return _z_array


def main(_band=True, _high_temporal_resolution=False, _offset_x=OFFSET_X, _offset_y_start=OFFSET_Y_START,
         _offset_y_end=OFFSET_Y_END, _offset_y_step=OFFSET_Y_STEP, _offset_z_start=OFFSET_Z_START,
         _offset_z_end=OFFSET_Z_END, _offset_z_step=OFFSET_Z_STEP):
    global _tuples, _experiments_fiber_densities, _z_array

    compute_z_array(_band=_band, _high_temporal_resolution=_high_temporal_resolution, _offset_x=OFFSET_X,
                    _offset_y_start=OFFSET_Y_START, _offset_y_end=OFFSET_Y_END, _offset_y_step=OFFSET_Y_STEP,
                    _offset_z_start=OFFSET_Z_START, _offset_z_end=OFFSET_Z_END, _offset_z_step=OFFSET_Z_STEP)

    # plot
    _offsets_y = np.arange(start=_offset_y_start, stop=_offset_y_end + _offset_y_step, step=_offset_y_step)
    _offsets_z = np.arange(start=_offset_z_start, stop=_offset_z_end + _offset_z_step, step=_offset_z_step)
    _colors_array = ['white', config.colors(1)]
    _fig = go.Figure(
        data=go.Heatmap(
            x=_offsets_z,
            y=_offsets_y,
            z=_z_array,
            colorscale=sns.color_palette(_colors_array).as_hex(),
            colorbar={
                'tickmode': 'array',
                'tickvals': [0, 0.35, 0.7],
                'ticktext': ['0.0', 'Z-score change', '0.7'],
                'tickangle': -90
            },
            showscale=True,
            zmin=0,
            zmax=0.7
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
        _filename='plot_high_time_' + str(_high_temporal_resolution) + '_band_' + str(_band)
    )


if __name__ == '__main__':
    main()

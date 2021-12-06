import os
from itertools import product

import numpy as np
import plotly.graph_objs as go
import seaborn as sns

from libs import compute_lib
from libs.experiments import filtering, load, compute, config, paths
from libs.experiments.config import all_experiments, QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, DERIVATIVE
from plotting import save

PAIR_DISTANCE_RANGE = [4, 10]

OFFSET_X = 0
OFFSET_Y = 0.5
OFFSET_Z = 0

BY = np.arange(start=0, stop=1 + 1 / 16, step=1 / 16)


def compute_data(_tuples, _arguments, _padding_y_by, _padding_z_by, _space_y_by, _space_z_by):
    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute, _subtract_border=True, _padding_y_by=_padding_y_by,
                                               _padding_z_by=_padding_z_by, _space_y_by=_space_y_by,
                                               _space_z_by=_space_z_by)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _correlations = []
    for _tuple in _tuples:
        _experiment, _series, _group = _tuple

        _left_cell_fiber_densities = _experiments_fiber_densities[(_experiment, _series, _group, 'left_cell')]
        _right_cell_fiber_densities = _experiments_fiber_densities[(_experiment, _series, _group, 'right_cell')]

        _properties = load.group_properties(_experiment, _series, _group)
        _left_cell_fiber_densities = compute.remove_blacklist(
            _experiment,
            _series,
            _properties['cells_ids']['left_cell'],
            _left_cell_fiber_densities
        )
        _right_cell_fiber_densities = compute.remove_blacklist(
            _experiment,
            _series,
            _properties['cells_ids']['right_cell'],
            _right_cell_fiber_densities
        )

        _left_cell_fiber_densities_filtered, _right_cell_fiber_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _left_cell_fiber_densities, _right_cell_fiber_densities
            )

        # ignore small arrays
        if len(_left_cell_fiber_densities_filtered) < compute.minimum_time_frames_for_correlation(_experiment):
            continue

        _correlation = compute_lib.correlation(
            compute_lib.derivative(_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
            compute_lib.derivative(_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
        )

        _correlations.append(_correlation)

    return np.mean(_correlations)


def main(_real_cells=True, _static=False, _band=True, _high_temporal_resolution=False):
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=_high_temporal_resolution,
        _is_bleb=False,
        _is_bleb_from_start=False,
        _is_dead_dead=False,
        _is_live_dead=False,
        _is_bead=False,
        _is_metastasis=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_pair_distance_range(_tuples, PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples, _real_pairs=_real_cells)
    _tuples = filtering.by_fake_static_pairs(_tuples, _fake_static_pairs=_static)
    _tuples = filtering.by_band(_tuples, _band=_band)
    print('Total tuples:', len(_tuples))

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _latest_time_frame = compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, OFFSET_X)
        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'length_z': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': _latest_time_frame
            })

    _z_array = np.zeros(shape=(len(BY), len(BY)))
    for (_padding_index, _padding_by), (_space_index, _space_by) in product(enumerate(BY), enumerate(BY)):
        print('Padding by: ', _padding_by, ', space by: ', _space_by)
        _correlation = compute_data(_tuples, _arguments, _padding_y_by=_padding_by, _padding_z_by=_padding_by,
                                    _space_y_by=_space_by, _space_z_by=_space_by)
        _z_array[_space_index, _padding_index] = _correlation

    # plot
    _colors_array = ['white', config.colors(1)]
    _fig = go.Figure(
        data=go.Heatmap(
            x=BY,
            y=BY,
            z=_z_array,
            colorscale=sns.color_palette(_colors_array).as_hex(),
            colorbar={
                'tickmode': 'array',
                'tickvals': [0, 0.5, 1],
                'ticktext': ['0', 'Correlation', '1'],
                'tickangle': -90
            },
            showscale=True,
            zmin=0,
            zmax=1
        ),
        layout={
            'xaxis': {
                'title': 'Border size (cell diameter)',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 0.5, 1, 1.5, 2]
            },
            'yaxis': {
                'title': 'Space from quantification window (cell diameter)',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 0.5, 1, 1.5, 2]
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

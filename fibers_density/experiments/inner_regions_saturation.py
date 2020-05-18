import os

import plotly.graph_objs as go
from scipy.stats import wilcoxon, pearsonr
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths, organize
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import save

# based on time resolution
EXPERIMENTS = ['SN16']
OFFSET_X = 0
OFFSET_Z = 0
DERIVATIVE = 1
CELLS_DISTANCE_RANGE = [4, 10]
BAND = True
REAL_CELLS = True
STATIC = False
MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50,
    'SN45': 50
}


def compute_fibers_densities(_offset_y):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_distance_range(_experiments, CELLS_DISTANCE_RANGE)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    _experiments = filtering.by_band(_experiments, _band=BAND)
    print('Total experiments:', len(_experiments))

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

        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': ROI_LENGTH,
                'length_y': ROI_HEIGHT,
                'length_z': ROI_WIDTH,
                'offset_x': OFFSET_X,
                'offset_y': _offset_y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': _latest_time_point
            })

    _rois_dictionary, _rois_to_compute = compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute, _saturation=True)

    _experiments_fibers_densities = {
        _key: [_fibers_densities[_tuple][:2] for _tuple in _rois_dictionary[_key]]
        for _key in _rois_dictionary
    }

    _experiments_saturation = {
        _key: [(_fibers_densities[_tuple][2], _fibers_densities[_tuple][1]) for _tuple in _rois_dictionary[_key]]
        for _key in _rois_dictionary
    }

    _tuples_by_experiment = organize.by_experiment(_experiments)

    _correlations_array = []
    _saturation_array = []
    for _experiment in _tuples_by_experiment:
        print('Experiment:', _experiment)
        _experiment_tuples = _tuples_by_experiment[_experiment]

        for _tuple in tqdm(_experiment_tuples, desc='Main loop'):
            _, _series, _group = _tuple

            _left_cell_fibers_densities = _experiments_fibers_densities[(_experiment, _series, _group, 'left_cell')]
            _right_cell_fibers_densities = _experiments_fibers_densities[(_experiment, _series, _group, 'right_cell')]

            _properties = load.group_properties(_experiment, _series, _group)
            _left_cell_fibers_densities = compute.remove_blacklist(
                _experiment,
                _series,
                _properties['cells_ids']['left_cell'],
                _left_cell_fibers_densities
            )
            _right_cell_fibers_densities = compute.remove_blacklist(
                _experiment,
                _series,
                _properties['cells_ids']['right_cell'],
                _right_cell_fibers_densities
            )

            _left_cell_fibers_densities_filtered, _right_cell_fibers_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _left_cell_fibers_densities, _right_cell_fibers_densities
                )

            # ignore small arrays
            if len(_left_cell_fibers_densities_filtered) < \
                    MINIMUM_CORRELATION_TIME_POINTS[_experiment]:
                continue

            _correlation = compute_lib.correlation(
                compute_lib.derivative(_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
            )

            # saturation
            _left_cell_saturation = _experiments_saturation[(_experiment, _series, _group, 'left_cell')]
            _right_cell_saturation = _experiments_saturation[(_experiment, _series, _group, 'right_cell')]

            _left_cell_saturation = compute.remove_blacklist(
                _experiment,
                _series,
                _properties['cells_ids']['left_cell'],
                _left_cell_saturation
            )
            _right_cell_saturation = compute.remove_blacklist(
                _experiment,
                _series,
                _properties['cells_ids']['right_cell'],
                _right_cell_saturation
            )

            _left_cell_saturation_filtered, _right_cell_saturation_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _left_cell_saturation, _right_cell_saturation
                )

            _left_saturation_fraction_last_time_point = _left_cell_saturation_filtered[-1]
            _right_saturation_fraction_last_time_point = _right_cell_saturation_filtered[-1]
            _group_saturation_fraction_mean = \
                (_left_saturation_fraction_last_time_point + _right_saturation_fraction_last_time_point) / 2

            _correlations_array.append(_correlation)
            _saturation_array.append(_group_saturation_fraction_mean)

    print('Total points:', len(_correlations_array))
    print('Wilcoxon of correlations around the zero:')
    print(wilcoxon(_correlations_array))
    print('Pearson correlation of correlations and saturation fraction mean:')
    print(pearsonr(_correlations_array, _saturation_array))

    return _correlations_array, _saturation_array


def main(_offset_y=0.5):
    _correlations_array, _saturation_array = compute_fibers_densities(_offset_y)

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_saturation_array,
            y=_correlations_array,
            mode='markers',
            marker={
                'size': 10,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'End mean saturation fraction',
                'zeroline': False,
                'range': [-0.01, 0.1],
                'tickmode': 'array',
                'tickvals': [0, 0.05, 0.1]
            },
            'yaxis': {
                'title': 'Inner regions correlation',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_offset_y_' + str(_offset_y)
    )


if __name__ == '__main__':
    main()

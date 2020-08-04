import os

import plotly.graph_objs as go
from scipy.stats import pearsonr
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import save

# based on time resolution
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41', 'SN44', 'SN45']
}
OFFSET_X = 0
OFFSET_Z = 0
DERIVATIVE = 1
OUT_OF_BOUNDARIES = False
MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 15,
    'SN41': 50,
    'SN44': 50,
    'SN45': 50
}


def main(_real_cells=True, _static=False, _band=None, _high_time_resolution=True, _offset_y=0.5):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_real_cells(_experiments, _real_cells=_real_cells)
    _experiments = filtering.by_static_cells(_experiments, _static=_static)
    _experiments = filtering.by_band(_experiments, _band=_band)
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
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {
        _key: [_fibers_densities[_tuple] for _tuple in _rois_dictionary[_key]]
        for _key in _rois_dictionary
    }

    _correlations = []
    _pair_distances = []
    for _tuple in tqdm(_experiments, desc='Experiments loop'):
        _experiment, _series_id, _group = _tuple
        _properties = load.group_properties(_experiment, _series_id, _group)

        _left_cell_fibers_densities = _experiments_fibers_densities[(_experiment, _series_id, _group, 'left_cell')]
        _right_cell_fibers_densities = _experiments_fibers_densities[(_experiment, _series_id, _group, 'right_cell')]

        _left_cell_fibers_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_cell_fibers_densities)
        _right_cell_fibers_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_cell_fibers_densities)

        if not OUT_OF_BOUNDARIES:
            _left_cell_fibers_densities, _right_cell_fibers_densities = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _left_cell_fibers_densities, _right_cell_fibers_densities
                )
        else:
            _left_cell_fibers_densities = [_fibers_density[0] for _fibers_density in _left_cell_fibers_densities]
            _right_cell_fibers_densities = [_fibers_density[0] for _fibers_density in _right_cell_fibers_densities]

        # ignore small arrays
        if len(_left_cell_fibers_densities) < MINIMUM_CORRELATION_TIME_POINTS[_experiment]:
            continue

        _correlations.append(compute_lib.correlation(
            compute_lib.derivative(_left_cell_fibers_densities, _n=DERIVATIVE),
            compute_lib.derivative(_right_cell_fibers_densities, _n=DERIVATIVE)
        ))
        _pair_distances.append(
            compute.cells_distance_in_cell_size_time_point(_experiment, _series_id, _group, _time_point=0)
        )

    print('Total pairs:', len(_correlations))
    print('Pearson:', pearsonr(_correlations, _pair_distances))

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_pair_distances,
            y=_correlations,
            mode='markers',
            marker={
                'size': 15,
                'color': 'black'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Pair distance (cell diameter)',
                'zeroline': False,
                # 'tickmode': 'array',
                # 'tickvals': [4, 6, 8, 10]
            },
            'yaxis': {
                'title': 'Correlation',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            # 'shapes': [
            #     {
            #         'type': 'line',
            #         'x0': -1,
            #         'y0': -1,
            #         'x1': -1,
            #         'y1': 1,
            #         'line': {
            #             'color': 'black',
            #             'width': 2
            #         }
            #     },
            #     {
            #         'type': 'line',
            #         'x0': -1,
            #         'y0': -1,
            #         'x1': 1,
            #         'y1': -1,
            #         'line': {
            #             'color': 'black',
            #             'width': 2
            #         }
            #     }
            # ]
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_real_' + str(_real_cells) + '_band_' + str(_band) + '_high_' + str(_high_time_resolution)
                  + '_y_' + str(_offset_y)
    )


if __name__ == '__main__':
    main()

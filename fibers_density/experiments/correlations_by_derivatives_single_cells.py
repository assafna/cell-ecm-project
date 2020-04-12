import os

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, organize, paths
from libs.experiments.config import SINGLE_CELL, ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
TIME_POINTS = 18
OUT_OF_BOUNDARIES = False
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['D', 'D\'', 'D\'\'']


def compute_single_cell_mean(_experiment, _series_id, _cell_tuples, _rois_dictionary, _fibers_densities):
    _cell_fibers_densities = []
    for _time_point in range(TIME_POINTS):
        _time_point_fibers_densities = []
        for _cell_tuple in _cell_tuples:
            _, _, _group = _cell_tuple
            for _direction in ['left', 'right']:
                _roi_tuple = _rois_dictionary[(_experiment, _series_id, _group, _direction)][_time_point]
                _fibers_density = _fibers_densities[_roi_tuple]

                if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                    continue

                _time_point_fibers_densities.append(_fibers_density)

        _cell_fibers_densities.append(np.mean(_time_point_fibers_densities))

    return _cell_fibers_densities


def main():
    _experiments = load.experiments_groups_as_tuples(SINGLE_CELL)
    _experiments = filtering.by_time_points_amount(_experiments, TIME_POINTS)
    _experiments = filtering.by_main_cell(_experiments)

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _direction in ['left', 'right']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': ROI_LENGTH,
                'length_y': ROI_HEIGHT,
                'length_z': ROI_WIDTH,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': 'cell',
                'direction': _direction,
                'time_points': TIME_POINTS
            })

    _rois_dictionary, _rois_to_compute = \
        compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'direction'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments = organize.by_single_cell_id(_experiments)
    print('Total experiments:', len(_experiments))
    _experiments_ids = list(_experiments.keys())

    _y_arrays = [[] for _i in DERIVATIVES]
    for _index_1 in tqdm(range(len(_experiments_ids)), desc='Main loop'):
        _tuple_1 = _experiments_ids[_index_1]
        _experiment_1, _series_id_1, _cell_id_1 = _tuple_1
        _fibers_densities_1 = compute_single_cell_mean(
            _experiment=_experiment_1,
            _series_id=_series_id_1,
            _cell_tuples=_experiments[_tuple_1],
            _rois_dictionary=_rois_dictionary,
            _fibers_densities=_fibers_densities
        )
        for _index_2 in range(_index_1 + 1, len(_experiments_ids)):
            _tuple_2 = _experiments_ids[_index_2]
            _experiment_2, _series_id_2, _cell_id_2 = _tuple_2
            _fibers_densities_2 = compute_single_cell_mean(
                _experiment=_experiment_2,
                _series_id=_series_id_2,
                _cell_tuples=_experiments[_tuple_2],
                _rois_dictionary=_rois_dictionary,
                _fibers_densities=_fibers_densities
            )
            for _derivative_index, _derivative in enumerate(DERIVATIVES):
                _y_arrays[_derivative_index].append(compute_lib.correlation(
                    compute_lib.derivative(_fibers_densities_1, _n=_derivative),
                    compute_lib.derivative(_fibers_densities_2, _n=_derivative)
                ))

    print('Total points:', len(_y_arrays[0]))
    print('Wilcoxon around the zero')
    for _y_array, _derivative in zip(_y_arrays, DERIVATIVES):
        print('Derivative:', _derivative, wilcoxon(_y_array))

    # plot
    _colors_array = ['#844b00', '#ea8500', '#edbc80']
    _fig = go.Figure(
        data=[
            go.Box(
                y=_y,
                name=_derivative,
                boxpoints='all',
                jitter=1,
                pointpos=0,
                line={
                    'width': 1
                },
                fillcolor='white',
                marker={
                    'size': 10,
                    'color': _color
                },
                opacity=0.7,
                showlegend=False
            ) for _y, _derivative, _color in zip(_y_arrays, DERIVATIVES_TEXT, _colors_array)
        ],
        layout={
            'xaxis': {
                'title': 'Density by derivatives',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Correlation',
                'range': [-1, 1],
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()

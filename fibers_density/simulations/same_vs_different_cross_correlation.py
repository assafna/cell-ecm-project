import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import compute, filtering, load, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import save

# by low connectivity
TIME_POINT = {
    False: 50,
    True: 35
}
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 2
CELLS_DISTANCE = 7
ALPHA = 1
BETA = 1
TIME_LAGS = [-2, -1, 0, 1, 2]


def compute_fibers_densities(_low_connectivity=False):
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, TIME_POINT[_low_connectivity])
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=True,
        _is_low_connectivity=_low_connectivity,
        _is_causality=True,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_causality(_simulations, _alpha=ALPHA, _beta=BETA)
    _simulations = filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)
    print('Total simulations:', len(_simulations))

    _arguments = []
    for _simulation in _simulations:
        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'simulation': _simulation,
                'length_x': ROI_WIDTH,
                'length_y': ROI_HEIGHT,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': TIME_POINT[_low_connectivity]
            })

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.roi_fibers_density_by_time, _arguments),
                total=len(_arguments), desc='Computing Rois & Fibers Densities'):
            _fibers_densities[
                (_keys['simulation'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    _same_correlation_vs_time_lag = {}
    _same_time_lags_arrays = [[] for _i in TIME_LAGS]
    _different_time_lags_arrays = [[] for _i in TIME_LAGS]
    _same_time_lags_highest = [0 for _i in TIME_LAGS]
    _different_time_lags_highest = [0 for _i in TIME_LAGS]
    for _same_index in tqdm(range(len(_simulations)), desc='Main loop'):
        _same_simulation = _simulations[_same_index]
        _same_left_cell_fibers_densities = _fibers_densities[(_same_simulation, 'left_cell')]
        _same_right_cell_fibers_densities = _fibers_densities[(_same_simulation, 'right_cell')]

        # time lag
        _same_highest_correlation = -1.1
        _same_highest_correlation_time_lag_index = 0
        _same_correlation_vs_time_lag[_same_simulation] = []
        for _time_lag_index, _time_lag in enumerate(TIME_LAGS):
            if _time_lag > 0:
                _same_left_cell_fibers_densities_time_lag = _same_left_cell_fibers_densities[:-_time_lag]
                _same_right_cell_fibers_densities_time_lag = _same_right_cell_fibers_densities[_time_lag:]
            elif _time_lag < 0:
                _same_left_cell_fibers_densities_time_lag = _same_left_cell_fibers_densities[-_time_lag:]
                _same_right_cell_fibers_densities_time_lag = _same_right_cell_fibers_densities[:_time_lag]
            else:
                _same_left_cell_fibers_densities_time_lag = _same_left_cell_fibers_densities
                _same_right_cell_fibers_densities_time_lag = _same_right_cell_fibers_densities

            _same_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fibers_densities_time_lag, _n=DERIVATIVE),
                compute_lib.derivative(_same_right_cell_fibers_densities_time_lag, _n=DERIVATIVE)
            )

            _same_time_lags_arrays[_time_lag_index].append(_same_correlation)
            _same_correlation_vs_time_lag[_same_simulation].append(_same_correlation)

            if _same_correlation > _same_highest_correlation:
                _same_highest_correlation = _same_correlation
                _same_highest_correlation_time_lag_index = _time_lag_index

        _same_time_lags_highest[_same_highest_correlation_time_lag_index] += 1

        for _different_index in range(len(_simulations)):
            if _same_index != _different_index:
                _different_simulation = _simulations[_different_index]
                for _same_cell_id, _different_cell_id in product(['left_cell', 'right_cell'],
                                                                 ['left_cell', 'right_cell']):
                    _same_fibers_densities = \
                        _fibers_densities[(_same_simulation, _same_cell_id)]
                    _different_fibers_densities = \
                        _fibers_densities[(_different_simulation, _different_cell_id)]

                    # time lag
                    _different_highest_correlation = -1.1
                    _different_highest_correlation_time_lag_index = 0
                    for _time_lag_index, _time_lag in enumerate(TIME_LAGS):
                        if _time_lag > 0:
                            _same_fibers_densities_time_lag = _same_fibers_densities[:-_time_lag]
                            _different_fibers_densities_time_lag = _different_fibers_densities[_time_lag:]
                        elif _time_lag < 0:
                            _same_fibers_densities_time_lag = _same_fibers_densities[-_time_lag:]
                            _different_fibers_densities_time_lag = _different_fibers_densities[:_time_lag]
                        else:
                            _same_fibers_densities_time_lag = _same_fibers_densities
                            _different_fibers_densities_time_lag = _different_fibers_densities

                        _different_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_fibers_densities_time_lag, _n=DERIVATIVE),
                            compute_lib.derivative(_different_fibers_densities_time_lag, _n=DERIVATIVE)
                        )

                        _different_time_lags_arrays[_time_lag_index].append(_different_correlation)

                        if _different_correlation > _different_highest_correlation:
                            _different_highest_correlation = _different_correlation
                            _different_highest_correlation_time_lag_index = _time_lag_index

                    _different_time_lags_highest[_different_highest_correlation_time_lag_index] += 1

    return _same_correlation_vs_time_lag, _same_time_lags_arrays, _different_time_lags_arrays, \
        _same_time_lags_highest, _different_time_lags_highest


def main(_low_connectivity=False):
    _same_correlation_vs_time_lag, _same_time_lags_arrays, _different_time_lags_arrays, _same_time_lags_highest, \
        _different_time_lags_highest = compute_fibers_densities(_low_connectivity)

    # individual plots
    for _same_simulation in _same_correlation_vs_time_lag:
        _fig = go.Figure(
            data=go.Scatter(
                x=TIME_LAGS,
                y=_same_correlation_vs_time_lag[_same_simulation],
                mode='markers',
                marker={
                    'size': 25,
                    'color': '#2e82bf'
                }
            ),
            layout={
                'xaxis': {
                    'title': 'Time lag',
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': TIME_LAGS
                },
                'yaxis': {
                    'title': 'Correlation',
                    'range': [-1, 1.1],
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': [-1, -0.5, 0, 0.5, 1]
                }
            }
        )

        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_' + _same_simulation
        )

    # box plots
    for _name, _arrays in zip(['same', 'different'], [_same_time_lags_arrays, _different_time_lags_arrays]):
        _fig = go.Figure(
            data=[
                go.Box(
                    y=_y,
                    name=_time_lag,
                    boxpoints=False,
                    line={
                        'width': 1
                    },
                    marker={
                        'size': 10,
                        'color': '#2e82bf'
                    },
                    showlegend=False
                ) for _y, _time_lag in zip(_arrays, TIME_LAGS)
            ],
            layout={
                'xaxis': {
                    'title': 'Time lag',
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': TIME_LAGS
                },
                'yaxis': {
                    'title': _name.capitalize() + ' network correlations',
                    'range': [-1, 1.1],
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': [-1, -0.5, 0, 0.5, 1]
                }
            }
        )

        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_box_low_con_' + str(_low_connectivity) + '_' + _name
        )

    # bar plot
    for _name, _sums in zip(['same', 'different'], [_same_time_lags_highest, _different_time_lags_highest]):
        _fig = go.Figure(
            data=go.Bar(
                x=TIME_LAGS,
                y=np.array(_sums) / sum(_sums),
                marker={
                    'color': '#2e82bf'
                }
            ),
            layout={
                'xaxis': {
                    'title': 'Time lag',
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': TIME_LAGS
                },
                'yaxis': {
                    'title': 'Highest correlations fraction',
                    'range': [0, 1.1],
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': [0, 0.5, 1]
                }
            }
        )

        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_bar_low_con_' + str(_low_connectivity) + '_' + _name
        )


if __name__ == '__main__':
    main()

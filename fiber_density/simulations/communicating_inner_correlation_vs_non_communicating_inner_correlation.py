import os
from itertools import product
from multiprocessing.pool import Pool

import plotly.graph_objs as go
from scipy.stats import ranksums
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import compute, filtering, load, paths, config
from libs.simulations.config import QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

TIME_POINT = 50
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 2
STDS = [0, 0.25, 0.5, 0.75]
PAIR_DISTANCE = 5


def compute_fiber_densities(_simulations):
    _arguments = []
    for _simulation in _simulations:
        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'simulation': _simulation,
                'length_x': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': TIME_POINT
            })

    _fiber_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.window_fiber_density_by_time, _arguments),
                total=len(_arguments), desc='Computing windows & fiber densities'):
            _fiber_densities[
                (_keys['simulation'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    return _fiber_densities


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, TIME_POINT)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=None,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False,
        _is_fibrin=False
    )
    _simulations = filtering.by_heterogeneities(_simulations, _stds=STDS)
    _simulations = filtering.by_pair_distance(_simulations, _distance=PAIR_DISTANCE)
    print('Total simulations:', len(_simulations))

    _fiber_densities = compute_fiber_densities(_simulations)

    # stds loop
    _std_communicating = [[] for _i in STDS]
    _std_non_communicating = [[] for _i in STDS]
    for _std_index, _std in enumerate(STDS):
        print('Std:', _std)
        _simulations_std = filtering.by_heterogeneity(_simulations, _std=_std)
        print('Total simulations:', len(_simulations_std))

        # communicating loop
        for _simulation in tqdm(_simulations_std, desc='Communicating pairs loop'):
            _left_cell_fiber_densities = _fiber_densities[(_simulation, 'left_cell')]
            _right_cell_fiber_densities = _fiber_densities[(_simulation, 'right_cell')]
            _correlation = compute_lib.correlation(
                compute_lib.derivative(_left_cell_fiber_densities, _n=DERIVATIVE),
                compute_lib.derivative(_right_cell_fiber_densities, _n=DERIVATIVE)
            )
            _std_communicating[_std_index].append(_correlation)

        # non-communicating loop
        _simulations_indices = range(len(_simulations_std))
        for _simulation_1_index in tqdm(_simulations_indices, desc='Non-communicating pairs loop'):
            _simulation_1 = _simulations_std[_simulation_1_index]
            for _simulation_2_index in _simulations_indices[_simulation_1_index + 1:]:
                _simulation_2 = _simulations_std[_simulation_2_index]
                for _simulation_1_cell_id, _simulation_2_cell_id in product(['left_cell', 'right_cell'],
                                                                            ['left_cell', 'right_cell']):
                    _simulation_1_fiber_densities = _fiber_densities[(_simulation_1, _simulation_1_cell_id)]
                    _simulation_2_fiber_densities = _fiber_densities[(_simulation_2, _simulation_2_cell_id)]
                    _correlation = compute_lib.correlation(
                        compute_lib.derivative(_simulation_1_fiber_densities, _n=DERIVATIVE),
                        compute_lib.derivative(_simulation_2_fiber_densities, _n=DERIVATIVE)
                    )
                    _std_non_communicating[_std_index].append(_correlation)

        print('Wilcoxon rank-sum test:', ranksums(_std_communicating[_std_index], _std_non_communicating[_std_index]))

    # histogram plot, once with and once without the legend
    _colors_array = config.colors(len(STDS))
    for _show_legend in [True, False]:
        _fig = go.Figure(
            data=[
                *[
                    go.Histogram(
                        x=_non_communicating,
                        name='STD ' + str(_std),
                        marker={
                            'color': _color
                        },
                        nbinsx=10,
                        histnorm='probability',
                        showlegend=_show_legend
                    ) for _std, _non_communicating, _color in zip(STDS, _std_non_communicating, _colors_array)
                ],
                *[
                    go.Scatter(
                        x=_communicating,
                        y=[0.5 - 0.02 * _std_index for _i in range(len(_communicating))],
                        mode='text',
                        text='●',
                        textfont={
                          'color': _color,
                          'size': 15
                        },
                        showlegend=False
                    ) for _std_index, _communicating, _color in zip(range(len(STDS)), _std_communicating, _colors_array)
                ]
            ],
            layout={
                'xaxis': {
                    'title': 'Correlation',
                    'zeroline': False
                },
                'yaxis': {
                    'title': 'Fraction',
                    'zeroline': False
                },
                'bargap': 0.1
            }
        )

        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_legend_' + str(_show_legend)
        )


if __name__ == '__main__':
    main()

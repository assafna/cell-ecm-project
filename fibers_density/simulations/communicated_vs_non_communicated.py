import os
import random
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import compute, filtering, load, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import scatter, save

TIME_POINT = {
    False: 50,
    True: 35
}
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 2
STD = 0.25
CELLS_DISTANCE = 7.0


def compute_fibers_densities(_simulations, _low_connectivity):
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

    return _fibers_densities


def main(_low_connectivity=True):
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, TIME_POINT[_low_connectivity])
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=True,
        _is_low_connectivity=_low_connectivity,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_heterogeneity(_simulations, _std=STD)
    _simulations = filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)

    _fibers_densities = compute_fibers_densities(_simulations, _low_connectivity)

    _communicated_correlations_array = []
    _non_communicated_correlations_array = []
    for _communicated_index in range(len(_simulations)):
        _communicated_simulation = _simulations[_communicated_index]
        _communicated_left_cell_fibers_densities = _fibers_densities[(_communicated_simulation, 'left_cell')]
        _communicated_right_cell_fibers_densities = _fibers_densities[(_communicated_simulation, 'right_cell')]
        _communicated_correlation = compute_lib.correlation(
            compute_lib.derivative(_communicated_left_cell_fibers_densities, _n=DERIVATIVE),
            compute_lib.derivative(_communicated_right_cell_fibers_densities, _n=DERIVATIVE)
        )
        for _non_communicated_index in range(len(_simulations)):
            if _communicated_index != _non_communicated_index:
                _non_communicated_simulation = _simulations[_non_communicated_index]
                for _communicated_cell_id, _non_communicated_cell_id in product(['left_cell', 'right_cell'],
                                                                                ['left_cell', 'right_cell']):
                    _communicated_fibers_densities = \
                        _fibers_densities[(_communicated_simulation, _communicated_cell_id)]
                    _non_communicated_fibers_densities = \
                        _fibers_densities[(_non_communicated_simulation, _non_communicated_cell_id)]
                    _non_communicated_correlations_array.append(compute_lib.correlation(
                        compute_lib.derivative(_communicated_fibers_densities, _n=DERIVATIVE),
                        compute_lib.derivative(_non_communicated_fibers_densities, _n=DERIVATIVE)
                    ))
                    _communicated_correlations_array.append(_communicated_correlation)

    # points plot
    _fig = scatter.create_plot(
        _x_array=[_communicated_correlations_array],
        _y_array=[_non_communicated_correlations_array],
        _names_array=['Distance ' + str(CELLS_DISTANCE)],
        _modes_array=['markers'],
        _show_legend_array=[False],
        _x_axis_title='Communicated Pair Correlation',
        _y_axis_title='Non-communicated Pair Correlation'
    )

    _fig = scatter.add_line(
        _fig=_fig,
        _x1=-1, _y1=-1, _x2=1, _y2=1,
        _name='y = x',
        _color='red',
        _showlegend=False
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_derivative_' + str(DERIVATIVE) + '_low_con_' + str(_low_connectivity)
    )

    _communicated_minus_non_communicated = \
        np.array(_communicated_correlations_array) - np.array(_non_communicated_correlations_array)
    _communicated_count = len(_communicated_minus_non_communicated[_communicated_minus_non_communicated > 0])
    _communicated_percentages = round(_communicated_count / len(_communicated_minus_non_communicated), 10)
    _wilcoxon_rank = wilcoxon(_communicated_minus_non_communicated)

    print('Communicated Pair:', str(_communicated_percentages * 100) + '%')
    print('Wilcoxon:', _wilcoxon_rank)

    # TODO: create bot plot ?


if __name__ == '__main__':
    main()

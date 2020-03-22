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
STD = 0.5
CELLS_DISTANCE = 5.0


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


def main(_low_connectivity=False):
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

    _master_correlations_array = []
    _slave_correlations_array = []
    for _master_index in range(len(_simulations)):
        _master_simulation = _simulations[_master_index]
        _master_left_cell_fibers_densities = _fibers_densities[(_master_simulation, 'left_cell')]
        _master_right_cell_fibers_densities = _fibers_densities[(_master_simulation, 'right_cell')]
        _master_correlation = compute_lib.correlation(
            compute_lib.derivative(_master_left_cell_fibers_densities, _n=DERIVATIVE),
            compute_lib.derivative(_master_right_cell_fibers_densities, _n=DERIVATIVE)
        )
        for _slave_index in range(len(_simulations)):
            if _master_index != _slave_index:
                _slave_simulation = _simulations[_slave_index]
                for _master_cell_id, _slave_cell_id in product(['left_cell', 'right_cell'],
                                                               ['left_cell', 'right_cell']):
                    _master_fibers_densities = _fibers_densities[(_master_simulation, _master_cell_id)]
                    _slave_fibers_densities = _fibers_densities[(_slave_simulation, _slave_cell_id)]
                    _slave_correlations_array.append(compute_lib.correlation(
                        compute_lib.derivative(_master_fibers_densities, _n=DERIVATIVE),
                        compute_lib.derivative(_slave_fibers_densities, _n=DERIVATIVE)
                    ))
                    _master_correlations_array.append(_master_correlation)

    # points plot
    _fig = scatter.create_plot(
        _x_array=[_master_correlations_array],
        _y_array=[_slave_correlations_array],
        _names_array=['Distance ' + str(CELLS_DISTANCE)],
        _modes_array=['markers'],
        _show_legend_array=[False],
        _x_axis_title='Master Network Correlation',
        _y_axis_title='Slave Network Correlation',
        _title='Master vs. Slave Network Correlations'
    )

    _fig = scatter.add_line(
        _fig=_fig,
        _x1=-1, _y1=-1, _x2=1, _y2=1,
        _name='y = x',
        _color='red',
        _showlegend=True
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_derivative_' + str(DERIVATIVE)
    )

    _master_minus_slave = np.array(_master_correlations_array) - np.array(_slave_correlations_array)
    _master_count = len(_master_minus_slave[_master_minus_slave > 0])
    _master_percentages = round(_master_count / len(_master_minus_slave), 10)
    _wilcoxon_rank = wilcoxon(_master_minus_slave)

    print('Master Network:', str(_master_percentages * 100) + '%')
    print('Wilcoxon:', _wilcoxon_rank)

    # TODO: create bot plot ?


if __name__ == '__main__':
    main()

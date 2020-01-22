import os

import numpy as np
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.simulations import load, filtering, compute, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import box, save

TIME_POINTS = 50
OFFSET_X = 0
OFFSET_Y = 0


def run(_simulation):
    _fibers_densities = []
    _normalization = load.normalization(_simulation)
    for _time_point in range(TIME_POINTS):
        _direction_fibers_densities = []
        for _direction in ['left', 'right', 'up', 'down']:
            _fibers_density = compute.roi_fibers_density_time_point(
                _simulation=_simulation,
                _length_x=ROI_HEIGHT if _direction in ['up', 'down'] else ROI_WIDTH,
                _length_y=ROI_WIDTH if _direction in ['up', 'down'] else ROI_HEIGHT,
                _offset_x=OFFSET_Y if _direction in ['up', 'down'] else OFFSET_X,
                _offset_y=OFFSET_X if _direction in ['up', 'down'] else OFFSET_Y,
                _cell_id='cell',
                _direction=_direction,
                _time_point=_time_point
            )
            _normalized_fibers_density = compute_lib.z_score(
                _fibers_density,
                _normalization['average'],
                _normalization['std']
            )
            _direction_fibers_densities.append(_normalized_fibers_density)

        _fibers_densities.append(np.mean(_direction_fibers_densities))

    return _fibers_densities


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINTS)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=True,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _correlations_array_derivative_0 = []
    _correlations_array_derivative_1 = []
    _correlations_array_derivative_2 = []
    for _index_1 in range(len(_simulations)):
        _fibers_densities_simulation_1 = run(_simulations[_index_1])
        for _index_2 in range(_index_1 + 1, len(_simulations)):
            print(_simulations[_index_1], _simulations[_index_2])
            _fibers_densities_simulation_2 = run(_simulations[_index_2])
            _correlations_array_derivative_0.append(compute_lib.correlation(
                compute_lib.derivative(_fibers_densities_simulation_1, _n=0),
                compute_lib.derivative(_fibers_densities_simulation_2, _n=0)
            ))
            _correlations_array_derivative_1.append(compute_lib.correlation(
                compute_lib.derivative(_fibers_densities_simulation_1, _n=1),
                compute_lib.derivative(_fibers_densities_simulation_2, _n=1)
            ))
            _correlations_array_derivative_2.append(compute_lib.correlation(
                compute_lib.derivative(_fibers_densities_simulation_1, _n=2),
                compute_lib.derivative(_fibers_densities_simulation_2, _n=2)
            ))

    # plot
    _fig = box.create_plot(
        _y_array=[_correlations_array_derivative_0, _correlations_array_derivative_1, _correlations_array_derivative_2],
        _names_array=['Fibers Densities', '1st Derivative', '2nd Derivative'],
        _x_axis_title='Fibers Densities by Derivatives',
        _y_axis_title='Correlation',
        _title='Single Cells Correlation by Fibers Densities Derivative'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )

    # wilcoxon
    print('Derivative 0:', wilcoxon(_correlations_array_derivative_0))
    print('Derivative 1:', wilcoxon(_correlations_array_derivative_1))
    print('Derivative 2:', wilcoxon(_correlations_array_derivative_2))


if __name__ == '__main__':
    main()

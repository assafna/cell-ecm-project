import os

import numpy as np
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.simulations import compute, filtering, load, organize, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import box, save

TIME_POINTS = 35
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 2
CELLS_DISTANCES = [4.0, 5.0, 7.0]


def run(_simulation, _cell_id):
    _inside_fibers_densities = compute.roi_fibers_density_by_time({
        'simulation': _simulation,
        'length_x': ROI_WIDTH,
        'length_y': ROI_HEIGHT,
        'offset_x': OFFSET_X,
        'offset_y': OFFSET_Y,
        'cell_id': _cell_id,
        'direction': 'inside',
        'time_points': TIME_POINTS
    })
    _outside_fibers_densities = compute.roi_fibers_density_by_time({
        'simulation': _simulation,
        'length_x': ROI_WIDTH,
        'length_y': ROI_HEIGHT,
        'offset_x': OFFSET_X,
        'offset_y': OFFSET_Y,
        'cell_id': _cell_id,
        'direction': 'outside',
        'time_points': TIME_POINTS
    })

    return compute_lib.correlation(
        compute_lib.derivative(_inside_fibers_densities, _n=DERIVATIVE),
        compute_lib.derivative(_outside_fibers_densities, _n=DERIVATIVE)
    )


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINTS)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=False,
        _is_low_connectivity=True,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_distances(_simulations, _distances=CELLS_DISTANCES)
    _simulations_by_distances = organize.by_distances(_simulations)
    _correlations_by_distances = {}
    for _distance in _simulations_by_distances:
        for _simulation in _simulations_by_distances[_distance]:
            print(_simulation)
            for _cell_id in ['left_cell', 'right_cell']:
                _cell_correlation = run(_simulation, _cell_id)
                if _distance in _correlations_by_distances:
                    _correlations_by_distances[_distance].append(_cell_correlation)
                else:
                    _correlations_by_distances[_distance] = [_cell_correlation]

    # plot
    _fig = box.create_plot(
        _y_array=[_values for _values in _correlations_by_distances.values()],
        _names_array=['Distance ' + str(_distance) for _distance in _correlations_by_distances.keys()],
        _x_axis_title='Cells Distances',
        _y_axis_title='Correlation',
        _title='Inside vs. Outside Correlation by Cell Distance'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_low_connectivity_derivative_' + str(DERIVATIVE)
    )


if __name__ == '__main__':
    main()

import os

import numpy as np
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.simulations import compute, filtering, load, paths, organize
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import save, box

MINIMUM_TIME_POINTS = 50
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 2
STD = 0.5
CELLS_DISTANCES = [4.0, 5.0, 7.0, 9.0]
DIRECTION = 'inside'


def run(_simulations):
    _communicated_correlations_array = []
    _non_communicated_correlations_array = []
    for _communicated_index in range(len(_simulations)):
        _communicated_simulation = _simulations[_communicated_index]
        _communicated_left_cell_fibers_densities = compute.roi_fibers_density_by_time({
            'simulation': _communicated_simulation,
            'length_x': ROI_WIDTH,
            'length_y': ROI_HEIGHT,
            'offset_x': OFFSET_X,
            'offset_y': OFFSET_Y,
            'cell_id': 'left_cell',
            'direction': DIRECTION,
            'time_points': MINIMUM_TIME_POINTS
        })
        _communicated_right_cell_fibers_densities = compute.roi_fibers_density_by_time({
            'simulation': _communicated_simulation,
            'length_x': ROI_WIDTH,
            'length_y': ROI_HEIGHT,
            'offset_x': OFFSET_X,
            'offset_y': OFFSET_Y,
            'cell_id': 'right_cell',
            'direction': DIRECTION,
            'time_points': MINIMUM_TIME_POINTS
        })
        _communicated_correlation = compute_lib.correlation(
            compute_lib.derivative(_communicated_left_cell_fibers_densities, _n=DERIVATIVE),
            compute_lib.derivative(_communicated_right_cell_fibers_densities, _n=DERIVATIVE)
        )
        for _non_communicated_index in range(_communicated_index + 1, len(_simulations)):
            _non_communicated_simulation = _simulations[_non_communicated_index]
            print(_communicated_simulation, _non_communicated_simulation, sep='\t')
            _non_communicated_left_cell_fibers_densities = compute.roi_fibers_density_by_time({
                'simulation': _non_communicated_simulation,
                'length_x': ROI_WIDTH,
                'length_y': ROI_HEIGHT,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'cell_id': 'left_cell',
                'direction': DIRECTION,
                'time_points': MINIMUM_TIME_POINTS
            })
            _non_communicated_correlations_array.append(compute_lib.correlation(
                compute_lib.derivative(_communicated_left_cell_fibers_densities, _n=DERIVATIVE),
                compute_lib.derivative(_non_communicated_left_cell_fibers_densities, _n=DERIVATIVE)
            ))
            _communicated_correlations_array.append(_communicated_correlation)

    return np.array(_communicated_correlations_array) - np.array(_non_communicated_correlations_array)


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, MINIMUM_TIME_POINTS)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=True,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_heterogeneity(_simulations, _std=STD)
    _simulations = filtering.by_distances(_simulations, _distances=CELLS_DISTANCES)
    _simulations_by_distance = organize.by_distances(_simulations)
    _y_arrays = [run(_simulations_by_distance[_distance]) for _distance in _simulations_by_distance]

    # plot
    _fig = box.create_plot(
        _y_array=_y_arrays,
        _names_array=CELLS_DISTANCES,
        _x_axis_title='Cells Distance',
        _y_axis_title='Communicated minus Non-communicated',
        _title='Communicated minus Non-communicated Pair Correlations by Cell Distance'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )

    # wilcoxon
    for _distance, _y in zip(CELLS_DISTANCES, _y_arrays):
        print('Wilcoxon distance', _distance, wilcoxon(_y))


if __name__ == '__main__':
    main()

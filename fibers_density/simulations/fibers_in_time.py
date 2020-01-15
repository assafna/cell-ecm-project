import os

import numpy as np

from libs import compute_lib
from libs.simulations import load, compute, filtering, paths
from libs.simulations.config import CELL_DIAMETER, ROI_WIDTH, ROI_HEIGHT
from plotting import scatter, save

TIME_POINTS = 50
CELLS_DISTANCE = 7.0
OFFSET_X = 0
OFFSET_Y = 0


def run(_simulations, _cell_id, _directions):
    _fibers_densities = list([None] * TIME_POINTS)
    for _simulation in _simulations:
        print(_simulation)
        _time_point_index = 0
        _normalization = load.normalization(_simulation)
        for _time_point in range(TIME_POINTS):
            _direction_fibers_densities = []
            for _direction in _directions:
                _fibers_density = compute.roi_fibers_density_time_point(
                    _simulation=_simulation,
                    _length_x=ROI_HEIGHT if _direction in ['up', 'down'] else ROI_WIDTH,
                    _length_y=ROI_WIDTH if _direction in ['up', 'down'] else ROI_HEIGHT,
                    _offset_x=OFFSET_Y if _direction in ['up', 'down'] else OFFSET_X,
                    _offset_y=OFFSET_X if _direction in ['up', 'down'] else OFFSET_Y,
                    _cell_id=_cell_id,
                    _direction=_direction,
                    _time_point=_time_point
                )
                _normalized_fibers_density = compute_lib.z_score(
                    _fibers_density,
                    _normalization['average'],
                    _normalization['std']
                )
                _direction_fibers_densities.append(_normalized_fibers_density)

            if _fibers_densities[_time_point_index] is None:
                _fibers_densities[_time_point_index] = [np.mean(_direction_fibers_densities)]
            else:
                _fibers_densities[_time_point_index].append(np.mean(_direction_fibers_densities))
            _time_point_index += 1

    return _fibers_densities


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINTS)

    # single cell
    _simulations_single_cells = filtering.by_categories(
        _simulations,
        _is_single_cell=True,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _fibers_densities_single_cells = run(
        _simulations=_simulations_single_cells,
        _cell_id='cell',
        _directions=['left', 'right', 'up', 'down']
    )

    # pairs
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations_pairs = filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)
    _fibers_densities_pairs_left = run(
        _simulations=_simulations_pairs,
        _cell_id='left_cell',
        _directions=['inside']
    )
    _fibers_densities_pairs_right = run(
        _simulations=_simulations_pairs,
        _cell_id='right_cell',
        _directions=['inside']
    )

    # plot
    _fig = scatter.create_error_bars_plot(
        _x_array=[list(range(TIME_POINTS))] * 3,
        _y_array=[_fibers_densities_single_cells, _fibers_densities_pairs_left, _fibers_densities_pairs_right],
        _names_array=['Single Cell', 'Pair Left Cell', 'Pair Right Cell'],
        _mode_array=['lines+markers'] * 3,
        _dash_array=['dash'] + ['solid'] * 2,
        _x_axis_title='Time (cell contraction percentages)',
        _y_axis_title='Fibers Density Z-score',
        _title='Fibers Densities vs. Time'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()

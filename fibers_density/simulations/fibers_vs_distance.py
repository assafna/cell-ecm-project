import os

import numpy as np

from libs import compute_lib
from libs.simulations import load, filtering, compute, paths
from libs.simulations.config import CELL_DIAMETER, ROI_WIDTH, ROI_HEIGHT
from plotting import scatter, save

TIME_POINT = 50
CELLS_DISTANCE = 7.0
OFFSET_X_END = CELLS_DISTANCE * CELL_DIAMETER - (CELL_DIAMETER / 2) - ROI_WIDTH
OFFSET_X_STEP = CELL_DIAMETER / 8
OFFSET_Y = 0


def run(_simulations, _cell_id, _directions):
    _fibers_densities = list([None] * len(np.arange(start=0, stop=OFFSET_X_END, step=OFFSET_X_STEP)))
    for _simulation in _simulations:
        print(_simulation)
        _offset_index = 0
        _normalization = load.normalization(_simulation)
        for _offset_x in np.arange(start=0, stop=OFFSET_X_END, step=OFFSET_X_STEP):
            _direction_fibers_densities = []
            for _direction in _directions:
                _fibers_density = compute.roi_fibers_density_time_point(
                    _simulation=_simulation,
                    _length_x=ROI_HEIGHT if _direction in ['up', 'down'] else ROI_WIDTH,
                    _length_y=ROI_WIDTH if _direction in ['up', 'down'] else ROI_HEIGHT,
                    _offset_x=OFFSET_Y if _direction in ['up', 'down'] else _offset_x,
                    _offset_y=_offset_x if _direction in ['up', 'down'] else OFFSET_Y,
                    _cell_id=_cell_id,
                    _direction=_direction,
                    _time_point=TIME_POINT
                )
                _normalized_fibers_density = compute_lib.z_score(
                    _fibers_density,
                    _normalization['average'],
                    _normalization['std']
                )
                _direction_fibers_densities.append(_normalized_fibers_density)

            if _fibers_densities[_offset_index] is None:
                _fibers_densities[_offset_index] = [np.mean(_direction_fibers_densities)]
            else:
                _fibers_densities[_offset_index].append(np.mean(_direction_fibers_densities))
            _offset_index += 1

    return _fibers_densities


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINT)

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
    _simulations_pairs = filtering.by_distances(_simulations, _distances=[CELLS_DISTANCE])
    _fibers_densities_pairs = run(
        _simulations=_simulations_pairs,
        _cell_id='left_cell',
        _directions=['inside']
    )

    # plot
    _fig = scatter.create_error_bars_plot(
        _x_array=[np.arange(start=0, stop=OFFSET_X_END, step=OFFSET_X_STEP) / CELL_DIAMETER] * 2,
        _y_array=[_fibers_densities_single_cells, _fibers_densities_pairs],
        _names_array=['Single Cell', 'Pair'],
        _mode_array=['lines+markers'] * 2,
        _dash_array=['dash', 'solid'],
        _x_axis_title='Distance from Left Cell (cell size)',
        _y_axis_title='Fibers Density Z-score',
        _title='Fibers Densities vs. Distance from Cell'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()

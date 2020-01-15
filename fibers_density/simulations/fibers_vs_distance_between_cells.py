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


def main():
    _simulations = load.structured()
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINT)
    _fibers_densities = list([None] * len(np.arange(start=0, stop=OFFSET_X_END, step=OFFSET_X_STEP)))
    for _simulation in _simulations:
        print(_simulation)
        _offset_index = 0
        _normalization = load.normalization(_simulation)
        for _offset_x in np.arange(start=0, stop=OFFSET_X_END, step=OFFSET_X_STEP):
            _fibers_density = compute.roi_fibers_density_time_point(
                _simulation=_simulation,
                _length_x=ROI_WIDTH,
                _length_y=ROI_HEIGHT,
                _offset_x=_offset_x,
                _offset_y=OFFSET_Y,
                _cell_id='left_cell',
                _direction='right',
                _time_point=TIME_POINT
            )
            _normalized_fibers_density = compute_lib.z_score(
                _fibers_density,
                _normalization['average'],
                _normalization['std']
            )

            if _fibers_densities[_offset_index] is None:
                _fibers_densities[_offset_index] = [_normalized_fibers_density]
            else:
                _fibers_densities[_offset_index].append(_normalized_fibers_density)
            _offset_index += 1

    _fig = scatter.create_error_bars_plot(
        _x_array=[np.arange(start=0, stop=OFFSET_X_END, step=OFFSET_X_STEP) / CELL_DIAMETER],
        _y_array=[_fibers_densities],
        _names_array=['Time-Point Last'],
        _modes_array=['lines+markers'],
        _dashes_array=['solid'],
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

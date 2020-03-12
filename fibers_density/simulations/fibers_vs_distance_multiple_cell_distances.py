import os
from multiprocessing.pool import Pool

import numpy as np
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import load, filtering, compute, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import scatter, save, edit

TIME_POINT = 50
CELLS_DISTANCES = [5, 7, 9]
OFFSET_X_STEP = 0.2
OFFSET_Y = 0


def compute_simulations_fibers_densities(_simulations, _offsets_x):
    _arguments = []
    for _simulation in _simulations:
        for _offset_x in _offsets_x:
            _arguments.append({
                'simulation': _simulation,
                'length_x': ROI_WIDTH,
                'length_y': ROI_HEIGHT,
                'offset_x': _offset_x,
                'offset_y': OFFSET_Y,
                'cell_id': 'left_cell',
                'direction': 'inside',
                'time_point': TIME_POINT
            })

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.roi_fibers_density_time_point, _arguments),
                total=len(_arguments), desc='Computing Rois & Fibers Densities'):
            _fibers_densities[
                (_keys['simulation'], _keys['offset_x'])] = _value
        _p.close()
        _p.join()

    return _fibers_densities


def main():
    _x_array = []
    _y_array = []
    _names_array = []
    for _distance in CELLS_DISTANCES:
        print('Cells Distance ' + str(_distance))
        _simulations = load.structured()
        _simulations = filtering.by_categories(
            _simulations,
            _is_single_cell=False,
            _is_heterogeneity=False,
            _is_low_connectivity=False,
            _is_causality=False,
            _is_dominant_passive=False
        )
        _simulations = filtering.by_distance(_simulations, _distance=_distance)
        _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINT)

        _offset_x_end = _distance - (TIME_POINT / 100) - ROI_WIDTH
        _offsets_x = np.arange(start=0, stop=_offset_x_end + OFFSET_X_STEP, step=OFFSET_X_STEP)

        _fibers_densities = compute_simulations_fibers_densities(_simulations, _offsets_x)

        _cells_distance_fibers_densities = list([None] * len(_offsets_x))
        for _simulation in _simulations:
            _offset_index = 0
            _normalization = load.normalization(_simulation)
            for _offset_x in _offsets_x:
                _fibers_density = _fibers_densities[(_simulation, _offset_x)]

                _normalized_fibers_density = compute_lib.z_score(
                    _fibers_density,
                    _normalization['average'],
                    _normalization['std']
                )

                if _cells_distance_fibers_densities[_offset_index] is None:
                    _cells_distance_fibers_densities[_offset_index] = [_normalized_fibers_density]
                else:
                    _cells_distance_fibers_densities[_offset_index].append(_normalized_fibers_density)
                _offset_index += 1

        _x_array.append(_offsets_x)
        _y_array.append(_cells_distance_fibers_densities)
        _names_array.append('Cells Distance ' + str(_distance))

    _fig = scatter.create_error_bars_plot(
        _x_array=_x_array,
        _y_array=_y_array,
        _names_array=_names_array,
        _modes_array=['lines+markers'] * len(CELLS_DISTANCES),
        _dashes_array=['solid'] * len(CELLS_DISTANCES),
        _x_axis_title='Distance from Left Cell (cell size)',
        _y_axis_title='Fibers Density Z-score'
    )

    _fig = edit.update_y_axis(
        _fig=_fig,
        _range=[-1.5, 17]
    )

    _fig = edit.update_x_axis(
        _fig=_fig,
        _range=[-0.25, 8]
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()

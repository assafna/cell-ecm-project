import os
import sys
from multiprocessing.pool import Pool

import numpy as np
from scipy.stats import stats

import libs.compute_lib
from libs import compute_lib
from libs.experiments import load, filtering, compute, paths, organize
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT, CELL_DIAMETER_IN_MICRONS
from plotting import scatter, save

MINIMUM_TIME_POINTS = sys.maxsize
OFFSET_X = 0
OFFSET_Y = CELL_DIAMETER_IN_MICRONS * 0
OFFSET_Z = 0
DERIVATIVE = 1
CELLS_DISTANCES = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
DIRECTION = 'inside'


def main():
    _experiments = load.experiments_groups_as_tuples(['SN16', 'SN41'])
    _experiments = filtering.by_distances(_experiments, CELLS_DISTANCES)
    # _experiments = filtering.by_band(_experiments)

    # prepare data in mp
    _fibers_densities = {}
    _arguments = []
    for _tuple in _experiments:
        _experiment, _series, _group = _tuple
        _arguments.append((_experiment, _series, _group, ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT,
                           OFFSET_X, OFFSET_Y, OFFSET_Z, DIRECTION, MINIMUM_TIME_POINTS))
    _p = Pool()
    _answers = _p.starmap(compute.roi_fibers_density_by_time_pairs, _arguments)
    _p.close()
    for _index in range(len(_arguments)):
        _tuple = _experiments[_index]
        _answer = _answers[_index]
        _fibers_densities[_tuple] = _answer

    _experiments_by_distance = organize.by_cells_distance(_experiments)
    _fibers_densities_by_distance = {}
    _change_in_fibers_densities_by_distance = {}
    for _distance in _experiments_by_distance:
        _fibers_densities_array = []
        _change_in_fibers_densities_array = []
        for _tuple in _experiments_by_distance[_distance]:
            print(_tuple)
            _experiment, _series_id, _ = _tuple
            _series_normalization = load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))
            for _cell_id in ['left_cell', 'right_cell']:
                _cell_fibers_densities = _fibers_densities[_tuple][_cell_id]

                # fix if found nan
                if True in np.isnan(_cell_fibers_densities):
                    _cell_fibers_densities = _cell_fibers_densities[:np.where(np.isnan(_cell_fibers_densities))[0][0]]

                _z_score_fibers_density = libs.compute_lib.z_score_fibers_densities_array(
                    _cell_fibers_densities, _series_normalization
                )
                _fibers_densities_array += _z_score_fibers_density
                _change_in_fibers_densities_array += [0] + compute_lib.derivative(
                    _z_score_fibers_density, _n=DERIVATIVE
                )
        _fibers_densities_by_distance[_distance] = _fibers_densities_array
        _change_in_fibers_densities_by_distance[_distance] = _change_in_fibers_densities_array

    _fig = scatter.create_plot(
        _x_array=_fibers_densities_by_distance.values(),
        _y_array=_change_in_fibers_densities_by_distance.values(),
        _names_array=['Distance ' + str(_distance) for _distance in _fibers_densities_by_distance.keys()],
        _modes_array=['markers'] * len(_fibers_densities_by_distance.keys()),
        _showlegend_array=[True] * len(_fibers_densities_by_distance.keys()),
        _x_axis_title='Fibers Densities Z-Score',
        _y_axis_title='Change in Fibers Densities Z-Score',
        _title='Fibers Densities vs. Change in Fibers Densities - ' + DIRECTION.capitalize()
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename=DIRECTION + '_points'
    )

    # line of best fit
    _best_fit_lines_x_array = []
    _best_fit_lines_y_array = []
    for _distance in _experiments_by_distance:
        _x_array = _fibers_densities_by_distance[_distance]
        _y_array = _change_in_fibers_densities_by_distance[_distance]
        _slope, _intercept, _r_value, _p_value, _std_err = stats.linregress(_x_array, _y_array)
        _x1, _x2 = min(_x_array), max(_x_array)
        _y1, _y2 = _slope * _x1 + _intercept, _slope * _x2 + _intercept
        _best_fit_lines_x_array.append([_x1, _x2])
        _best_fit_lines_y_array.append([_y1, _y2])

    _fig = scatter.create_plot(
        _x_array=_best_fit_lines_x_array,
        _y_array=_best_fit_lines_y_array,
        _names_array=['Distance ' + str(_distance) for _distance in _experiments_by_distance],
        _modes_array=['lines'] * len(_experiments_by_distance),
        _showlegend_array=[True] * len(_experiments_by_distance),
        _x_axis_title='Fibers Densities Z-Score',
        _y_axis_title='Change in Fibers Densities Z-Score',
        _title='Fibers Densities vs. Change in Fibers Densities - ' + DIRECTION.capitalize() + ' - Line of Best Fit'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename=DIRECTION + '_best_fit'
    )


if __name__ == '__main__':
    main()

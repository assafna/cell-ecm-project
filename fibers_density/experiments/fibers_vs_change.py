import os
import sys
from multiprocessing.pool import Pool

import numpy as np
import seaborn
from scipy.stats import stats, pearsonr
import plotly.graph_objs as go

import libs.compute_lib
from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import load, filtering, compute, paths, organize
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import scatter, save, heatmap

MINIMUM_TIME_POINTS = sys.maxsize
OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
DERIVATIVE = 1
CELLS_DISTANCES = range(6, 10)
DIRECTION = 'inside'


def main():
    _experiments = load.experiments_groups_as_tuples(['SN16'])
    _experiments = filtering.by_distances(_experiments, CELLS_DISTANCES)
    _experiments = filtering.by_real_cells(_experiments)
    _experiments = filtering.by_band(_experiments)

    # _experiments.remove(('SN41', 6, 'cells_1_2'))

    # prepare data in mp
    _fibers_densities = {}
    _arguments = []
    for _tuple in _experiments:
        _experiment, _series, _group = _tuple
        _arguments.append((_experiment, _series, _group, ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH,
                           OFFSET_X, OFFSET_Y, OFFSET_Z, DIRECTION, MINIMUM_TIME_POINTS))
    _p = Pool(CPUS_TO_USE)
    _answers = _p.starmap(compute.roi_fibers_density_by_time_pairs, _arguments)
    _p.close()
    for _index in range(len(_arguments)):
        _tuple = _experiments[_index]
        _answer = _answers[_index]
        _fibers_densities[_tuple] = _answer

    _experiments_by_distance = organize.by_cells_distance(_experiments)
    _fibers_densities_by_distance = {}
    _change_in_fibers_densities_by_distance = {}
    _heatmap_fibers = []
    _heatmap_fibers_change = []
    for _distance in _experiments_by_distance:
        _fibers_densities_array = []
        _change_in_fibers_densities_array = []
        for _tuple in _experiments_by_distance[_distance]:
            print(_tuple)
            _experiment, _series_id, _group = _tuple
            _series_normalization = load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))
            _series_normalization = [_series_normalization['average'], _series_normalization['std']]
            for _cell_id in ['left_cell', 'right_cell']:
                _cell_fibers_densities = _fibers_densities[_tuple][_cell_id]
                _cell_fibers_densities = compute.remove_blacklist(
                    _experiment, _series_id, _group, _cell_fibers_densities)
                _cell_fibers_densities = compute.longest_fibers_densities_ascending_sequence(_cell_fibers_densities)

                # fix if found nan
                if True in np.isnan(_cell_fibers_densities):
                    _cell_fibers_densities = _cell_fibers_densities[:np.where(np.isnan(_cell_fibers_densities))[0][0]]

                # not enough data
                if len(_cell_fibers_densities) < DERIVATIVE + 1:
                    continue

                _z_score_fibers_density = libs.compute_lib.z_score_fibers_densities_array(
                    _cell_fibers_densities, _series_normalization
                )
                if _experiment == 'SN41':
                    for _start_index in [0, 1, 2]:
                        _fibers_densities_array += _z_score_fibers_density[_start_index::3][DERIVATIVE:]
                        _change_in_fibers_densities_array += compute_lib.derivative(
                            _z_score_fibers_density[_start_index::3], _n=DERIVATIVE
                        )
                else:
                    _fibers_densities_array += _z_score_fibers_density[DERIVATIVE:]
                    _change_in_fibers_densities_array += compute_lib.derivative(
                        _z_score_fibers_density, _n=DERIVATIVE
                    )
        _fibers_densities_by_distance[_distance] = _fibers_densities_array
        _change_in_fibers_densities_by_distance[_distance] = _change_in_fibers_densities_array
        _heatmap_fibers += _fibers_densities_array
        _heatmap_fibers_change += _change_in_fibers_densities_array

    print(pearsonr(_heatmap_fibers, _heatmap_fibers_change))

    _x_labels_start = -3
    _x_labels_end = 15
    _y_labels_start = -2
    _y_labels_end = 3
    _x_bins = 10
    _y_bins = 5
    _x_shape = int(round((_x_labels_end - _x_labels_start) * _x_bins))
    _y_shape = int(round((_y_labels_end - _y_labels_start) * _y_bins))
    _total_points = 0
    _z_array = np.zeros(shape=(_x_shape, _y_shape))
    for _x, _y in zip(_heatmap_fibers, _heatmap_fibers_change):
        _x_rounded, _y_rounded = int(round(_x * _x_bins)), int(round(_y * _y_bins))
        _x_index, _y_index = int(_x_rounded - _x_labels_start * _x_bins), int(_y_rounded - _y_labels_start * _y_bins)
        if 0 <= _x_index < _z_array.shape[0] and 0 <= _y_index < _z_array.shape[1]:
            _z_array[_x_index][_y_index] += 1
            _total_points += 1
    _z_array = _z_array / _total_points

    # _z_array = list(zip(*_z_array))

    # second normalization
    _z_array_normalized_by_cols = np.zeros(shape=np.array(_z_array).shape)
    for _fibers_index, _fibers_density_z_score in enumerate(_z_array):
        _sum = np.sum(_fibers_density_z_score)
        for _change_index, _change_z_score in enumerate(_fibers_density_z_score):
            _z_array_normalized_by_cols[_fibers_index][_change_index] = (_change_z_score / _sum) if _sum != 0 else 0

    _z_array_normalized_by_cols = list(zip(*_z_array_normalized_by_cols))

    _fig = heatmap.create_plot(
        _x_labels=np.arange(start=_x_labels_start, stop=_x_labels_end, step=1 / _x_bins),
        _y_labels=np.arange(start=_y_labels_start, stop=_y_labels_end, step=1 / _y_bins),
        _z_array=_z_array_normalized_by_cols,
        _x_axis_title='Fibers Densities Z-Score',
        _y_axis_title='Change in Fibers Densities Z-Score',
        # _color_scale=seaborn.light_palette("navy", reverse=True).as_hex(),
        _color_scale='Viridis',
        _zmax=0.015
    )

    # line of best fit
    # _best_fit_lines_x_array = []
    # _best_fit_lines_y_array = []
    # _x_array = _heatmap_fibers
    # _y_array = _heatmap_fibers_change
    # _slope, _intercept, _r_value, _p_value, _std_err = stats.linregress(_x_array, _y_array)
    # _x1, _x2 = max(_x_labels_start, min(_x_array)), min(_x_labels_end, max(_x_array))
    # _y1, _y2 = _slope * _x1 + _intercept, _slope * _x2 + _intercept
    # _best_fit_lines_x_array.append([_x1, _x2])
    # _best_fit_lines_y_array.append([_y1, _y2])
    #
    # _fig.add_trace(go.Scatter(
    #     x=_best_fit_lines_x_array[0],
    #     y=_best_fit_lines_y_array[0],
    #     name=None,
    #     mode='lines',
    #     showlegend=False
    # ))

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename=DIRECTION + '_points'
    )

    # _fig = scatter.create_plot(
    #     _x_array=_fibers_densities_by_distance.values(),
    #     _y_array=_change_in_fibers_densities_by_distance.values(),
    #     _names_array=['Distance ' + str(_distance) for _distance in _fibers_densities_by_distance.keys()],
    #     _modes_array=['markers'] * len(_fibers_densities_by_distance.keys()),
    #     _showlegend_array=[True] * len(_fibers_densities_by_distance.keys()),
    #     _x_axis_title='Fibers Densities Z-Score',
    #     _y_axis_title='Change in Fibers Densities Z-Score',
    #     _title='Fibers Densities vs. Change in Fibers Densities - ' + DIRECTION.capitalize()
    # )
    #
    # save.to_html(
    #     _fig=_fig,
    #     _path=os.path.join(paths.PLOTS, save.get_module_name()),
    #     _filename=DIRECTION + '_points'
    # )

    # line of best fit
    # _best_fit_lines_x_array = []
    # _best_fit_lines_y_array = []
    # for _distance in _experiments_by_distance:
    #     _x_array = _fibers_densities_by_distance[_distance]
    #     _y_array = _change_in_fibers_densities_by_distance[_distance]
    #     _slope, _intercept, _r_value, _p_value, _std_err = stats.linregress(_x_array, _y_array)
    #     _x1, _x2 = min(_x_array), max(_x_array)
    #     _y1, _y2 = _slope * _x1 + _intercept, _slope * _x2 + _intercept
    #     _best_fit_lines_x_array.append([_x1, _x2])
    #     _best_fit_lines_y_array.append([_y1, _y2])
    #
    # _fig = scatter.create_plot(
    #     _x_array=_best_fit_lines_x_array,
    #     _y_array=_best_fit_lines_y_array,
    #     _names_array=['Distance ' + str(_distance) for _distance in _experiments_by_distance],
    #     _modes_array=['lines'] * len(_experiments_by_distance),
    #     _showlegend_array=[True] * len(_experiments_by_distance),
    #     _x_axis_title='Fibers Densities Z-Score',
    #     _y_axis_title='Change in Fibers Densities Z-Score',
    #     _title='Fibers Densities vs. Change in Fibers Densities - ' + DIRECTION.capitalize() + ' - Line of Best Fit'
    # )
    #
    # save.to_html(
    #     _fig=_fig,
    #     _path=os.path.join(paths.PLOTS, save.get_module_name()),
    #     _filename=DIRECTION + '_best_fit'
    # )


if __name__ == '__main__':
    main()

import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly
import seaborn as sns
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT, CELL_DIAMETER_IN_MICRONS
from methods.experiments import export_video
from plotting import scatter, save, heatmap, contour

EXPERIMENTS = ['SN16']
EXPERIMENTS_STR = '_'.join(EXPERIMENTS)
REAL_CELLS = False
STATIC = True
BAND = False
MINIMUM_TIME_POINTS = 20
VALUES = [-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
VALUES_BY_CELL_DIAMETER = np.array(VALUES) * CELL_DIAMETER_IN_MICRONS
_OFFSET_X = 0
DERIVATIVE = 2
CELLS_DISTANCES = [6, 7, 8, 9]
DIRECTION = 'inside'
PRINT = False


def main():
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_distances(_experiments, CELLS_DISTANCES)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    if BAND:
        _experiments = filtering.by_band(_experiments)
    _experiments = filtering.by_time_points_amount(_experiments, MINIMUM_TIME_POINTS)

    # _experiments.remove(('SN41', 6, 'cells_1_2'))
    # _experiments.remove(('SN41', 2, 'cells_0_2'))

    _z_array = np.zeros(shape=(len(VALUES), len(VALUES)))
    for (_offset_y_index, _offset_y), (_offset_z_index, _offset_z) in \
            product(enumerate(VALUES_BY_CELL_DIAMETER), enumerate(VALUES_BY_CELL_DIAMETER)):

        # prepare data in mp
        _fibers_densities = {}
        _arguments = []
        _answers_keys = []
        for _tuple in _experiments:
            _experiment, _series, _group = _tuple
            _arguments.append((_experiment, _series, _group, ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT,
                               _OFFSET_X, _offset_y, _offset_z, DIRECTION, MINIMUM_TIME_POINTS, PRINT))
            _answers_keys.append((_experiment, _series, _group))

        _p = Pool(CPUS_TO_USE)
        _answers = _p.starmap(compute.roi_fibers_density_by_time_pairs, _arguments)
        _p.close()
        for _answer_key, _answer in zip(_answers_keys, _answers):
            _fibers_densities[_answer_key] = _answer

        _master_correlations_array = []
        _slave_correlations_array = []
        for _master_index in range(len(_experiments)):
            _master_tuple = _experiments[_master_index]
            _master_experiment, _master_series, _master_group = _master_tuple
            _master_left_cell_fibers_densities = \
                _fibers_densities[(_master_experiment, _master_series, _master_group)]['left_cell']
            _master_right_cell_fibers_densities = \
                _fibers_densities[(_master_experiment, _master_series, _master_group)]['right_cell']

            _master_left_cell_fibers_densities = compute.remove_blacklist(
                _master_experiment, _master_series, _master_group, _master_left_cell_fibers_densities)
            _master_right_cell_fibers_densities = compute.remove_blacklist(
                _master_experiment, _master_series, _master_group, _master_right_cell_fibers_densities)

            _master_left_cell_fibers_densities_filtered, _master_right_cell_fibers_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _master_left_cell_fibers_densities, _master_right_cell_fibers_densities
                )

            _master_correlation = compute_lib.correlation(
                compute_lib.derivative(_master_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_master_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
            )
            for _slave_index in range(len(_experiments)):
                if _master_index != _slave_index:
                    _slave_tuple = _experiments[_slave_index]
                    _slave_experiment, _slave_series, _slave_group = _slave_tuple
                    for _master_cell_id, _slave_cell_id in product(['left_cell', 'right_cell'],
                                                                   ['left_cell', 'right_cell']):
                        _master_fibers_densities = \
                            _fibers_densities[(_master_experiment, _master_series, _master_group)][_master_cell_id]
                        _slave_fibers_densities = \
                            _fibers_densities[(_slave_experiment, _slave_series, _slave_group)][_slave_cell_id]

                        _master_fibers_densities_filtered, _slave_fibers_densities_filtered = \
                            compute.longest_same_indices_shared_in_borders_sub_array(
                                _master_fibers_densities, _slave_fibers_densities
                            )

                        _slave_correlation = compute_lib.correlation(
                            compute_lib.derivative(_master_fibers_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_slave_fibers_densities_filtered, _n=DERIVATIVE)
                        )
                        _master_correlations_array.append(_master_correlation)
                        _slave_correlations_array.append(_slave_correlation)

        # compute percentages
        _master_minus_slave = np.array(_master_correlations_array) - np.array(_slave_correlations_array)
        _master_count = len(_master_minus_slave[_master_minus_slave > 0])
        _master_percentages = round(_master_count / len(_master_minus_slave), 10)

        print('z', _offset_y / CELL_DIAMETER_IN_MICRONS, 'xy', _offset_z / CELL_DIAMETER_IN_MICRONS,
              _master_percentages, sep='\t')
        _z_array[_offset_z_index, _offset_y_index] = _master_percentages

    # plot
    _fig = heatmap.create_plot(
        _x_labels=VALUES,
        _y_labels=VALUES,
        _z_array=_z_array,
        _x_axis_title='Offset in Z axis',
        _y_axis_title='Offset in XY axis',
        _color_scale=sns.color_palette('BrBG').as_hex(),
        _zmin=0,
        _zmax=1
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_' + EXPERIMENTS_STR + '_real_' + str(REAL_CELLS) + '_static_' + str(STATIC) + '_band_' +
                  str(BAND) + '_oob_heatmap'
    )

    _fig = contour.create_plot(
        _x_labels=VALUES,
        _y_labels=VALUES,
        _z_array=_z_array,
        _x_axis_title='Offset in Z axis',
        _y_axis_title='Offset in XY axis',
        _color_scale=sns.color_palette('BrBG').as_hex(),
        _zmin=0,
        _zmax=1
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_' + EXPERIMENTS_STR + '_real_' + str(REAL_CELLS) + '_static_' + str(STATIC) + '_band_' +
                  str(BAND) + '_oob_contour'
    )


if __name__ == '__main__':
    main()

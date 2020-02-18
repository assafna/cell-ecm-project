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
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from methods.experiments import export_video
from plotting import scatter, save, heatmap, contour

EXPERIMENTS = ['SN16', 'SN41']
EXPERIMENTS_STR = '_'.join(EXPERIMENTS)
REAL_CELLS = False
STATIC = False
BAND = False
VALUES_BY_CELL_DIAMETER = np.array(
    [-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
_OFFSET_X = 0
DERIVATIVE = 2
CELLS_DISTANCES = [6, 7, 8, 9]
DIRECTION = 'inside'


def main():
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    # _experiments = filtering.by_distances(_experiments, CELLS_DISTANCES)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    # _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    # if BAND:
    #     _experiments = filtering.by_band(_experiments)

    # _experiments.remove(('SN41', 6, 'cells_1_2'))
    # _experiments.remove(('SN41', 2, 'cells_0_2'))

    _z_array = np.zeros(shape=(len(VALUES_BY_CELL_DIAMETER), len(VALUES_BY_CELL_DIAMETER)))
    for (_offset_y_index, _offset_y), (_offset_z_index, _offset_z) in \
            product(enumerate(VALUES_BY_CELL_DIAMETER), enumerate(VALUES_BY_CELL_DIAMETER)):

        # prepare data in mp
        _fibers_densities = {}
        _arguments = []
        _answers_keys = []
        for _tuple in _experiments:
            _experiment, _series, _group = _tuple
            _arguments.append((_experiment, _series, _group, ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH,
                               _OFFSET_X, _offset_y, _offset_z, DIRECTION))
            _answers_keys.append((_experiment, _series, _group))

        _p = Pool(CPUS_TO_USE)
        _answers = _p.starmap(compute.roi_fibers_density_by_time_pairs, _arguments)
        _p.close()

        # stop if needed
        if os.path.isfile(os.path.join(paths.EXPERIMENTS, 'stop.txt')):
            print('STOPPED!')
            return

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

            _master_properties = load.group_properties(_master_experiment, _master_series, _master_group)
            _master_left_cell_fibers_densities = compute.remove_blacklist(
                _master_experiment, _master_series, _master_properties['cells_ids']['left_cell'],
                _master_left_cell_fibers_densities)
            _master_right_cell_fibers_densities = compute.remove_blacklist(
                _master_experiment, _master_series, _master_properties['cells_ids']['right_cell'],
                _master_right_cell_fibers_densities)

            _master_left_cell_fibers_densities_filtered, _master_right_cell_fibers_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _master_left_cell_fibers_densities, _master_right_cell_fibers_densities
                )

            # ignore small arrays
            if _master_experiment in ['SN16', 'SN18']:
                if len(_master_left_cell_fibers_densities_filtered) < 15:
                    continue
            elif _master_experiment in ['SN41', 'SN44']:
                if len(_master_left_cell_fibers_densities_filtered) < 50:
                    continue
            else:
                raise Exception('No such experiment!')

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

                        _slave_properties = load.group_properties(_slave_experiment, _slave_series, _slave_group)
                        _master_fibers_densities = compute.remove_blacklist(
                            _master_experiment, _master_series, _master_properties['cells_ids'][_master_cell_id],
                            _master_fibers_densities)
                        _slave_fibers_densities = compute.remove_blacklist(
                            _slave_experiment, _slave_series, _slave_properties['cells_ids'][_slave_cell_id],
                            _slave_fibers_densities)

                        _master_fibers_densities_filtered, _slave_fibers_densities_filtered = \
                            compute.longest_same_indices_shared_in_borders_sub_array(
                                _master_fibers_densities, _slave_fibers_densities
                            )

                        # ignore small arrays
                        if _slave_experiment in ['SN16', 'SN18']:
                            if len(_master_fibers_densities_filtered) < 15:
                                continue
                        elif _slave_experiment in ['SN41', 'SN44']:
                            if len(_master_fibers_densities_filtered) < 50:
                                continue
                        else:
                            raise Exception('No such experiment!')

                        _slave_correlation = compute_lib.correlation(
                            compute_lib.derivative(_master_fibers_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_slave_fibers_densities_filtered, _n=DERIVATIVE)
                        )
                        _master_correlations_array.append(_master_correlation)
                        _slave_correlations_array.append(_slave_correlation)

        # compute percentages
        _master_minus_slave = np.array(_master_correlations_array) - np.array(_slave_correlations_array)
        _master_count = len(_master_minus_slave[_master_minus_slave > 0])
        if len(_master_minus_slave) > 0:
            _master_percentages = round(_master_count / len(_master_minus_slave), 10)
        else:
            _master_percentages = None

        print('z:', _offset_y, 'xy:', _offset_z, 'Master:', _master_percentages, 'N:', len(_master_minus_slave), 'Wilcoxon:', wilcoxon(_master_minus_slave), sep='\t')
        _z_array[_offset_z_index, _offset_y_index] = _master_percentages

    # plot
    _fig = heatmap.create_plot(
        _x_labels=VALUES_BY_CELL_DIAMETER,
        _y_labels=VALUES_BY_CELL_DIAMETER,
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
                  str(BAND) + '_smooth_heatmap'
    )

    _fig = contour.create_plot(
        _x_labels=VALUES_BY_CELL_DIAMETER,
        _y_labels=VALUES_BY_CELL_DIAMETER,
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
                  str(BAND) + '_smooth_contour'
    )


if __name__ == '__main__':
    main()

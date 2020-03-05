import os
import sys
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly
import seaborn as sns
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE, USE_MULTIPROCESSING
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT, NO_RETURN
from methods.experiments import export_video
from plotting import scatter, save, heatmap, contour

EXPERIMENTS = ['SN41']
EXPERIMENTS_STR = '_'.join(EXPERIMENTS)
REAL_CELLS = False
STATIC = True
BAND = False
DOMAIN_ABSOLUTE_VALUE = 5
VALUES_BY_CELL_DIAMETER = [
    round(_value, 10) for _value in
    np.arange(start=0, stop=DOMAIN_ABSOLUTE_VALUE * 2 + 0.1, step=0.1) - DOMAIN_ABSOLUTE_VALUE]
OFFSET_X = 0
DERIVATIVE = 2
CELLS_DISTANCES = [6, 7, 8, 9]
DIRECTION = 'inside'
MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50
}

# globals
_experiments = []
_arguments = None
_rois_dictionary = None
_rois_to_compute = None
_fibers_densities = None
_experiments_fibers_densities = None
_z_array = None
_annotations_array = None


def compute_data(_arguments):
    _offset_y_index, _offset_y, _offset_z_index, _offset_z = _arguments
    _master_correlations_array = []
    _slave_correlations_array = []
    for _master_index in range(len(_experiments)):
        _master_tuple = _experiments[_master_index]
        _master_experiment, _master_series, _master_group = _master_tuple

        _master_left_cell_fibers_densities = \
            _experiments_fibers_densities[
                (_master_experiment, _master_series, _master_group, _offset_y, _offset_z, 'left_cell')]
        _master_right_cell_fibers_densities = \
            _experiments_fibers_densities[
                (_master_experiment, _master_series, _master_group, _offset_y, _offset_z, 'right_cell')]

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
        if len(_master_left_cell_fibers_densities_filtered) < MINIMUM_CORRELATION_TIME_POINTS[_master_experiment]:
            continue

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
                        _experiments_fibers_densities[
                            (_master_experiment, _master_series, _master_group, _offset_y, _offset_z,
                             _master_cell_id)]
                    _slave_fibers_densities = \
                        _experiments_fibers_densities[
                            (_slave_experiment, _slave_series, _slave_group, _offset_y, _offset_z, _slave_cell_id)]

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
                    if len(_master_fibers_densities_filtered) < MINIMUM_CORRELATION_TIME_POINTS[_slave_experiment]:
                        continue

                    _slave_correlation = compute_lib.correlation(
                        compute_lib.derivative(_master_fibers_densities_filtered, _n=DERIVATIVE),
                        compute_lib.derivative(_slave_fibers_densities_filtered, _n=DERIVATIVE)
                    )
                    _master_correlations_array.append(_master_correlation)
                    _slave_correlations_array.append(_slave_correlation)

    # compute percentages
    _annotation = None
    _master_minus_slave = np.array(_master_correlations_array) - np.array(_slave_correlations_array)
    _master_count = len(_master_minus_slave[_master_minus_slave > 0])
    if len(_master_minus_slave) > 0:
        _master_percentages = round(_master_count / len(_master_minus_slave), 10)
        _wilcoxon = wilcoxon(_master_minus_slave)
        _p_value = _wilcoxon[1]
        if _p_value > 0.05:
            _annotation = {
                'text': 'X',
                'show_arrow': False,
                'x': _offset_y,
                'y': _offset_z,
                'font': {
                    'size': 6,
                    'color': 'red'
                }
            }
    else:
        _master_percentages = None

    return _offset_z_index, _offset_y_index, _master_percentages, _annotation


def main():
    global _experiments, _arguments, _rois_dictionary, _rois_to_compute, _fibers_densities,\
        _experiments_fibers_densities, _z_array, _annotations_array

    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_distances(_experiments, CELLS_DISTANCES)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    if BAND:
        _experiments = filtering.by_band(_experiments)

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _offset_y, _offset_z, _cell_id in \
                product(VALUES_BY_CELL_DIAMETER, VALUES_BY_CELL_DIAMETER, ['left_cell', 'right_cell']):
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': ROI_LENGTH,
                'length_y': ROI_HEIGHT,
                'length_z': ROI_WIDTH,
                'offset_x': OFFSET_X,
                'offset_y': _offset_y,
                'offset_z': _offset_z,
                'cell_id': _cell_id,
                'direction': DIRECTION
            })

    _rois_dictionary, _rois_to_compute =\
        compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'offset_y', 'offset_z', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {}
    for _key in tqdm(_rois_dictionary, desc='Organizing Fibers Densities'):
        _experiments_fibers_densities[_key] = [_fibers_densities[_tuple] for _tuple in _rois_dictionary[_key]]

    _arguments = []
    for (_offset_y_index, _offset_y), (_offset_z_index, _offset_z) in \
            product(enumerate(VALUES_BY_CELL_DIAMETER), enumerate(VALUES_BY_CELL_DIAMETER)):
        _arguments.append(
            (_offset_y_index, _offset_y, _offset_z_index, _offset_z))

    _z_array = np.zeros(shape=(len(VALUES_BY_CELL_DIAMETER), len(VALUES_BY_CELL_DIAMETER)))
    _annotations_array = []
    with Pool(CPUS_TO_USE) as _p:
        for _answer in tqdm(_p.imap_unordered(compute_data, _arguments), total=len(_arguments),
                            desc='Computing Heatmap'):
            _offset_z_index, _offset_y_index, _master_percentages, _annotation = _answer
            _z_array[_offset_z_index, _offset_y_index] = _master_percentages
            if _annotation is not None:
                _annotations_array.append(_annotation)
        _p.close()
        _p.join()

    # plot
    _fig = heatmap.create_plot(
        _x_labels=VALUES_BY_CELL_DIAMETER,
        _y_labels=VALUES_BY_CELL_DIAMETER,
        _z_array=_z_array,
        _x_axis_title='Offset in Z axis',
        _y_axis_title='Offset in XY axis',
        _annotations_array=_annotations_array,
        _color_scale=sns.color_palette('BrBG').as_hex(),
        _zmin=0,
        _zmax=1
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_' + EXPERIMENTS_STR + '_real_' + str(REAL_CELLS) + '_static_' + str(STATIC) + '_band_' +
                  str(BAND) + '_heatmap'
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
                  str(BAND) + '_contour'
    )


if __name__ == '__main__':
    main()

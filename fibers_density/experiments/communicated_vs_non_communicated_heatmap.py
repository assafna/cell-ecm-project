import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import seaborn as sns
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import save, heatmap, contour

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
_experiments_fibers_densities = None
_z_array = None
_annotations_array = None


def compute_data(_arguments):
    _offset_y_index, _offset_y, _offset_z_index, _offset_z = _arguments
    _communicated_correlations_array = []
    _non_communicated_correlations_array = []
    for _communicated_index in range(len(_experiments)):
        _communicated_tuple = _experiments[_communicated_index]
        _communicated_experiment, _communicated_series, _communicated_group = _communicated_tuple

        _communicated_left_cell_fibers_densities = _experiments_fibers_densities[(
            _communicated_experiment,
            _communicated_series,
            _communicated_group,
            _offset_y,
            _offset_z,
            'left_cell'
        )]
        _communicated_right_cell_fibers_densities = _experiments_fibers_densities[(
            _communicated_experiment,
            _communicated_series,
            _communicated_group,
            _offset_y,
            _offset_z,
            'right_cell'
        )]

        _communicated_properties = \
            load.group_properties(_communicated_experiment, _communicated_series, _communicated_group)
        _communicated_left_cell_fibers_densities = compute.remove_blacklist(
            _communicated_experiment, _communicated_series, _communicated_properties['cells_ids']['left_cell'],
            _communicated_left_cell_fibers_densities)
        _communicated_right_cell_fibers_densities = compute.remove_blacklist(
            _communicated_experiment, _communicated_series, _communicated_properties['cells_ids']['right_cell'],
            _communicated_right_cell_fibers_densities)

        _communicated_left_cell_fibers_densities_filtered, _communicated_right_cell_fibers_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _communicated_left_cell_fibers_densities, _communicated_right_cell_fibers_densities
            )

        # ignore small arrays
        if len(_communicated_left_cell_fibers_densities_filtered) < \
                MINIMUM_CORRELATION_TIME_POINTS[_communicated_experiment]:
            continue

        _communicated_correlation = compute_lib.correlation(
            compute_lib.derivative(_communicated_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
            compute_lib.derivative(_communicated_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
        )
        for _non_communicated_index in range(len(_experiments)):
            if _communicated_index != _non_communicated_index:
                _non_communicated_tuple = _experiments[_non_communicated_index]
                _non_communicated_experiment, _non_communicated_series, _non_communicated_group = \
                    _non_communicated_tuple

                for _communicated_cell_id, _non_communicated_cell_id in product(['left_cell', 'right_cell'],
                                                                                ['left_cell', 'right_cell']):
                    _communicated_fibers_densities = _experiments_fibers_densities[(
                        _communicated_experiment,
                        _communicated_series,
                        _communicated_group,
                        _offset_y,
                        _offset_z,
                        _communicated_cell_id
                    )]
                    _non_communicated_fibers_densities = _experiments_fibers_densities[(
                        _non_communicated_experiment,
                        _non_communicated_series,
                        _non_communicated_group,
                        _offset_y,
                        _offset_z,
                        _non_communicated_cell_id
                    )]

                    _non_communicated_properties = load.group_properties(
                        _non_communicated_experiment,
                        _non_communicated_series,
                        _non_communicated_group
                    )
                    _communicated_fibers_densities = compute.remove_blacklist(
                        _communicated_experiment,
                        _communicated_series,
                        _communicated_properties['cells_ids'][_communicated_cell_id],
                        _communicated_fibers_densities
                    )
                    _non_communicated_fibers_densities = compute.remove_blacklist(
                        _non_communicated_experiment,
                        _non_communicated_series,
                        _non_communicated_properties['cells_ids'][_non_communicated_cell_id],
                        _non_communicated_fibers_densities
                    )

                    _communicated_fibers_densities_filtered, _non_communicated_fibers_densities_filtered = \
                        compute.longest_same_indices_shared_in_borders_sub_array(
                            _communicated_fibers_densities, _non_communicated_fibers_densities
                        )

                    # ignore small arrays
                    if len(_communicated_fibers_densities_filtered) < \
                            MINIMUM_CORRELATION_TIME_POINTS[_non_communicated_experiment]:
                        continue

                    _non_communicated_correlation = compute_lib.correlation(
                        compute_lib.derivative(_communicated_fibers_densities_filtered, _n=DERIVATIVE),
                        compute_lib.derivative(_non_communicated_fibers_densities_filtered, _n=DERIVATIVE)
                    )
                    _communicated_correlations_array.append(_communicated_correlation)
                    _non_communicated_correlations_array.append(_non_communicated_correlation)

    # compute percentages
    _annotation = None
    _communicated_minus_non_communicated = \
        np.array(_communicated_correlations_array) - np.array(_non_communicated_correlations_array)
    _communicated_count = len(_communicated_minus_non_communicated[_communicated_minus_non_communicated > 0])
    if len(_communicated_minus_non_communicated) > 0:
        _communicated_percentages = round(_communicated_count / len(_communicated_minus_non_communicated), 10)
        _wilcoxon = wilcoxon(_communicated_minus_non_communicated)
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
        _communicated_percentages = None

    return _offset_z_index, _offset_y_index, _communicated_percentages, _annotation


def main():
    global _experiments, _experiments_fibers_densities, _z_array, _annotations_array

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

    # clean
    _fibers_densities = None
    _rois_dictionary = None
    _rois_to_compute = None

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
            _offset_z_index, _offset_y_index, _communicated_percentages, _annotation = _answer
            _z_array[_offset_z_index, _offset_y_index] = _communicated_percentages
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

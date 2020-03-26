import os

import numpy as np
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, config, filtering, compute, organize, paths
from plotting import box, save

OFFSET_X = 0
OFFSET_Y = 0

# experiments
TIME_POINTS = 18
OFFSET_Z = 0
OUT_OF_BOUNDARIES = False


def compute_single_cell_mean(_experiment, _series_id, _cell_tuples, _rois_dictionary, _fibers_densities):
    _normalization = load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))
    _cell_fibers_densities = []
    for _time_point in range(TIME_POINTS):
        _time_point_fibers_densities = []
        for _cell_tuple in _cell_tuples:
            _, _, _group = _cell_tuple
            for _direction in ['left', 'right']:
                _roi_tuple = _rois_dictionary[(_experiment, _series_id, _group, _direction)][_time_point]
                _fibers_density = _fibers_densities[_roi_tuple]

                if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                    continue

                _normalized_fibers_density = compute_lib.z_score(
                    _x=_fibers_density[0],
                    _average=_normalization['average'],
                    _std=_normalization['std']
                )

                if not np.isnan(_normalized_fibers_density):
                    _time_point_fibers_densities.append(_normalized_fibers_density)

        _cell_fibers_densities.append(np.mean(_time_point_fibers_densities))

    return _cell_fibers_densities


def main():
    _experiments = load.experiments_groups_as_tuples(config.SINGLE_CELL)
    _experiments = filtering.by_time_points_amount(_experiments, TIME_POINTS)
    _experiments = filtering.by_main_cell(_experiments)

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _direction in ['left', 'right']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': config.ROI_LENGTH,
                'length_y': config.ROI_HEIGHT,
                'length_z': config.ROI_WIDTH,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': 'cell',
                'direction': _direction,
                'time_points': TIME_POINTS
            })

    _rois_dictionary, _rois_to_compute = \
        compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'direction'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments = organize.by_single_cell_id(_experiments)
    _experiments_ids = list(_experiments.keys())

    _correlations_array_derivative_0 = []
    _correlations_array_derivative_1 = []
    _correlations_array_derivative_2 = []
    for _index_1 in tqdm(range(len(_experiments_ids)), desc='Main Loop'):
        _tuple_1 = _experiments_ids[_index_1]
        _experiment_1, _series_id_1, _cell_id_1 = _tuple_1
        _fibers_densities_1 = compute_single_cell_mean(
            _experiment=_experiment_1,
            _series_id=_series_id_1,
            _cell_tuples=_experiments[_tuple_1],
            _rois_dictionary=_rois_dictionary,
            _fibers_densities=_fibers_densities
        )
        for _index_2 in range(_index_1 + 1, len(_experiments_ids)):
            _tuple_2 = _experiments_ids[_index_2]
            _experiment_2, _series_id_2, _cell_id_2 = _tuple_2
            _fibers_densities_2 = compute_single_cell_mean(
                _experiment=_experiment_2,
                _series_id=_series_id_2,
                _cell_tuples=_experiments[_tuple_2],
                _rois_dictionary=_rois_dictionary,
                _fibers_densities=_fibers_densities
            )
            _correlations_array_derivative_0.append(compute_lib.correlation(
                compute_lib.derivative(_fibers_densities_1, _n=0),
                compute_lib.derivative(_fibers_densities_2, _n=0)
            ))
            _correlations_array_derivative_1.append(compute_lib.correlation(
                compute_lib.derivative(_fibers_densities_1, _n=1),
                compute_lib.derivative(_fibers_densities_2, _n=1)
            ))
            _correlations_array_derivative_2.append(compute_lib.correlation(
                compute_lib.derivative(_fibers_densities_1, _n=2),
                compute_lib.derivative(_fibers_densities_2, _n=2)
            ))

    # plot
    _fig = box.create_plot(
        _y_array=[_correlations_array_derivative_0, _correlations_array_derivative_1,
                  _correlations_array_derivative_2],
        _names_array=['Fibers Densities', '1st Derivative', '2nd Derivative'],
        _x_axis_title='Fibers Densities by Derivatives',
        _y_axis_title='Correlation'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )

    # wilcoxon
    print('Derivative 0:', wilcoxon(_correlations_array_derivative_0))
    print('Derivative 1:', wilcoxon(_correlations_array_derivative_1))
    print('Derivative 2:', wilcoxon(_correlations_array_derivative_2))


if __name__ == '__main__':
    main()

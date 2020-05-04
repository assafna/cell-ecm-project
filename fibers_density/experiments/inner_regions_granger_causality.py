from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import grangercausalitytests, adfuller

from libs import compute_lib
from libs.experiments import load, filtering, compute
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT

# based on time resolution
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41', 'SN44']
}
TIME_RESOLUTION = {
    False: 15,
    True: 5
}
OFFSET_X = 0
# TODO: set the offset in y according to the angle in the original Z slices of the cells
OFFSET_Y = 0.5
OFFSET_Z = 0
CELLS_DISTANCE_RANGE = [4, 10]
REAL_CELLS = True
STATIC = False
DIRECTION = 'inside'
MINIMUM_TIME_POINTS = 30


def main(_band=True, _high_time_resolution=True):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_time_points_amount(_experiments, MINIMUM_TIME_POINTS)
    _experiments = filtering.by_distance_range(_experiments, CELLS_DISTANCE_RANGE)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    _experiments = filtering.by_band(_experiments, _band=_band)
    print('Total experiments:', len(_experiments))

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple

        # stop when windows are overlapping
        _properties = load.group_properties(_experiment, _series_id, _group)
        _latest_time_point = len(_properties['time_points'])
        if DIRECTION == 'inside':
            for _time_point in range(len(_properties['time_points'])):
                _cells_distance = \
                    compute.cells_distance_in_cell_size_time_point(_experiment, _series_id, _group, _time_point)
                if _cells_distance - 1 - OFFSET_X * 2 < ROI_LENGTH * 2:
                    _latest_time_point = _time_point - 1
                    break

        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': ROI_LENGTH,
                'length_y': ROI_HEIGHT,
                'length_z': ROI_WIDTH,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': DIRECTION,
                'time_points': _latest_time_point
            })

    _rois_dictionary, _rois_to_compute = compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {
        _key: [_fibers_densities[_tuple] for _tuple in _rois_dictionary[_key]]
        for _key in _rois_dictionary
    }

    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple

        _left_cell_fibers_densities = _experiments_fibers_densities[(_experiment, _series_id, _group, 'left_cell')]
        _right_cell_fibers_densities = _experiments_fibers_densities[(_experiment, _series_id, _group, 'right_cell')]

        _properties = load.group_properties(_experiment, _series_id, _group)
        _left_cell_fibers_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_cell_fibers_densities)
        _right_cell_fibers_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_cell_fibers_densities)

        _left_cell_fibers_densities_filtered, _right_cell_fibers_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _left_cell_fibers_densities, _right_cell_fibers_densities)

        # ignore small arrays
        if len(_left_cell_fibers_densities_filtered) < MINIMUM_TIME_POINTS:
            continue

        # stationary test
        _left_cell_fibers_densities_derivative = None
        _right_cell_fibers_densities_derivative = None
        _derivative = None
        for _derivative in range(10):
            _left_cell_fibers_densities_derivative = \
                compute_lib.derivative(_left_cell_fibers_densities_filtered, _n=_derivative)
            _, _adf_p_value, _, _, _, _ = adfuller(_left_cell_fibers_densities_derivative)
            if _adf_p_value > 0.05:
                continue
            _right_cell_fibers_densities_derivative = \
                compute_lib.derivative(_right_cell_fibers_densities_filtered, _n=_derivative)
            _, _adf_p_value, _, _, _, _ = adfuller(_right_cell_fibers_densities_derivative)

            # both are stationary
            if _adf_p_value < 0.05:
                break

        # causality
        try:
            for _causality, _cell_1_array, _cell_2_array in zip(
                    [['left', 'right'], ['right', 'left']],
                    [_left_cell_fibers_densities_derivative, _right_cell_fibers_densities_derivative],
                    [_right_cell_fibers_densities_derivative, _left_cell_fibers_densities_derivative]
            ):
                _x = [[_cell_1_value, _cell_2_value] for _cell_1_value, _cell_2_value in
                      zip(_cell_1_array, _cell_2_array)]

                # var model to retrieve lag
                _var_model = VAR(_x)
                _lag_order_results = _var_model.select_order()
                _estimators_lags = [
                    _lag_order_results.aic,
                    _lag_order_results.bic,
                    _lag_order_results.fpe,
                    _lag_order_results.hqic
                ]
                _min_estimator_lag = min(_estimators_lags)

                # found a lag
                if _min_estimator_lag > 0:

                    # granger causality
                    _granger_causality_results = grangercausalitytests(x=_x, maxlag=_min_estimator_lag, verbose=False)
                    _f_test_p_value = _granger_causality_results[_min_estimator_lag][0]['ssr_ftest'][1]

                    if _f_test_p_value < 0.05:
                        print(_tuple, _causality[0].capitalize() + ' causality ' + _causality[1] + '!',
                              'time-points:' + str(len(_left_cell_fibers_densities_derivative)),
                              'stationary derivative: ' + str(_derivative),
                              'p-value: ' + str(round(_f_test_p_value, 4)),
                              'lag: ' + str(_min_estimator_lag), sep='\t')

        # not enough time points
        except ValueError:
            continue


if __name__ == '__main__':
    main()

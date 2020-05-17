import warnings

import pandas as pd
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.stats.stattools import durbin_watson
from statsmodels.tools.sm_exceptions import InterpolationWarning
from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import grangercausalitytests, adfuller, kpss, acf

from libs import compute_lib
from libs.experiments import load, filtering, compute
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT

# based on time resolution
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41', 'SN44', 'SN45']
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
MAXIMUM_LAG = 10

# stationary tests
ADF_TEST = True
KPSS_TEST = True


def main(_band=None, _high_time_resolution=True):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_time_points_amount(_experiments, _time_points=MINIMUM_TIME_POINTS)
    _experiments = filtering.by_distance_range(_experiments, _distance_range=CELLS_DISTANCE_RANGE)
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

        _left_cell_fibers_densities = \
            _experiments_fibers_densities[(_experiment, _series_id, _group, 'left_cell')]
        _right_cell_fibers_densities = \
            _experiments_fibers_densities[(_experiment, _series_id, _group, 'right_cell')]

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
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=InterpolationWarning)

            # find derivative for stationary
            for _derivative in range(10):
                _left_cell_fibers_densities_derivative = \
                    compute_lib.derivative(_left_cell_fibers_densities_filtered, _n=_derivative)
                _right_cell_fibers_densities_derivative = \
                    compute_lib.derivative(_right_cell_fibers_densities_filtered, _n=_derivative)

                if ADF_TEST:
                    _, _left_cell_adf_p_value, _, _, _, _ = adfuller(_left_cell_fibers_densities_derivative)
                    _, _right_cell_adf_p_value, _, _, _, _ = adfuller(_right_cell_fibers_densities_derivative)
                    if _left_cell_adf_p_value > 0.05 or _right_cell_adf_p_value > 0.05:
                        continue

                if KPSS_TEST:
                    _, _left_cell_kpss_p_value, _, _ = kpss(_left_cell_fibers_densities_derivative, nlags='legacy')
                    _, _right_cell_kpss_p_value, _, _ = kpss(_right_cell_fibers_densities_derivative, nlags='legacy')
                    if _left_cell_kpss_p_value < 0.05 or _right_cell_kpss_p_value < 0.05:
                        continue

                # stationary
                break

        # causality
        try:
            _x = pd.DataFrame(
                data=[[_left_value, _right_value] for _left_value, _right_value in
                      zip(_left_cell_fibers_densities_derivative, _right_cell_fibers_densities_derivative)],
                columns=['left', 'right'])

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
            if 0 < _min_estimator_lag <= MAXIMUM_LAG:
                _var_model_results = _var_model.fit(maxlags=_min_estimator_lag, ic=None)

                # _acf_result_1 = \
                #     acf(x=_var_model_results.resid[:, 0], fft=True, nlags=_min_estimator_lag, qstat=True)
                # _acf_result_2 = \
                #     acf(x=_var_model_results.resid[:, 1], fft=True, nlags=_min_estimator_lag, qstat=True)
                # if _acf_result_1[2][0] > 0.05 and _acf_result_2[2][0] > 0.05:
                #     _residuals_autocorrelation = False
                _whiteness = _var_model_results.test_whiteness(nlags=_min_estimator_lag + 1)

                # no autocorrelation in the residuals
                if _whiteness.pvalue > 0.05:

                    # granger causality
                    for _caused, _causing in zip(['left', 'right'], ['right', 'left']):
                        _granger = _var_model_results.test_causality(caused=_caused, causing=_causing)

                        if _granger.pvalue < 0.05:
                            _normality = _var_model_results.test_normality()
                            _inst_granger = _var_model_results.test_inst_causality (causing=_causing)

                            print(_tuple, _causing.capitalize() + ' causes ' + _caused + '!',
                                  'time-points: ' + str(len(_left_cell_fibers_densities_derivative)),
                                  'stationary derivative: ' + str(_derivative),
                                  'p-value: ' + str(round(_granger.pvalue, 4)),
                                  'lag: ' + str(_min_estimator_lag),
                                  'normality p-value: ' + str(round(_normality.pvalue, 4)),
                                  'inst p-value: ' + str(round(_inst_granger.pvalue, 4)),
                                  sep='\t')

        # not enough time points
        except ValueError:
            continue


if __name__ == '__main__':
    main()

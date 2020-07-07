import os
import warnings

import numpy as np
import pandas as pd
import plotly.graph_objs as go
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests
from statsmodels.tools.sm_exceptions import InterpolationWarning
from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import adfuller, kpss

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import save

# based on time resolution
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41', 'SN44', 'SN45']
}
TIME_POINT = {
    False: 18,
    True: 52
}
OFFSET_X = 0
# TODO: set the offset in y according to the angle in the original Z slices of the cells
OFFSET_Y = 0.5
OFFSET_Z = 0
TIME_RESOLUTION = 5
CELLS_DISTANCE_RANGE = [4, 10]
REAL_CELLS = True
STATIC = False
DIRECTION = 'inside'
MINIMUM_TIME_POINTS = 30
MAXIMUM_LAG = 10

# stationary tests
ADF_TEST = True
KPSS_TEST = True


def main(_band=None, _high_time_resolution=True, _tuples_to_mark=None, _tuples_to_plot=None, _plots=None):
    if _plots is None:
        _plots = ['whiteness', 'granger']

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

    _n_pairs = 0
    _n_pairs_with_band = 0
    _whiteness_p_values = []
    _n_passed_whiteness_with_band = 0
    _granger_causality_p_values = []
    _n_passed_granger_causality_with_band = 0
    _correlations = []
    _time_lag_correlations = []
    _end_fiber_densities = []
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

        _n_pairs += 1
        if _properties['band']:
            _n_pairs_with_band += 1

        _start_time_point = 0
        for _left in _left_cell_fibers_densities:
            if _left[0] == _left_cell_fibers_densities_filtered[0]:
                break
            _start_time_point += 1

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

                _whiteness = _var_model_results.test_whiteness(nlags=_min_estimator_lag + 1)
                _whiteness_p_values.append(_whiteness.pvalue)

                if _tuples_to_mark is not None and _tuple in _tuples_to_mark and _whiteness.pvalue > 0.05:
                    print(_tuple, 'marked whiteness p-value:', _whiteness.pvalue)

                # no autocorrelation in the residuals
                if _whiteness.pvalue > 0.05:
                    if _properties['band']:
                        _n_passed_whiteness_with_band += 1

                    # granger causality
                    for _caused, _causing in zip(['left', 'right'], ['right', 'left']):
                        _granger = _var_model_results.test_causality(caused=_caused, causing=_causing)
                        _granger_causality_p_values.append(_granger.pvalue)

                        # time lag = 0
                        _correlation = compute_lib.correlation(
                            _left_cell_fibers_densities_derivative, _right_cell_fibers_densities_derivative)
                        _correlations.append(_correlation)

                        # time lag = min estimator
                        if _causing == 'left':
                            _left_fibers_densities_time_lag = \
                                _left_cell_fibers_densities_derivative[:-_min_estimator_lag]
                            _right_fibers_densities_time_lag = \
                                _right_cell_fibers_densities_derivative[_min_estimator_lag:]
                        else:
                            _left_fibers_densities_time_lag = \
                                _left_cell_fibers_densities_derivative[_min_estimator_lag:]
                            _right_fibers_densities_time_lag = \
                                _right_cell_fibers_densities_derivative[:-_min_estimator_lag]
                        _time_lag_correlation = compute_lib.correlation(
                            _left_fibers_densities_time_lag, _right_fibers_densities_time_lag
                        )
                        _time_lag_correlations.append(_time_lag_correlation)

                        # end fiber density
                        if len(_left_cell_fibers_densities_filtered) > TIME_POINT[_high_time_resolution]:
                            _end_fiber_density = \
                                (_left_cell_fibers_densities_filtered[TIME_POINT[_high_time_resolution]] +
                                 _right_cell_fibers_densities_filtered[TIME_POINT[_high_time_resolution]]) / 2
                        else:
                            _end_fiber_density = \
                                (_left_cell_fibers_densities_filtered[-1] +
                                 _right_cell_fibers_densities_filtered[-1]) / 2
                        _end_fiber_densities.append(_end_fiber_density)

                        # marking
                        if _tuples_to_mark is not None and _tuple in _tuples_to_mark and _granger.pvalue < 0.05:
                            print(_tuple, 'causing:', _causing, 'marked granger p-value:', _granger.pvalue)

                        if _granger.pvalue < 0.05:
                            if _properties['band']:
                                _n_passed_granger_causality_with_band += 1

                            _normality = _var_model_results.test_normality()
                            _inst_granger = _var_model_results.test_inst_causality(causing=_causing)

                            print(_tuple, _causing.capitalize() + ' causes ' + _caused + '!',
                                  'time-points: ' + str(len(_left_cell_fibers_densities_derivative)),
                                  'stationary derivative: ' + str(_derivative),
                                  'band:' + str(_properties['band']),
                                  'p-value: ' + str(round(_granger.pvalue, 4)),
                                  'lag: ' + str(_min_estimator_lag),
                                  'normality p-value: ' + str(round(_normality.pvalue, 4)),
                                  'inst p-value: ' + str(round(_inst_granger.pvalue, 4)),
                                  sep='\t')

                            # lag = 0
                            print('Time lag = 0 correlation:', _correlation)

                            # rest of lags
                            for _lag in range(1, _min_estimator_lag + 1):
                                if _causing == 'left':
                                    _left_fibers_densities_time_lag = _left_cell_fibers_densities_derivative[:-_lag]
                                    _right_fibers_densities_time_lag = _right_cell_fibers_densities_derivative[_lag:]
                                else:
                                    _left_fibers_densities_time_lag = _left_cell_fibers_densities_derivative[_lag:]
                                    _right_fibers_densities_time_lag = _right_cell_fibers_densities_derivative[:-_lag]

                                _correlation = compute_lib.correlation(
                                    _left_fibers_densities_time_lag, _right_fibers_densities_time_lag
                                )
                                print('Time lag = ' + str(_lag) + ' correlation:', _correlation)

                            # plots
                            if _tuples_to_plot is not None and _tuple in _tuples_to_plot:
                                _y_arrays = [_left_cell_fibers_densities_derivative,
                                             _right_cell_fibers_densities_derivative]
                                _names_array = ['Left cell', 'Right cell']
                                _colors_array = ['#844b00', '#ea8500']
                                _fig = go.Figure(
                                    data=[
                                        go.Scatter(
                                            x=np.arange(
                                                start=_start_time_point,
                                                stop=_start_time_point + len(_left_cell_fibers_densities_derivative),
                                                step=1) * TIME_RESOLUTION,
                                            y=_y,
                                            name=_name,
                                            mode='lines',
                                            line={
                                                'color': _color,
                                                'width': 1
                                            }
                                        ) for _y, _name, _color in zip(_y_arrays, _names_array, _colors_array)
                                    ],
                                    layout={
                                        'xaxis': {
                                            'title': 'Time (minutes)',
                                            'zeroline': False
                                        },
                                        'yaxis': {
                                            'title': 'Fiber density (z-score)' + '\'' * _derivative,
                                            'zeroline': False
                                        },
                                        'legend': {
                                            'xanchor': 'left',
                                            'x': 0.1,
                                            'yanchor': 'top',
                                            'bordercolor': 'black',
                                            'borderwidth': 2,
                                            'bgcolor': 'white'
                                        },
                                    }
                                )

                                _experiment, _series_id, _group = _tuple
                                save.to_html(
                                    _fig=_fig,
                                    _path=os.path.join(paths.PLOTS, save.get_module_name()),
                                    _filename='plot_' + _experiment + '_' + str(_series_id) + '_' + _group
                                )

                                # residuals
                                _y_arrays = \
                                    [_var_model_results.resid.values[:, 0], _var_model_results.resid.values[:, 1]]
                                _fig = go.Figure(
                                    data=[
                                        go.Scatter(
                                            x=np.arange(
                                                start=_start_time_point,
                                                stop=_start_time_point + len(_y),
                                                step=1) * TIME_RESOLUTION,
                                            y=_y,
                                            name=_name,
                                            mode='lines',
                                            line={
                                                'color': _color,
                                                'width': 1
                                            }
                                        ) for _y, _name, _color in zip(_y_arrays, _names_array, _colors_array)
                                    ],
                                    layout={
                                        'xaxis': {
                                            'title': 'Time (minutes)',
                                            'zeroline': False
                                        },
                                        'yaxis': {
                                            'title': 'Residual',
                                            'zeroline': False
                                        },
                                        'legend': {
                                            'xanchor': 'left',
                                            'x': 0.1,
                                            'yanchor': 'top',
                                            'bordercolor': 'black',
                                            'borderwidth': 2,
                                            'bgcolor': 'white'
                                        },
                                    }
                                )

                                _experiment, _series_id, _group = _tuple
                                save.to_html(
                                    _fig=_fig,
                                    _path=os.path.join(paths.PLOTS, save.get_module_name()),
                                    _filename='plot_residuals_' + _experiment + '_' + str(_series_id) + '_' + _group
                                )

        # not enough time points
        except ValueError:
            continue

    print('Total pairs:', _n_pairs)
    print('Total pairs with band:', _n_pairs_with_band)
    print('Total pairs passed whiteness:', (np.array(_whiteness_p_values) > 0.05).sum())
    print('Total pairs passed whiteness with band:', _n_passed_whiteness_with_band)
    print('Total cells passed granger causality:', (np.array(_granger_causality_p_values) < 0.05).sum())
    print('Total cells passed granger causality with band:', _n_passed_granger_causality_with_band)

    # p-value correction
    print('Corrections of GC p-value < 0.05:')
    _granger_causality_p_values_corrected = multipletests(pvals=_granger_causality_p_values, method='fdr_bh')
    for _p_value, _p_value_corrected in zip(_granger_causality_p_values, _granger_causality_p_values_corrected[1]):
        if _p_value < 0.05:
            print('Original GC p-value:', _p_value, 'corrected:', _p_value_corrected)

    # plots
    for _test_name, _y_title, _y_array in \
            zip(
                ['whiteness', 'granger'],
                ['Whiteness p-value', 'Granger causality p-value'],
                [_whiteness_p_values, _granger_causality_p_values]
            ):
        if _test_name in _plots:
            _fig = go.Figure(
                data=go.Box(
                    y=_y_array,
                    boxpoints='all',
                    jitter=1,
                    pointpos=0,
                    line={
                        'width': 1
                    },
                    fillcolor='white',
                    marker={
                        'size': 10,
                        'color': '#ea8500'
                    },
                    opacity=0.7,
                    showlegend=False
                ),
                layout={
                    'xaxis': {
                        'zeroline': False
                    },
                    'yaxis': {
                        'title': _y_title,
                        'zeroline': False,
                        'range': [-0.1, 1.1],
                        'tickmode': 'array',
                        'tickvals': [0.05, 1]
                    },
                    'shapes': [
                        {
                            'type': 'line',
                            'x0': -0.75,
                            'y0': 0.05,
                            'x1': 0.75,
                            'y1': 0.05,
                            'line': {
                                'color': 'red',
                                'width': 2,
                                'dash': 'dash'
                            }
                        }
                    ]
                }
            )

            save.to_html(
                _fig=_fig,
                _path=os.path.join(paths.PLOTS, save.get_module_name()),
                _filename='plot_' + _test_name
            )

    # granger versus correlation
    print('GC vs. correlation pearson correlation:', pearsonr(_granger_causality_p_values, _correlations))
    _fig = go.Figure(
        data=go.Scatter(
            x=_granger_causality_p_values,
            y=_correlations,
            mode='markers',
            marker={
                'size': 10,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Granger causality p-value',
                'zeroline': False,
            },
            'yaxis': {
                'title': 'Inner correlation',
                'zeroline': False,
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_gc_vs_correlation'
    )

    # granger versus time lag correlation
    print('GC vs. time lag correlation pearson correlation:', pearsonr(_granger_causality_p_values,
                                                                       _time_lag_correlations))
    _fig = go.Figure(
        data=go.Scatter(
            x=_granger_causality_p_values,
            y=_time_lag_correlations,
            mode='markers',
            marker={
                'size': 10,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Granger causality p-value',
                'zeroline': False,
            },
            'yaxis': {
                'title': 'GC lag inner correlation',
                'zeroline': False,
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_gc_vs_time_lag_correlation'
    )

    # granger versus end fiber density
    print('GC vs. end fiber density pearson correlation:', pearsonr(_granger_causality_p_values, _end_fiber_densities))
    _fig = go.Figure(
        data=go.Scatter(
            x=_granger_causality_p_values,
            y=_end_fiber_densities,
            mode='markers',
            marker={
                'size': 10,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Granger causality p-value',
                'zeroline': False,
            },
            'yaxis': {
                'title': 'End fiber density (z-score)',
                'zeroline': False,
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_gc_vs_end_density'
    )


if __name__ == '__main__':
    main()

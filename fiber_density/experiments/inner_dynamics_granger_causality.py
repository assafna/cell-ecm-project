import os
import warnings

import numpy as np
import pandas as pd
import plotly.graph_objs as go
from statsmodels.stats.multitest import multipletests
from statsmodels.tools.sm_exceptions import InterpolationWarning
from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import adfuller, kpss

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths, config
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0.5
OFFSET_Z = 0

PAIR_DISTANCE_RANGE = [4, 10]
MINIMUM_TIME_FRAMES = 30
MAXIMUM_LAG = 10

# stationary tests
ADF_TEST = True
KPSS_TEST = True


def main(_band=None, _high_temporal_resolution=True, _tuples_to_mark=None, _tuples_to_plot=None, _plots=None):
    if _plots is None:
        _plots = ['whiteness', 'granger']

    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=_high_temporal_resolution,
        _is_bleb=False,
        _is_bleb_from_start=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_time_frames_amount(_tuples, _time_frames=MINIMUM_TIME_FRAMES)
    _tuples = filtering.by_pair_distance_range(_tuples, _distance_range=PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples, _band=_band)
    print('Total tuples:', len(_tuples))

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple

        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'length_z': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, OFFSET_X)
            })

    _windows_dictionary, _windows_to_compute = compute.windows(_arguments,
                                                               _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
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
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple

        _left_cell_fiber_densities = \
            _experiments_fiber_densities[(_experiment, _series_id, _group, 'left_cell')]
        _right_cell_fiber_densities = \
            _experiments_fiber_densities[(_experiment, _series_id, _group, 'right_cell')]

        _properties = load.group_properties(_experiment, _series_id, _group)
        _left_cell_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_cell_fiber_densities)
        _right_cell_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_cell_fiber_densities)

        _left_cell_fiber_densities_filtered, _right_cell_fiber_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _left_cell_fiber_densities, _right_cell_fiber_densities)

        # ignore small arrays
        if len(_left_cell_fiber_densities_filtered) < MINIMUM_TIME_FRAMES:
            continue

        _n_pairs += 1
        if _properties['band']:
            _n_pairs_with_band += 1

        _start_time_frame = 0
        for _left in _left_cell_fiber_densities:
            if _left[0] == _left_cell_fiber_densities_filtered[0]:
                break
            _start_time_frame += 1

        # stationary test
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=InterpolationWarning)

            # find derivative for stationary
            for _derivative in range(10):
                _left_cell_fiber_densities_derivative = \
                    compute_lib.derivative(_left_cell_fiber_densities_filtered, _n=_derivative)
                _right_cell_fiber_densities_derivative = \
                    compute_lib.derivative(_right_cell_fiber_densities_filtered, _n=_derivative)

                if ADF_TEST:
                    _, _left_cell_adf_p_value, _, _, _, _ = adfuller(_left_cell_fiber_densities_derivative)
                    _, _right_cell_adf_p_value, _, _, _, _ = adfuller(_right_cell_fiber_densities_derivative)
                    if _left_cell_adf_p_value > 0.05 or _right_cell_adf_p_value > 0.05:
                        continue

                if KPSS_TEST:
                    _, _left_cell_kpss_p_value, _, _ = kpss(_left_cell_fiber_densities_derivative, nlags='legacy')
                    _, _right_cell_kpss_p_value, _, _ = kpss(_right_cell_fiber_densities_derivative, nlags='legacy')
                    if _left_cell_kpss_p_value < 0.05 or _right_cell_kpss_p_value < 0.05:
                        continue

                # stationary
                break

        # causality
        try:
            _x = pd.DataFrame(
                data=[[_left_value, _right_value] for _left_value, _right_value in
                      zip(_left_cell_fiber_densities_derivative, _right_cell_fiber_densities_derivative)],
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

                    # time lag = 0
                    _correlation = compute_lib.correlation(
                        _left_cell_fiber_densities_derivative, _right_cell_fiber_densities_derivative)

                    # if _correlation < 0.5:
                    #     continue

                    # granger causality
                    for _caused, _causing in zip(['left', 'right'], ['right', 'left']):
                        _granger = _var_model_results.test_causality(caused=_caused, causing=_causing)
                        _granger_causality_p_values.append(_granger.pvalue)

                        # time lag = 0
                        _correlations.append(_correlation)

                        # time lag = min estimator
                        if _causing == 'left':
                            _left_fiber_densities_time_lag = \
                                _left_cell_fiber_densities_derivative[:-_min_estimator_lag]
                            _right_fiber_densities_time_lag = \
                                _right_cell_fiber_densities_derivative[_min_estimator_lag:]
                        else:
                            _left_fiber_densities_time_lag = \
                                _left_cell_fiber_densities_derivative[_min_estimator_lag:]
                            _right_fiber_densities_time_lag = \
                                _right_cell_fiber_densities_derivative[:-_min_estimator_lag]
                        _time_lag_correlation = compute_lib.correlation(
                            _left_fiber_densities_time_lag, _right_fiber_densities_time_lag
                        )
                        _time_lag_correlations.append(_time_lag_correlation)

                        # end fiber density
                        _time_frame = compute.density_time_frame(_experiment)
                        if len(_left_cell_fiber_densities_filtered) > _time_frame:
                            _end_fiber_density = \
                                (_left_cell_fiber_densities_filtered[_time_frame] +
                                 _right_cell_fiber_densities_filtered[_time_frame]) / 2
                        else:
                            _end_fiber_density = \
                                (_left_cell_fiber_densities_filtered[-1] +
                                 _right_cell_fiber_densities_filtered[-1]) / 2
                        _normalization = load.normalization_series_file_data(_experiment, _series_id)
                        _normalized_fiber_density = compute_lib.z_score(
                            _end_fiber_density,
                            _normalization['average'],
                            _normalization['std']
                        )
                        _end_fiber_densities.append(_normalized_fiber_density)

                        # marking
                        if _tuples_to_mark is not None and _tuple in _tuples_to_mark and _granger.pvalue < 0.05:
                            print(_tuple, 'causing:', _causing, 'marked granger p-value:', _granger.pvalue)

                        if _granger.pvalue < 0.05:
                            if _properties['band']:
                                _n_passed_granger_causality_with_band += 1

                            _normality = _var_model_results.test_normality()
                            _inst_granger = _var_model_results.test_inst_causality(causing=_causing)

                            print(_tuple, _causing.capitalize() + ' causes ' + _caused + '!',
                                  'time-points: ' + str(len(_left_cell_fiber_densities_derivative)),
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
                                    _left_fiber_densities_time_lag = _left_cell_fiber_densities_derivative[:-_lag]
                                    _right_fiber_densities_time_lag = _right_cell_fiber_densities_derivative[_lag:]
                                else:
                                    _left_fiber_densities_time_lag = _left_cell_fiber_densities_derivative[_lag:]
                                    _right_fiber_densities_time_lag = _right_cell_fiber_densities_derivative[:-_lag]

                                _correlation = compute_lib.correlation(
                                    _left_fiber_densities_time_lag, _right_fiber_densities_time_lag
                                )
                                print('Time lag = ' + str(_lag) + ' correlation:', _correlation)

                            # plots
                            if _tuples_to_plot is not None and _tuple in _tuples_to_plot:
                                _y_arrays = [_left_cell_fiber_densities_derivative,
                                             _right_cell_fiber_densities_derivative]
                                _names_array = ['Left cell', 'Right cell']
                                _colors_array = config.colors(2)
                                _temporal_resolution = compute.temporal_resolution_in_minutes(_experiment)
                                _fig = go.Figure(
                                    data=[
                                        go.Scatter(
                                            x=np.arange(
                                                start=_start_time_frame,
                                                stop=_start_time_frame + len(_left_cell_fiber_densities_derivative),
                                                step=1) * _temporal_resolution,
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
                                                start=_start_time_frame,
                                                stop=_start_time_frame + len(_y),
                                                step=1) * _temporal_resolution,
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
    print('GC vs. correlation pearson correlation:',
          compute_lib.correlation(_granger_causality_p_values, _correlations, _with_p_value=True))
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
    print('GC vs. time lag correlation pearson correlation:',
          compute_lib.correlation(_granger_causality_p_values, _time_lag_correlations, _with_p_value=True))
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
    print('GC vs. end fiber density pearson correlation:',
          compute_lib.correlation(_granger_causality_p_values, _end_fiber_densities, _with_p_value=True))
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

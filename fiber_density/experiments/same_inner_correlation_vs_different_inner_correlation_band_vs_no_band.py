import os

import plotly.graph_objs as go

from fiber_density.experiments import same_inner_correlation_vs_different_inner_correlation
from libs import compute_lib
from libs.experiments import paths
from plotting import save


def main(_high_time_resolution=False):
    _y_arrays = [[], []]
    for _band_index, _band in enumerate([True, False]):
        print('Band:', _band)
        _same_correlations_array, _different_correlations_array = \
            same_inner_correlation_vs_different_inner_correlation.compute_fiber_densities(_band=_band,
                                                                                          _high_time_resolution=_high_time_resolution)
        for _same, _different in zip(_same_correlations_array, _different_correlations_array):
            _point_distance = compute_lib.distance_from_a_point_to_a_line(
                _line=[-1, -1, 1, 1],
                _point=[_same, _different]
            )
            if _same > _different:
                _y_arrays[_band_index].append(_point_distance)
            else:
                _y_arrays[_band_index].append(-_point_distance)

    # plot
    _colors_array = ['#ea8500', '#844b00']
    _names_array = ['Band', 'No Band']
    _fig = go.Figure(
        data=[
            go.Box(
                y=_y_array,
                name=_name,
                boxpoints=False,
                line={
                    'width': 1,
                    'color': _color
                },
                showlegend=False
            ) for _y_array, _name, _color in zip(_y_arrays, _names_array, _colors_array)
        ],
        layout={
            'xaxis': {
                'zeroline': False
            },
            'yaxis': {
                'title': 'Same minus different correlation',
                'zeroline': False
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_high_time_res_' + str(_high_time_resolution)
    )


if __name__ == '__main__':
    main()

import os

import plotly.graph_objs as go

from fibers_density.experiments import same_vs_different
from libs import compute_lib
from libs.experiments import paths
from plotting import save


def main(_high_time_resolution=False):
    _y_arrays = [[], []]
    for _band_index, _band in enumerate([True, False]):
        print('Band:', _band)
        _same_correlations_array, _different_correlations_array = \
            same_vs_different.compute_fibers_densities(_band=_band, _high_time_resolution=_high_time_resolution)
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
    _colors_array = ['#844b00', '#edbc80']
    _fig = go.Figure(
        data=[
            go.Box(
                y=_y_array,
                name=_band,
                boxpoints=False,
                line={
                    'width': 1
                },
                showlegend=False
            ) for _y_array, _band, _color in zip(_y_arrays, [True, False], _colors_array)
        ],
        layout={
            'xaxis': {
                'title': 'Pairs with band',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Distance from y = x',
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

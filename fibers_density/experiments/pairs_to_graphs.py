import os
import numpy as np
import plotly
import plotly.graph_objs as go

from libs.experiments import paths
from plotting.save import div_to_website


def plot_best_combined(_experiment, _normalized_pairs):
    _best = [
        ('series_1', 'pair_2', 'green'),
        ('series_17', 'pair_2', 'red'),
        ('series_19', 'pair_3', 'blue')
    ]

    _plot_file_name = 'best_pairs'
    _experiment_path = paths.plots(_experiment)
    _plot_path = os.path.join(_experiment_path, _plot_file_name)

    _data = []

    for _series_pair in _best:
        _series, _pair_name, _color = _series_pair
        _pair = _normalized_pairs[_series][_pair_name]
        _tp_first, _tp_last = (_pair[0], _pair[1]) if _pair[0]['tp'] == '1' else (_pair[1], _pair[0])
        for _tp in [_tp_first, _tp_last]:
            _data.append(
                (go.Scatter(
                    x=np.arange(start=0.0, stop=100.0, step=0.125)[:len(_tp['normalized_data'])],
                    y=_tp['normalized_data'],
                    mode='lines',
                    line=dict(
                        color=_color,
                        dash='dash' if _tp['tp'] == '1' else 'solid'
                    ),
                    legendgroup=_series + '_' + _pair_name,
                    name='Time-Point ' + _tp['tp']
                ))
            )

    _fig = go.Figure(data=_data)

    _fig.update_layout(
        title='best pairs',
        xaxis_title='Distance from Left Cell (cell size)',
        yaxis_title='Normalized Fibers Density Change (%)'
    )

    _fig.update_layout(
        yaxis=dict(
            zerolinecolor='black',
            zerolinewidth=2,
            tickformat=',.0%'
        )
    )

    _fig.update_yaxes(range=[-0.15, 1.0])

    _div = plotly.offline.plot(figure_or_data=_fig, output_type='div', include_plotlyjs=False)
    div_to_website(_div, _plot_path)


if __name__ == '__main__':
    experiment = 'SN16_CZI'
    # plot_best_combined(experiment, normalized_pairs)

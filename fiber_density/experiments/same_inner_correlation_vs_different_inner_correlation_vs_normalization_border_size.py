import os
from itertools import product

import numpy as np
import plotly.graph_objs as go
import seaborn as sns

from fiber_density.experiments import same_inner_correlation_vs_different_inner_correlation
from libs.experiments import config, paths
from plotting import save

BY = np.arange(start=0, stop=1 + 1 / 16, step=1 / 16)


def main(_real_cells=True, _static=False, _band=True, _high_temporal_resolution=False):
    _z_array = np.zeros(shape=(len(BY), len(BY)))
    for (_padding_index, _padding_by), (_space_index, _space_by) in product(enumerate(BY), enumerate(BY)):
        print('Padding by: ', _padding_by, ', space by: ', _space_by)
        _same_correlations_array, _different_correlations_array = \
            same_inner_correlation_vs_different_inner_correlation.compute_fiber_densities(
                _real_cells=_real_cells, _static=_static, _band=_band, _high_temporal_resolution=_high_temporal_resolution,
                _padding_y_by=_padding_by, _padding_z_by=_padding_by, _space_y_by=_space_by, _space_z_by=_space_by
            )
        _same_minus_different = np.array(_same_correlations_array) - np.array(_different_correlations_array)
        _same_greater_than_different = (_same_minus_different > 0).sum() / len(_same_minus_different)
        _z_array[_space_index, _padding_index] = _same_greater_than_different

    # plot
    _colors_array = ['black', 'white', config.colors(1)]
    _fig = go.Figure(
        data=go.Heatmap(
            x=BY,
            y=BY,
            z=_z_array,
            colorscale=sns.color_palette(_colors_array).as_hex(),
            colorbar={
                'tickmode': 'array',
                'tickvals': [0.5, 0.75, 1],
                'ticktext': ['0.5', 'Same > different', '1'],
                'tickangle': -90
            },
            showscale=True,
            zmin=0.5,
            zmax=1
        ),
        layout={
            'xaxis': {
                'title': 'Border size (cell diameter)',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 0.5, 1, 1.5, 2]
            },
            'yaxis': {
                'title': 'Space from window (cell diameter)',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 0.5, 1, 1.5, 2]
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_real_' + str(_real_cells) + '_static_' + str(_static) + '_band_' + str(_band) +
                  '_high_time_' + str(_high_temporal_resolution)
    )


if __name__ == '__main__':
    main()

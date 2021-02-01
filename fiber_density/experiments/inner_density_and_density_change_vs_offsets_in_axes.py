import os

import numpy as np
import plotly.graph_objs as go
import seaborn as sns

from fiber_density.experiments import inner_density_vs_offsets_in_axes, inner_density_change_vs_offsets_in_axes
from fiber_density.experiments.inner_density_vs_offsets_in_axes import OFFSET_Y_START, OFFSET_Y_END, OFFSET_Y_STEP, \
    OFFSET_Z_START, OFFSET_Z_END, OFFSET_Z_STEP
from libs.experiments import paths, config
from plotting import save


def main(_band=True, _high_temporal_resolution=False):
    print('Computing fiber densities vs. offsets in axes:')
    _fiber_densities_z_array = inner_density_vs_offsets_in_axes.compute_z_array(
        _band=_band, _high_temporal_resolution=_high_temporal_resolution)

    print('Computing fiber density changes vs. offsets in axes:')
    _fiber_density_changes_z_array = inner_density_change_vs_offsets_in_axes.compute_z_array(
        _band=_band, _high_temporal_resolution=_high_temporal_resolution)

    _z_array = _fiber_density_changes_z_array / _fiber_densities_z_array

    # plot
    _offsets_y = np.arange(start=OFFSET_Y_START, stop=OFFSET_Y_END + OFFSET_Y_STEP, step=OFFSET_Y_STEP)
    _offsets_z = np.arange(start=OFFSET_Z_START, stop=OFFSET_Z_END + OFFSET_Z_STEP, step=OFFSET_Z_STEP)
    _colors_array = ['white', config.colors(1)]
    _fig = go.Figure(
        data=go.Heatmap(
            x=_offsets_z,
            y=_offsets_y,
            z=_z_array,
            colorscale=sns.color_palette(_colors_array).as_hex(),
            colorbar={
                'tickmode': 'array',
                'tickvals': [0, 0.25, 0.5],
                'ticktext': ['0.0', 'Z-score change / z-score', '0.5'],
                'tickangle': -90
            },
            showscale=True,
            zmin=-0,
            zmax=0.5
        ),
        layout={
            'xaxis': {
                'title': 'Offset in XY axis (cell diameter)',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [-4, -2, 0, 2, 4]
            },
            'yaxis': {
                'title': 'Offset in Z axis (cell diameter)',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [-1, 0, 1, 2]
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_high_time_' + str(_high_temporal_resolution) + '_band_' + str(_band)
    )


if __name__ == '__main__':
    main()

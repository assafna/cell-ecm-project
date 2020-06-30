import os

import plotly.graph_objs as go
from scipy.stats import pearsonr

from fibers_density.experiments import fibers_vs_offsets_in_axes, same_vs_different_offsets_in_axes
from libs.experiments import paths
from plotting import save

OFFSET_Y_START = -1
OFFSET_Y_END = 1
OFFSET_Z_START = -1
OFFSET_Z_END = 1


def main():
    print('Computing fibers vs. offsets in axes:')
    _fibers_z_array = fibers_vs_offsets_in_axes.compute_z_array(
        _offset_y_start=OFFSET_Y_START, _offset_y_end=OFFSET_Y_END,
        _offset_z_start=OFFSET_Z_START, _offset_z_end=OFFSET_Z_END).flatten()

    print('Computing "same vs. different" vs. offset in axes:')
    _same_vs_different_z_array = \
        same_vs_different_offsets_in_axes.compute_z_array(
            _offset_y_start=OFFSET_Y_START, _offset_y_end=OFFSET_Y_END,
            _offset_z_start=OFFSET_Z_START, _offset_z_end=OFFSET_Z_END).flatten()

    print('Correlation:', pearsonr(_fibers_z_array, _same_vs_different_z_array))

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_fibers_z_array,
            y=_same_vs_different_z_array,
            mode='markers',
            marker={
                'size': 5,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Fiber density (z-score)',
                'zeroline': False
            },
            'yaxis': {
                'title': '"same" > "different" (%)',
                'zeroline': False
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()

import os

import numpy as np

import plotting.update
from plotting import save, scatter
from libs.experiments.compute import fibers_density_cut_edges
from libs.compute_lib import z_score_fibers_densities_array
from libs.experiments import config, load, paths


def get_first_and_last_tps(_normalized_fibers_density):
    _tp_first, _tp_last = min(_normalized_fibers_density.keys()), max(_normalized_fibers_density.keys())
    return (
        {'tp': _tp_first, 'normalized_data': _normalized_fibers_density[_tp_first]},
        {'tp': _tp_last, 'normalized_data': _normalized_fibers_density[_tp_last]},
    )


def plot(_experiment, _series, _group, _z_group, _fibers_density, _normalization):
    _normalized_fibers_density = z_score_fibers_densities_array(_fibers_density, _normalization)
    _normalized_fibers_density_cut = fibers_density_cut_edges(_normalized_fibers_density)
    _tp_first, _tp_last = get_first_and_last_tps(_normalized_fibers_density_cut)
    _fig = scatter.create_plot(
        _x_array=[np.arange(start=0.0, stop=100.0, step=0.125)[:len(_tp_first['normalized_data'])]] * 2,
        _y_array=[
            _tp_last['normalized_data'],
            _tp_first['normalized_data']
        ],
        _names_array=[
            'Time-Point ' + str(_tp_last['tp']),
            'Time-Point ' + str(_tp_first['tp'])
        ],
        _modes_array=[
            'lines+markers',
            'lines+markers'
        ],
        _show_legend_array=[True] * 2,
        _x_axis_title='Distance from left cell (cell size)',
        _y_axis_title='Fiber density (z-score)',
        _title=_experiment + ' - ' + _series + ' - ' + _group + ' - ' + _z_group + ' - TP First vs. Last'
    )

    _fig = plotting.update.y_axis(
        _fig=_fig,
        _range=[0, 25],
        _color='black',
        _width=2
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(os.path.join(os.path.join(paths.plots(_experiment), _series), _group), _z_group),
        _filename='TP First vs Last'
    )


def main():
    for _experiment in paths.folders(paths.FIBERS_DENSITY):
        print('Experiment', _experiment)
        for _series in paths.folders(paths.fibers_density(_experiment)):
            print(_series)
            _normalization = load.normalization_series_file_data(_experiment, _series)
            for _group in paths.folders(paths.fibers_density(_experiment, _series)):
                print('Group', _group)
                for _z_group in paths.folders(paths.fibers_density(_experiment, _series, _group)):
                    print('Z Group', _z_group)
                    _fibers_density = load.fibers_density_z_group_file_data(_experiment, _series, _group, _z_group)
                    plot(_experiment, _series, _group, _z_group, _fibers_density, _normalization)


if __name__ == '__main__':
    main()

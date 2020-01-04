import os

import numpy as np

from plotting import save, scatter
from libs.experiments.compute import normalized_fibers_density, fibers_density_cut_edges
from libs.experiments import config, load, paths


def get_first_and_last_tps(_normalized_fibers_density):
    _tp_first, _tp_last = min(_normalized_fibers_density.keys()), max(_normalized_fibers_density.keys())
    return (
        {'tp': _tp_first, 'normalized_data': _normalized_fibers_density[_tp_first]},
        {'tp': _tp_last, 'normalized_data': _normalized_fibers_density[_tp_last]},
    )


def plot(_experiment, _series, _group, _z_group, _fibers_density, _normalization):
    _normalized_fibers_density = normalized_fibers_density(_fibers_density, _normalization)
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
        _x_axis_title='Distance from Left Cell (cell size)',
        _y_axis_title='Normalized Fibers Density Change (%)',
        _title=_experiment + ' - ' + _series + ' - ' + _group + ' - ' + _z_group + ' - TP First vs. Last'
    )

    _fig = scatter.update_y_axis(
        _fig=_fig,
        _color='black',
        _width=2,
        _format=',.0%',
        _range=[-0.15, 1.5]
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(os.path.join(os.path.join(paths.plots(_experiment), _series), _group), _z_group),
        _filename='TP First vs Last'
    )


def plot_all(_experiment, _experiment_fibers_density, _experiment_normalization):
    for _series in _experiment_fibers_density:
        _series_fibers_density = _experiment_fibers_density[_series]
        _series_normalization = experiment_normalization[_series]
        for _group in _series_fibers_density:
            _group_fibers_density = _series_fibers_density[_group]
            _group_normalization = _series_normalization[_group]
            for _z_group in _group_fibers_density:
                print('Plotting', _series, _group, _z_group)
                _fibers_density = _group_fibers_density[_z_group]
                plot(_experiment, _series, _group, _z_group, _fibers_density, _group_normalization)


if __name__ == '__main__':
    for experiment in config.SINGLE_CELL:
        experiment_fibers_density = load.fibers_density_experiment_file_data(experiment)
        experiment_normalization = load.experiment_normalization(experiment)
        plot_all(experiment, experiment_fibers_density, experiment_normalization)

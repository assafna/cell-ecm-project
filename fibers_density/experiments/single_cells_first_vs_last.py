import numpy as np

import libs.experiments.compute
from plotting import save, scatter
from libs.experiments.compute import normalized_fibers_density, fibers_density_cut_left_edge
from libs.experiments import config, load, paths


def per_cell(_experiment_fibers_density, _experiment_normalization):
    _experiment_fibers_density_per_cell = {}
    for _series in _experiment_fibers_density:
        _series_fibers_density = _experiment_fibers_density[_series]
        _series_normalization = _experiment_normalization[_series]
        _series_fibers_density_per_cell = {}
        for _group in _series_fibers_density:
            _group_fibers_density = _series_fibers_density[_group]
            _group_normalization = _series_normalization[_group]
            _cell_name = str(_group.split('_degree')[0])
            if _cell_name not in _series_fibers_density_per_cell:
                _series_fibers_density_per_cell[_cell_name] = []
            for _z_group in _group_fibers_density:
                _z_group_fibers_density = _group_fibers_density[_z_group]
                _fibers_density_normalized = normalized_fibers_density(_z_group_fibers_density, _group_normalization)
                _cut_left = fibers_density_cut_left_edge(_fibers_density_normalized)
                _series_fibers_density_per_cell[_cell_name].append(_cut_left)
        _experiment_fibers_density_per_cell[_series] = _series_fibers_density_per_cell

    return _experiment_fibers_density_per_cell


def all_directions_average(_per_cell):
    _all_averages = {}
    for _series in _per_cell:
        _series_fibers_density = _per_cell[_series]
        _series_averages = {}
        for _cell in _series_fibers_density:
            _cell_directions = _series_fibers_density[_cell]
            _tps_values = {}
            for _direction in _cell_directions:
                for _tp in _direction:
                    _tp_direction_fibers_density = _direction[_tp]
                    if _tp not in _tps_values:
                        _tps_values[_tp] = [_tp_direction_fibers_density]
                    else:
                        _tps_values[_tp].append(_tp_direction_fibers_density)
            # average
            _tps_averages = {}
            for _tp in _tps_values:
                _tp_fibers_density = _tps_values[_tp]
                _max_distances = max([len(_direction) for _direction in _tp_fibers_density])
                _tps_averages[_tp] = list(np.zeros(_max_distances))
                for _distance_index in range(_max_distances):
                    _num = 0
                    _sum = 0
                    for _direction in _tp_fibers_density:
                        if len(_direction) > _distance_index:
                            _sum += _direction[_distance_index]
                            _num += 1
                    _tps_averages[_tp][_distance_index] = _sum / _num
            _series_averages[_cell] = _tps_averages
        _all_averages[_series] = _series_averages

    return _all_averages


def plot_average(_experiment, _fibers_density_per_cell_averages):
    _tp_1_data = []
    _tp_latest_data = []
    for _series in _fibers_density_per_cell_averages:
        _series_fibers_density = _fibers_density_per_cell_averages[_series]
        for _cell in _series_fibers_density:
            _cell_fibers_density = _series_fibers_density[_cell]
            _tp_1_fibers_density = _cell_fibers_density[1]
            _tp_latest_fibers_density = _cell_fibers_density[max(_cell_fibers_density.keys())]
            for _distance_index, _distance in enumerate(_tp_1_fibers_density):
                if len(_tp_1_data) < _distance_index + 1:
                    _tp_1_data.append([_distance])
                else:
                    _tp_1_data[_distance_index].append(_distance)
            for _distance_index, _distance in enumerate(_tp_latest_fibers_density):
                if len(_tp_latest_data) < _distance_index + 1:
                    _tp_latest_data.append([_distance])
                else:
                    _tp_latest_data[_distance_index].append(_distance)

    _fig = scatter.create_error_bars_plot(
        _x_array=[
            np.arange(start=0.0, stop=100.0, step=0.125)[:len(_tp_latest_data)],
            np.arange(start=0.025, stop=100.0, step=0.125)[:len(_tp_1_data)]
        ],
        _y_array=[
            _tp_latest_data,
            _tp_1_data
        ],
        _names_array=[
            'Time-Point Last',
            'Time-Point 1'
        ],
        _mode_array=[
            'lines+markers',
            'lines+markers'
        ],
        _dash_array=[
            'solid',
            'dash'
        ],
        _x_axis_title='Distance from Left Cell (cell size)',
        _y_axis_title='Normalized Fibers Density Change (%)',
        _title=_experiment + ' - TP First vs. Last'
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
        _path=paths.plots(_experiment),
        _filename='TP First vs Last'
    )


def plot_deltas_average(_experiment, _fibers_density_per_cell_averages):
    _delta_data = []
    for _series in _fibers_density_per_cell_averages:
        _series_fibers_density = _fibers_density_per_cell_averages[_series]
        for _cell in _series_fibers_density:
            _cell_fibers_density = _series_fibers_density[_cell]
            _tp_1_fibers_density = _cell_fibers_density[1]
            _tp_latest_fibers_density = _cell_fibers_density[max(_cell_fibers_density.keys())]
            for (_distance_index, _tp_1_distance), _tp_latest_distance in zip(enumerate(_tp_1_fibers_density), _tp_latest_fibers_density):
                _delta = _tp_latest_distance - _tp_1_distance
                if len(_delta_data) < _distance_index + 1:
                    _delta_data.append([_delta])
                else:
                    _delta_data[_distance_index].append(_delta)

    _fig = scatter.create_error_bars_plot(
        _x_array=[
            np.arange(start=0.0, stop=100.0, step=0.125)[:len(_delta_data)]
        ],
        _y_array=[
            _delta_data
        ],
        _names_array=[
            'Deltas'
        ],
        _mode_array=[
            'lines+markers'
        ],
        _dash_array=[
            'solid'
        ],
        _x_axis_title='Distance from Left Cell (cell size)',
        _y_axis_title='Normalized Fibers Density Change (%)',
        _title=_experiment + ' - TP First vs. Last - Deltas'
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
        _path=paths.plots(_experiment),
        _filename='TP First vs Last - Deltas'
    )


if __name__ == '__main__':
    for experiment in config.SINGLE_CELL:
        experiment_fibers_density = load.experiment(experiment)
        experiment_normalization = libs.experiments.compute.experiment(experiment)
        fibers_density_per_cell = per_cell(experiment_fibers_density, experiment_normalization)
        fibers_density_per_cell_averages = all_directions_average(fibers_density_per_cell)
        plot_average(experiment, fibers_density_per_cell_averages)
        # plot_deltas_average(experiment, fibers_density_per_cell_averages)

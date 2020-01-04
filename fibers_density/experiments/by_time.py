import os

import numpy as np
import seaborn

from libs.experiments import load, paths
from libs.experiments.compute import z_score_fibers_density_array, fibers_density_cut_edges
from plotting import heatmap, save


def plot(_experiment, _series, _group, _z_group, _fibers_density, _normalization):
    _normalized_fibers_density = z_score_fibers_density_array(_fibers_density, _normalization)
    _normalized_fibers_density_cut = fibers_density_cut_edges(_normalized_fibers_density)
    _sorted_by_time = sorted(_normalized_fibers_density_cut)
    _max_distances = max([len(_tp) for _tp in _normalized_fibers_density_cut.values()])
    _z_array = []
    for _tp in _sorted_by_time:
        _data = list(_normalized_fibers_density_cut[_tp])
        for _i in range(len(_data), _max_distances):
            _data.append(None)
        _z_array.append(_data)

    _z_array = list(zip(*_z_array))

    _fig = heatmap.create_plot(
        _x_labels=np.arange(start=0, stop=10000, step=15),
        _y_labels=np.arange(start=0.0, stop=100.0, step=0.125),
        _z_array=_z_array,
        _x_axis_title='Time Since Imaging Start (minutes)',
        _y_axis_title='Distance from Left Cell Center (cell size)',
        _color_scale=seaborn.light_palette("navy", reverse=True).as_hex(),
        _show_scale=True,
        _title=_experiment + ' - ' + _series + ' - ' + _group + ' - ' + _z_group + ' - All TPs'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(os.path.join(os.path.join(paths.plots(_experiment), _series), _group), _z_group),
        _filename='All TPs'
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

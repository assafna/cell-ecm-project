import libs.experiments.compute
from libs.experiments.compute import normalized_fibers_density, fibers_density_cut_edges
from libs.experiments import config, load


def plot(_experiment, _series, _group, _z_group, _fibers_density, _normalization):
    _normalized_fibers_density = normalized_fibers_density(_fibers_density, _normalization)
    _normalized_fibers_density_cut = fibers_density_cut_edges(_normalized_fibers_density)
    _sorted_by_time = sorted(_normalized_fibers_density_cut)
    _max_distances = max([len(_tp) for _tp in _normalized_fibers_density_cut.values()])
    _z_array = []
    for _tp in _sorted_by_time:
        _data = list(_normalized_fibers_density_cut[_tp])
        if _series == 'Series 1' and _group == 'cells_3_4':
            print(_tp, _data[0], _data[-1])
        for _i in range(len(_data), _max_distances):
            _data.append(None)
        _z_array.append(_data)

    _z_array = list(zip(*_z_array))

    # _fig = heatmap.create_plot(
    #     _x_labels=np.arange(start=0, stop=10000, step=15),
    #     _y_labels=np.arange(start=0.0, stop=100.0, step=0.125),
    #     _z_array=_z_array,
    #     _x_axis_title='Time Since Imaging Start (minutes)',
    #     _y_axis_title='Distance from Left Cell Center (cell size)',
    #     _color_scale=seaborn.cubehelix_palette().as_hex(),
    #     _show_scale=True,
    #     _title=_experiment + ' - ' + _series + ' - ' + _group + ' - ' + _z_group + ' - All TPs'
    # )

    # save.save_plot(
    #     _fig=_fig,
    #     _path=os.path.join(os.path.join(os.path.join(paths.plots(_experiment), _series), _group), _z_group),
    #     _filename='All TPs'
    # )


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
    for experiment in config.PAIRS:
        experiment_fibers_density = load.experiment(experiment)
        experiment_normalization = libs.experiments.compute.experiment(experiment)
        plot_all(experiment, experiment_fibers_density, experiment_normalization)

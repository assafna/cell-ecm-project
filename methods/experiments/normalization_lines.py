import os

from libs.experiments import config, paths, load


def main():
    for _experiment in config.PAIRS:
        print('Experiment', _experiment)
        os.mkdir(paths.normalization_lines(_experiment)) if not os.path.isdir(
            paths.normalization_lines(_experiment)) else None
        _experiment_objects_path = os.path.join(paths.objects(_experiment, _is_z=False))
        _serieses = load.serieses(_experiment_objects_path)
        # for _series in _serieses:
        #     _tp_first_path =
        # coordinates = cell_coordinates_by_time(objects, _is_z=_is_z) if objects is not None else None
        # save_to_file(coordinates, experiment, _is_z=_is_z) if coordinates is not None else None



if __name__ == '__main__':
    main()

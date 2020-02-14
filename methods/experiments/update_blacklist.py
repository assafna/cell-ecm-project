import os

from libs.experiments import load, save, paths


def add_to_blacklist(_experiment, _series_id, _group, _time_point_start, _time_point_end, _reason):
    if os.path.isfile(paths.blacklist(_experiment, _series_id, _group)):
        _blacklist = load.blacklist(_experiment, _series_id, _group)
    else:
        _blacklist = {}

    for _time_point in range(_time_point_start, _time_point_end + 1):
        _blacklist[_time_point] = _reason
        print('Added:', _experiment, _series_id, _group, _reason, sep='\t')

    save.blacklist(_experiment, _series_id, _group, _blacklist)
    print('Saved:', _experiment, _series_id, _group)


def remove_from_blacklist(_experiment, _series_id, _group, _time_point_start, _time_point_end):
    if os.path.isfile(paths.blacklist(_experiment, _series_id, _group)):
        _blacklist = load.blacklist(_experiment, _series_id, _group)

        for _time_point in range(_time_point_start, _time_point_end + 1):
            _blacklist = {}
            _reason = _blacklist.pop(_time_point, None)
            print('Removed:', _experiment, _series_id, _group, _reason, sep='\t')

        save.blacklist(_experiment, _series_id, _group, _blacklist)
        print('Saved:', _experiment, _series_id, _group)
    else:
        print('No such blacklist')


if __name__ == '__main__':
    add_to_blacklist(_experiment='SN16',
                     _series_id=1,
                     _group='cells_0_1',
                     _time_point_start=0,
                     _time_point_end=10,
                     _reason='FAKE REASON'
                     )

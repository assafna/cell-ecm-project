import os

from libs.experiments import load, save, paths


def add_to_blacklist(_experiment, _series_id, _cell_id, _time_point_start, _time_point_end, _reason):
    if os.path.isfile(paths.blacklist(_experiment, _series_id)):
        _blacklist = load.blacklist(_experiment, _series_id)
    else:
        _blacklist = {}

    for _time_point in range(_time_point_start, _time_point_end + 1):
        if _cell_id in _blacklist:
            if _time_point in _blacklist[_cell_id]:
                _blacklist[_cell_id][_time_point].append(_reason)
            else:
                _blacklist[_cell_id].append({_time_point: [_reason]})
        else:
            _blacklist[_cell_id] = [{_time_point: [_reason]}]
    print('Added:', _experiment, _series_id, _cell_id, (_time_point_start, _time_point_end, _reason), sep='\t')

    save.blacklist(_experiment, _series_id, _blacklist)
    print('Saved:', _experiment, _series_id)


def remove_from_blacklist(_experiment, _series_id, _cell_id, _time_point_start, _time_point_end):
    if os.path.isfile(paths.blacklist(_experiment, _series_id)):
        _blacklist = load.blacklist(_experiment, _series_id)

        for _time_point in range(_time_point_start, _time_point_end + 1):
            if _cell_id in _blacklist:
                if _time_point in _blacklist[_cell_id]:
                    _reason = _blacklist[_cell_id].pop(_time_point, None)
                    print('Removed:', _experiment, _series_id, _cell_id, _time_point, _reason, sep='\t')

        save.blacklist(_experiment, _series_id, _blacklist)
        print('Saved:', _experiment, _series_id)
    else:
        print('No such blacklist')


if __name__ == '__main__':
    add_to_blacklist(_experiment='SN41',
                     _series_id=3,
                     _cell_id=None,
                     _time_point_start=34,
                     _time_point_end=34,
                     _reason='Black image starting from Z slice #22'
                     )
    add_to_blacklist(_experiment='SN41',
                     _series_id=3,
                     _cell_id=None,
                     _time_point_start=2,
                     _time_point_end=4,
                     _reason='Wave of light'
                     )
    add_to_blacklist(_experiment='SN16',
                     _series_id=1,
                     _cell_id=4,
                     _time_point_start=21,
                     _time_point_end=31,
                     _reason='Cell expands'
                     )
    add_to_blacklist(_experiment='SN16',
                     _series_id=2,
                     _cell_id=1,
                     _time_point_start=21,
                     _time_point_end=31,
                     _reason='Cell expands'
                     )
    add_to_blacklist(_experiment='SN16',
                     _series_id=3,
                     _cell_id=2,
                     _time_point_start=18,
                     _time_point_end=31,
                     _reason='Cell expands'
                     )
    add_to_blacklist(_experiment='SN16',
                     _series_id=7,
                     _cell_id=1,
                     _time_point_start=10,
                     _time_point_end=31,
                     _reason='Cell expands'
                     )
    add_to_blacklist(_experiment='SN16',
                     _series_id=9,
                     _cell_id=1,
                     _time_point_start=23,
                     _time_point_end=31,
                     _reason='Cell expands'
                     )
    add_to_blacklist(_experiment='SN16',
                     _series_id=17,
                     _cell_id=1,
                     _time_point_start=18,
                     _time_point_end=31,
                     _reason='Cell expands'
                     )
    add_to_blacklist(_experiment='SN16',
                     _series_id=17,
                     _cell_id=3,
                     _time_point_start=20,
                     _time_point_end=31,
                     _reason='Cell expands'
                     )
    add_to_blacklist(_experiment='SN16',
                     _series_id=17,
                     _cell_id=3,
                     _time_point_start=20,
                     _time_point_end=31,
                     _reason='Cell attaches to another cell'
                     )
    add_to_blacklist(_experiment='SN16',
                     _series_id=21,
                     _cell_id=1,
                     _time_point_start=21,
                     _time_point_end=31,
                     _reason='Cell expands'
                     )
    add_to_blacklist(_experiment='SN16',
                     _series_id=21,
                     _cell_id=2,
                     _time_point_start=23,
                     _time_point_end=31,
                     _reason='Cell expands'
                     )
    add_to_blacklist(_experiment='SN16',
                     _series_id=21,
                     _cell_id=3,
                     _time_point_start=19,
                     _time_point_end=31,
                     _reason='Cell expands'
                     )

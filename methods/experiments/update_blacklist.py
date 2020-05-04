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
                _blacklist[_cell_id][_time_point] = [_reason]
        else:
            _blacklist[_cell_id] = {_time_point: [_reason]}
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
    for _cell_id in [0, 1, 2, 3]:
        add_to_blacklist(
            _experiment='SN44',
            _series_id=2,
            _cell_id=_cell_id,
            _time_point_start=1,
            _time_point_end=40,
            _reason='Light wave'
        )

    for _cell_id in [0, 1, 2, 3, 4]:
        add_to_blacklist(
            _experiment='SN44',
            _series_id=3,
            _cell_id=_cell_id,
            _time_point_start=1,
            _time_point_end=115,
            _reason='Light wave'
        )

    for _cell_id in [0, 1, 2]:
        add_to_blacklist(
            _experiment='SN44',
            _series_id=4,
            _cell_id=_cell_id,
            _time_point_start=9,
            _time_point_end=11,
            _reason='Light wave'
        )

    for _cell_id in [0, 1, 2, 3]:
        add_to_blacklist(
            _experiment='SN44',
            _series_id=5,
            _cell_id=_cell_id,
            _time_point_start=8,
            _time_point_end=10,
            _reason='Light wave'
        )

    for _cell_id in [1, 3, 4, 6]:
        add_to_blacklist(
            _experiment='SN44',
            _series_id=6,
            _cell_id=_cell_id,
            _time_point_start=24,
            _time_point_end=27,
            _reason='Light wave'
        )
        add_to_blacklist(
            _experiment='SN44',
            _series_id=6,
            _cell_id=_cell_id,
            _time_point_start=31,
            _time_point_end=33,
            _reason='Light wave'
        )

    for _cell_id in [0, 1, 2, 3]:
        add_to_blacklist(
            _experiment='SN44',
            _series_id=7,
            _cell_id=_cell_id,
            _time_point_start=27,
            _time_point_end=27,
            _reason='Light wave'
        )

    add_to_blacklist(
        _experiment='SN44',
        _series_id=7,
        _cell_id=3,
        _time_point_start=111,
        _time_point_end=247,
        _reason='Cell expands'
    )

    for _cell_id in [0, 1, 2, 3]:
        add_to_blacklist(
            _experiment='SN44',
            _series_id=8,
            _cell_id=_cell_id,
            _time_point_start=27,
            _time_point_end=27,
            _reason='Light wave'
        )

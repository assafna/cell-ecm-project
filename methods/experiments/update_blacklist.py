import os

from libs.experiments import load, save, paths


def add_to_blacklist(_experiment, _series_id, _cell_id, _time_frame_start, _time_frame_end, _reason):
    if os.path.isfile(paths.blacklist(_experiment, _series_id)):
        _blacklist = load.blacklist(_experiment, _series_id)
    else:
        _blacklist = {}

    for _time_frame in range(_time_frame_start, _time_frame_end + 1):
        if _cell_id in _blacklist:
            if _time_frame in _blacklist[_cell_id]:
                _blacklist[_cell_id][_time_frame].append(_reason)
            else:
                _blacklist[_cell_id][_time_frame] = [_reason]
        else:
            _blacklist[_cell_id] = {_time_frame: [_reason]}
    print('Added:', _experiment, _series_id, _cell_id, (_time_frame_start, _time_frame_end, _reason), sep='\t')

    save.blacklist(_experiment, _series_id, _blacklist)
    print('Saved:', _experiment, _series_id)


def remove_from_blacklist(_experiment, _series_id, _cell_id, _time_frame_start, _time_frame_end):
    if os.path.isfile(paths.blacklist(_experiment, _series_id)):
        _blacklist = load.blacklist(_experiment, _series_id)

        for _time_frame in range(_time_frame_start, _time_frame_end + 1):
            if _cell_id in _blacklist:
                if _time_frame in _blacklist[_cell_id]:
                    _reason = _blacklist[_cell_id].pop(_time_frame, None)
                    print('Removed:', _experiment, _series_id, _cell_id, _time_frame, _reason, sep='\t')

        save.blacklist(_experiment, _series_id, _blacklist)
        print('Saved:', _experiment, _series_id)
    else:
        print('No such blacklist')


if __name__ == '__main__':
    _tuples = [
        ('SN45', 6, 0, 4, 7, 'Light wave'),
        ('SN45', 6, 1, 4, 7, 'Light wave'),
        ('SN45', 6, 2, 4, 7, 'Light wave'),
        ('SN45', 6, 3, 4, 7, 'Light wave'),
        ('SN45', 7, 0, 10, 15, 'Light wave'),
        ('SN45', 7, 1, 10, 15, 'Light wave'),
        ('SN45', 8, 0, 10, 20, 'Light wave'),
        ('SN45', 8, 1, 10, 20, 'Light wave'),
        ('SN45', 8, 2, 10, 20, 'Light wave'),
        ('SN45', 8, 3, 10, 20, 'Light wave'),
        ('SN45', 9, 0, 4, 9, 'Light wave'),
        ('SN45', 9, 1, 4, 9, 'Light wave'),
        ('SN45', 9, 2, 4, 9, 'Light wave'),
        ('SN45', 1, 0, 9, 24, 'Light wave'),
        ('SN45', 1, 2, 9, 24, 'Light wave'),
        ('SN45', 4, 0, 0, 5, 'Light wave'),
        ('SN45', 4, 3, 0, 5, 'Light wave'),
        ('SN45', 5, 1, 0, 10, 'Light wave'),
        ('SN45', 5, 4, 0, 10, 'Light wave')
    ]
    for _tuple in _tuples:
        add_to_blacklist(
            _experiment=_tuple[0],
            _series_id=_tuple[1],
            _cell_id=_tuple[2],
            _time_frame_start=_tuple[3],
            _time_frame_end=_tuple[4],
            _reason=_tuple[5]
        )

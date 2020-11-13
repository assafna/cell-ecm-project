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
        ('SN26_BlebAdded', 1, 2, 11, 12, 'Light wave'),
        ('SN26_BlebAdded', 1, 3, 11, 12, 'Light wave'),
        ('SN26_BlebAdded', 2, 1, 3, 6, 'Light wave'),
        ('SN26_BlebAdded', 2, 2, 3, 6, 'Light wave'),
        ('SN26_BlebAdded', 2, 3, 3, 6, 'Light wave'),
        ('SN26_BlebAdded', 2, 4, 3, 6, 'Light wave'),
        ('SN26_BlebAdded', 4, 3, 17, 22, 'Bad segmentation'),
        ('SN26_BlebAdded', 7, 0, 2, 4, 'Light wave'),
        ('SN26_BlebAdded', 7, 1, 2, 4, 'Light wave'),
        ('SN26_BlebAdded', 7, 2, 2, 4, 'Light wave'),
        ('SN26_BlebAdded', 7, 3, 2, 4, 'Light wave'),
        ('SN26_BlebAdded', 8, 1, 2, 3, 'Light wave'),
        ('SN26_BlebAdded', 8, 2, 2, 3, 'Light wave'),
        ('SN26_BlebAdded', 9, 0, 2, 3, 'Light wave'),
        ('SN26_BlebAdded', 9, 1, 2, 3, 'Light wave'),
        ('SN26_BlebAdded', 9, 3, 2, 3, 'Light wave'),
        ('SN26_BlebAdded', 13, 0, 2, 2, 'Light wave'),
        ('SN26_BlebAdded', 13, 3, 2, 2, 'Light wave'),
        ('SN26_BlebAdded', 13, 4, 2, 2, 'Light wave'),
        ('SN26_BlebAdded', 13, 4, 1, 1, 'Bad segmentation'),
        ('SN26_BlebAdded', 15, 0, 0, 2, 'Light wave'),
        ('SN26_BlebAdded', 15, 1, 0, 2, 'Light wave'),
        ('SN26_BlebAdded', 15, 2, 0, 2, 'Light wave'),
        ('SN26_BlebAdded', 16, 2, 11, 22, 'Bad segmentation'),
        ('SN26_BlebAdded', 22, 1, 0, 4, 'Light wave'),
        ('SN26_BlebAdded', 22, 2, 0, 4, 'Light wave'),
        ('SN26_BlebAdded', 22, 3, 0, 4, 'Light wave'),
        ('SN26_BlebAdded', 23, 0, 0, 0, 'Bad segmentation'),
        ('SN26_BlebAdded', 23, 0, 1, 1, 'Light wave'),
        ('SN26_BlebAdded', 23, 1, 1, 1, 'Light wave'),
        ('SN26_BlebAdded', 23, 3, 1, 1, 'Light wave'),
        ('SN26_BlebAdded', 23, 4, 1, 1, 'Light wave'),
        ('SN26_BlebAdded', 23, 0, 10, 12, 'Light wave'),
        ('SN26_BlebAdded', 23, 1, 10, 12, 'Light wave'),
        ('SN26_BlebAdded', 23, 3, 10, 12, 'Light wave'),
        ('SN26_BlebAdded', 23, 4, 10, 12, 'Light wave'),
        ('SN26_BlebAdded', 23, 3, 13, 13, 'Bad segmentation')
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

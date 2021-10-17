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
        ('Metastasis_Airyscan_210401_sgAXL', 5, 3, 12, 19, 'Bad segmentation'),
        ('Metastasis_Airyscan_210401_sgAXL', 9, 1, 0, 0, 'Bad segmentation'),
        ('Metastasis_Airyscan_210401_sgAXL', 12, 2, 19, 19, 'Bad segmentation')
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

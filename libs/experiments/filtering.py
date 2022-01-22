from libs.experiments import load, compute
from libs.experiments.config import SINGLE_CELL, BLEB, BLEB_FROM_START, HIGH_TEMPORAL_RESOLUTION_IN_MINUTES, \
    LIVE_DEAD, BEAD, METASTASIS, DEAD_DEAD, BLEB_AMOUNT_UM


def is_single_cell(_experiment):
    return _experiment in SINGLE_CELL


def is_high_temporal_resolution(_experiment):
    return compute.temporal_resolution_in_minutes(_experiment) == HIGH_TEMPORAL_RESOLUTION_IN_MINUTES


def is_bleb(_experiment):
    return _experiment in BLEB


def is_dead_dead(_experiment):
    return _experiment in DEAD_DEAD


def is_live_dead(_experiment):
    return _experiment in LIVE_DEAD


def is_bead(_experiment):
    return _experiment in BEAD


def is_metastasis(_experiment):
    return _experiment in METASTASIS


def is_real(_tuple, _real=True):
    _experiment, _series_id, _group = _tuple
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    return _real == ('real_fake' not in _group_properties)


def is_band(_tuple, _band=True):
    _experiment, _series_id, _group = _tuple
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    return _band == _group_properties['band']


def by_categories(_experiments, _is_single_cell=None, _is_high_temporal_resolution=None, _is_bleb=None,
                  _is_dead_dead=None, _is_live_dead=None, _is_bead=None, _is_metastasis=None):
    return [_experiment for _experiment in _experiments if
            (_is_single_cell is None or _is_single_cell == is_single_cell(_experiment)) and
            (_is_high_temporal_resolution is None or _is_high_temporal_resolution == is_high_temporal_resolution(
                _experiment)) and
            (_is_bleb is None or (_is_bleb and is_bleb(_experiment)) or (not _is_bleb and not is_bleb(_experiment))) and
            (_is_dead_dead is None or _is_dead_dead == is_dead_dead(_experiment)) and
            (_is_live_dead is None or _is_live_dead == is_live_dead(_experiment)) and
            (_is_bead is None or _is_bead == is_bead(_experiment)) and
            (_is_metastasis is None or _is_metastasis == is_metastasis(_experiment))
            ]


def by_pair_distance(_experiments_tuples, _distance, _time_frame=0, _range=1):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        if len(_group_properties['time_points']) <= _time_frame:
            continue
        _left_cell_coordinates = _group_properties['time_points'][_time_frame]['left_cell']['coordinates']
        _right_cell_coordinates = _group_properties['time_points'][_time_frame]['right_cell']['coordinates']
        _pair_distance = compute.pair_distance_in_cell_size(
            _experiment=_experiment,
            _series_id=_series_id,
            _cell_1_coordinates=
            [(_left_cell_coordinates['x'], _left_cell_coordinates['y'], _left_cell_coordinates['z'])],
            _cell_2_coordinates=
            [(_right_cell_coordinates['x'], _right_cell_coordinates['y'], _right_cell_coordinates['z'])]
        )
        if _distance - _range <= _pair_distance <= _distance + _range:
            _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered


def by_pair_distance_range(_experiments_tuples, _distance_range, _time_frame=0):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        if len(_group_properties['time_points']) <= _time_frame:
            continue
        _left_cell_coordinates = _group_properties['time_points'][_time_frame]['left_cell']['coordinates']
        _right_cell_coordinates = _group_properties['time_points'][_time_frame]['right_cell']['coordinates']
        _pair_distance = compute.pair_distance_in_cell_size(
            _experiment=_experiment,
            _series_id=_series_id,
            _cell_1_coordinates=
            [(_left_cell_coordinates['x'], _left_cell_coordinates['y'], _left_cell_coordinates['z'])],
            _cell_2_coordinates=
            [(_right_cell_coordinates['x'], _right_cell_coordinates['y'], _right_cell_coordinates['z'])]
        )
        if _distance_range[0] <= _pair_distance <= _distance_range[1]:
            _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered


def by_pair_distances(_experiments_tuples, _distances, _range=1):
    # TODO: improve efficiency
    _experiments_tuples_filtered = []
    for _distance in _distances:
        _distance_tuples = by_pair_distance(_experiments_tuples, _distance, _range)
        for _tuple in _distance_tuples:
            if _tuple not in _experiments_tuples_filtered:
                _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered


def by_time_frames_amount(_experiments_tuples, _time_frames, _exactly=False):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        _time_frames_amount = len(_group_properties['time_points'])
        if (_exactly and _time_frames_amount == _time_frames) or \
                (not _exactly and _time_frames_amount >= _time_frames):
            _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered


def by_band(_experiments_tuples, _band=True):
    if _band is None:
        return _experiments_tuples

    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        _experiments_tuples_filtered.append(_tuple) if _group_properties['band'] == _band else None

    return _experiments_tuples_filtered


def by_series_id(_experiments_tuples, _series_id):
    return [_tuple for _tuple in _experiments_tuples if _tuple[1] == _series_id]


def by_real_pairs(_experiments_tuples, _real_pairs=True):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        if _group_properties['fake'] != _real_pairs and \
                ('real_fake' not in _group_properties or _group_properties['real_fake'] != _real_pairs):
            _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered


def by_real_fake_pairs(_experiments_tuples, _real_fake_pairs=True):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        if _real_fake_pairs and 'real_fake' in _group_properties and _group_properties['real_fake']:
            _experiments_tuples_filtered.append(_tuple)
        elif 'real_fake' not in _group_properties or not _group_properties['real_fake']:
            _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered


def by_fake_static_pairs(_experiments_tuples, _fake_static_pairs=True):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        if _group_properties['static'] == _fake_static_pairs:
            _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered


def by_main_cell(_experiments_tuples):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _cell_id = int(_group.split('_')[1])
        _series_properties = load.image_properties(_experiment, _series_id)
        if _cell_id == _series_properties['main_cell_id']:
            _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered


def by_triplets(_experiments_tuples):
    _experiments_triplets = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _ = _tuple
        # TODO: recognize how many cells and change the '10'
        for _cell_id_1 in range(0, 10):
            for _cell_id_2 in range(_cell_id_1 + 1, 10):
                _group_1 = 'cells_' + str(_cell_id_1) + '_' + str(_cell_id_2)
                for _cell_id_3 in range(_cell_id_2 + 1, 10):
                    _group_2 = 'cells_' + str(_cell_id_1) + '_' + str(_cell_id_3)
                    _group_3 = 'cells_' + str(_cell_id_2) + '_' + str(_cell_id_3)
                    if (_experiment, _series_id, _group_1) in _experiments_tuples and \
                            (_experiment, _series_id, _group_2) in _experiments_tuples and \
                            (_experiment, _series_id, _group_3) in _experiments_tuples:
                        _triplet = [
                            (_experiment, _series_id, _group_1),
                            (_experiment, _series_id, _group_2),
                            (_experiment, _series_id, _group_3)
                        ]
                        if _triplet not in _experiments_triplets:
                            _experiments_triplets.append([
                                (_experiment, _series_id, _group_1),
                                (_experiment, _series_id, _group_2),
                                (_experiment, _series_id, _group_3)
                            ])

    return _experiments_triplets


def by_dead_live(_experiments_tuples, _dead=None, _live=None):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_split = _group.split('_')
        _cell_id_1 = int(_group_split[1])
        _cell_id_2 = int(_group_split[2])
        _properties = load.image_properties(_experiment, _series_id)

        _cell_1_dead = _cell_id_1 in _properties['dead_cell_ids']
        _cell_2_dead = _cell_id_2 in _properties['dead_cell_ids']

        # both dead
        if _dead and not _live and _cell_1_dead and _cell_2_dead:
            _experiments_tuples_filtered.append(_tuple)
        # both live
        elif not _dead and _live and not _cell_1_dead and not _cell_2_dead:
            _experiments_tuples_filtered.append(_tuple)
        # at least one is dead
        elif _dead and _live is None and (_cell_1_dead or _cell_2_dead):
            _experiments_tuples_filtered.append(_tuple)
        # at least one is alive
        elif _dead is None and _live and (not _cell_1_dead or not _cell_2_dead):
            _experiments_tuples_filtered.append(_tuple)
        # one is dead, one is alive
        elif _dead and _live and ((_cell_1_dead and not _cell_2_dead) or (not _cell_1_dead and _cell_2_dead)):
            _experiments_tuples_filtered.append(_tuple)
        # return all
        elif _dead is None and _live is None:
            return _experiments_tuples

    return _experiments_tuples_filtered


def by_bleb_from_start(_experiments_tuples, _from_start=None):
    if _from_start is None:
        return _experiments_tuples

    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _, _ = _tuple
        if _from_start == (_experiment in BLEB_FROM_START):
            _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered


def by_bleb_amount_um(_experiments_tuples, _amount_um=None):
    if _amount_um is None:
        return _experiments_tuples

    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _, _ = _tuple
        if _amount_um == BLEB_AMOUNT_UM[_experiment]:
            _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered

from libs.experiments import compute, load, organize, filtering, config

SINGLE_CELLS = ['Single_Cell_Ortal']

CELL_PAIRS_REGULAR = ['SN16', 'LiveLive_201207', 'SN41', 'SN44', 'SN45']

CELL_PAIRS_BLEB_FROM_START = ['SN20_Bleb_fromStart', 'Bleb_150']
CELL_PAIRS_BLEB_IN_THE_MIDDLE = ['SN26_BlebAdded']
CELL_PAIRS_BLEB = CELL_PAIRS_BLEB_FROM_START + CELL_PAIRS_BLEB_IN_THE_MIDDLE

CELL_PAIRS_LIVE_DEAD = ['LiveDead_201117', 'LiveDead_201220']
CELL_PAIRS_DEAD_ONLY = ['DeadDead_201208']
CELL_PAIRS_WITH_DEAD_CELLS = CELL_PAIRS_LIVE_DEAD + CELL_PAIRS_DEAD_ONLY

CELL_PAIRS_BEADS = ['BeadBead_201209']

CELL_PAIRS_METASTASIS = ['Metastasis_Airyscan_210207']

CELL_PAIRS = CELL_PAIRS_REGULAR + CELL_PAIRS_BLEB + CELL_PAIRS_WITH_DEAD_CELLS + CELL_PAIRS_BEADS + \
             CELL_PAIRS_METASTASIS

EXPERIMENTS = SINGLE_CELLS + CELL_PAIRS

DISTANCE_RANGES = [(4, 6), (6, 8), (8, 10)]
DISTANCE_RANGE = (4, 10)


def get_type(_tuple):
    _experiment, _, _ = _tuple
    if _experiment in CELL_PAIRS_REGULAR:
        return 'Regular'
    elif _experiment in CELL_PAIRS_BLEB:
        return 'Bleb'
    elif _experiment in CELL_PAIRS_WITH_DEAD_CELLS:
        return 'With dead cells'
    elif _experiment in CELL_PAIRS_BEADS:
        return 'Beads'
    elif _experiment in CELL_PAIRS_METASTASIS:
        return 'Metastasis'
    else:
        raise Exception('No such experiment')


def get_pair_distance_range(_tuple):
    _experiment, _series_id, _group = _tuple
    try:
        _distance = compute.pair_distance_in_cell_size_time_frame(_experiment, _series_id, _group, _time_frame=0)
        if _distance < 4:
            return '<4'
        elif 4 <= _distance < 6:
            return '4-6'
        elif 6 <= _distance < 8:
            return '6-8'
        elif 8 <= _distance < 10:
            return '8-10'
        else:
            return '>10'
    except IndexError:
        return ''


def is_bleb_from_start(_tuple):
    _experiment, _, _ = _tuple
    return _experiment in CELL_PAIRS_BLEB_FROM_START


def bleb_amount(_tuple):
    _experiment, _, _ = _tuple
    if _experiment in config.BLEB_AMOUNT_UM:
        return config.BLEB_AMOUNT_UM[_experiment]
    return ''


def is_dead_cells_only(_tuple):
    _experiment, _, _ = _tuple
    return _experiment in CELL_PAIRS_DEAD_ONLY


def main():
    print(
        'Single cell / cell pair',
        'Regular / bleb / with dead cells / beads / metastasis',
        'Is real?',
        'Is band?',
        'Temporal resolution',
        'Pair distance',
        'Bleb from the beginning?',
        'Bleb amount',
        'is dead cells only?',
        'N',
        sep='\t'
    )

    _rows = []

    # single cells
    _single_cells = load.experiments_groups_as_tuples(SINGLE_CELLS)
    _single_cells = filtering.by_main_cell(_single_cells)
    _single_cells = organize.by_single_cell_id(_single_cells)
    print(
        'Single cell',
        'Regular',
        'Real',
        '',
        compute.temporal_resolution_in_minutes(SINGLE_CELLS[0]),
        '',
        '',
        '',
        '',
        len(_single_cells),
        sep='\t'
    )

    # cell pairs
    _cell_pairs = load.experiments_groups_as_tuples(CELL_PAIRS)
    for _tuple in _cell_pairs:
        _experiment, _series_id, _group = _tuple
        print(
            'Cell pair',
            get_type(_tuple),
            filtering.is_real(_tuple),
            filtering.is_band(_tuple),
            compute.temporal_resolution_in_minutes(_experiment),
            get_pair_distance_range(_tuple),
            is_bleb_from_start(_tuple),
            bleb_amount(_tuple),
            is_dead_cells_only(_tuple),
            1,
            sep='\t'
        )


if __name__ == '__main__':
    main()

import pandas as pd
from tqdm import tqdm

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
            return 'Less than 4'
        elif 4 <= _distance < 6:
            return 'Between 4 and 6'
        elif 6 <= _distance < 8:
            return 'Between 6 and 8'
        elif 8 <= _distance < 10:
            return 'Between 8 and 10'
        else:
            return 'More than 10'
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


def print_df(_df):
    print('\t'.join(map(str, _df.columns)))
    for _i in range(len(_df.index)):
        _row = list(_df.iloc[_i])
        print('\t'.join(map(str, _row)))


def main():
    _headers = [
        'Single cell / cell pair',
        'Regular / bleb / with dead cells / beads / metastasis',
        'Is real?',
        'Is band?',
        'Temporal resolution',
        'Pair distance',
        'Bleb from the beginning?',
        'Bleb amount',
        'is dead cells only?',
        'N'
    ]

    _rows = []

    # single cells
    print('Single cells')
    _single_cells = load.experiments_groups_as_tuples(SINGLE_CELLS)
    _single_cells = filtering.by_main_cell(_single_cells)
    _single_cells = organize.by_single_cell_id(_single_cells)
    _rows.append([
        'Single cell',
        'Regular',
        True,
        '',
        compute.temporal_resolution_in_minutes(SINGLE_CELLS[0]),
        '',
        '',
        '',
        '',
        len(_single_cells)
    ])

    # cell pairs
    _cell_pairs = load.experiments_groups_as_tuples(CELL_PAIRS)
    for _tuple in tqdm(_cell_pairs, desc='Cell pairs'):
        _experiment, _series_id, _group = _tuple
        _rows.append([
            'Cell pair',
            get_type(_tuple),
            filtering.is_real(_tuple),
            filtering.is_band(_tuple),
            compute.temporal_resolution_in_minutes(_experiment),
            get_pair_distance_range(_tuple),
            is_bleb_from_start(_tuple),
            bleb_amount(_tuple),
            is_dead_cells_only(_tuple),
            1
        ])

    _df = pd.DataFrame(data=_rows, columns=_headers)

    # summarize
    _summarize = _df.groupby(by=[
        'Single cell / cell pair',
        'Regular / bleb / with dead cells / beads / metastasis',
        'Is real?',
        'Is band?',
        'Temporal resolution',
        'Pair distance',
        'Bleb from the beginning?',
        'Bleb amount',
        'is dead cells only?'
    ]).sum()

    # print(_summarize.reset_index())
    print_df(_summarize.reset_index())


if __name__ == '__main__':
    main()

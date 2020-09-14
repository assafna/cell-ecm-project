from libs.experiments import load, organize, filtering

SINGLE_CELLS = ['Single_Cell_Ortal']
CELL_PAIRS_EXPERIMENTS_LOW_TEMPORAL_RESOLUTION = ['SN16']
CELL_PAIRS_EXPERIMENTS_HIGH_TEMPORAL_RESOLUTION = ['SN41', 'SN44', 'SN45']
CELL_PAIRS = CELL_PAIRS_EXPERIMENTS_LOW_TEMPORAL_RESOLUTION + CELL_PAIRS_EXPERIMENTS_HIGH_TEMPORAL_RESOLUTION
DISTANCE_RANGES = [(4, 6), (6, 8), (8, 10)]


def main():
    # single cells
    print('Single cells')
    _single_cells = load.experiments_groups_as_tuples(SINGLE_CELLS)
    _single_cells = filtering.by_main_cell(_single_cells)
    _single_cells = organize.by_single_cell_id(_single_cells)
    print('\tTotal single cell experiments:', len(_single_cells))

    # cell pairs
    _cell_pairs = load.experiments_groups_as_tuples(CELL_PAIRS)
    print('\nTotal cell pairs:', len(_cell_pairs))

    for _distance_range in DISTANCE_RANGES:
        print('Cell pairs distance range:', _distance_range)

        # total
        _distance_range_cell_pairs = filtering.by_pair_distance_range(_cell_pairs, _distance_range=_distance_range)
        print('\tTotal cell pairs:', len(_distance_range_cell_pairs))

        for _experiment in CELL_PAIRS:
            print('\tExperiment:', _experiment)

            _experiment_cell_pairs = load.experiment_groups_as_tuples(_experiment)
            _distance_range_experiment_cell_pairs = \
                filtering.by_pair_distance_range(_experiment_cell_pairs, _distance_range=_distance_range)
            print('\t\tTotal cell pairs:', len(_distance_range_experiment_cell_pairs))

            # real cell pairs
            _real_cell_pairs = filtering.by_real_pairs(_distance_range_experiment_cell_pairs)
            print('\t\tTotal real cell pairs:', len(_real_cell_pairs))

            _real_cell_pairs_with_band = filtering.by_band(_real_cell_pairs)
            print('\t\tTotal real cell pairs with band:', len(_real_cell_pairs_with_band))

            # fake cell pairs
            _fake_cell_pairs = filtering.by_real_pairs(_distance_range_experiment_cell_pairs, _real_pairs=False)
            print('\t\tTotal fake cell pairs:', len(_fake_cell_pairs))

            _fake_static_cell_pairs = filtering.by_fake_static_pairs(_fake_cell_pairs)
            print('\t\tTotal fake static cell pairs:', len(_fake_static_cell_pairs))

            _fake_following_cell_pairs = filtering.by_fake_static_pairs(_fake_cell_pairs, _fake_static_pairs=False)
            print('\t\tTotal fake following cell pairs:', len(_fake_following_cell_pairs))

            _fake_following_cell_pairs_with_band = filtering.by_band(_fake_following_cell_pairs)
            print('\t\tTotal fake following cell pairs with band:', len(_fake_following_cell_pairs_with_band))


if __name__ == '__main__':
    main()

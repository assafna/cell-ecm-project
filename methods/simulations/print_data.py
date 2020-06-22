from itertools import product

from libs.simulations import load, filtering

DISTANCES = [5, 7, 9]
STDS = [0, 0.25, 0.5, 0.75]
ALPHAS = [0, 0.25, 0.5, 0.75, 1]
BETAS = [1, 1.05, 1.1, 1.2]


def main():
    _simulations = load.structured()

    # single cell
    print('Single cells')
    _single_cell_simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=True,
        _is_heterogeneity=None,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    print('\tTotal single cell simulations:', len(_single_cell_simulations))

    # heterogeneity
    _heterogeneity_single_cell_simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=True,
        _is_heterogeneity=True,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    print('\tTotal single cell heterogeneity:', len(_heterogeneity_single_cell_simulations))

    for _std in STDS:
        _std_simulations = filtering.by_heterogeneity(_heterogeneity_single_cell_simulations, _std=_std)
        print('\t\tTotal with std. ' + str(_std) + ':', len(_std_simulations))

    # cell pairs
    _cell_pairs_simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=None,
        _is_low_connectivity=False,
        _is_causality=None,
        _is_dominant_passive=False
    )

    for _distance in DISTANCES:
        print('\nCell pairs distance:', _distance)

        # total
        _distance_simulations = filtering.by_distance(_cell_pairs_simulations, _distance=_distance)
        print('\tTotal cell pairs simulations:', len(_distance_simulations))

        # no heterogeneity
        _no_heterogeneity_simulations = filtering.by_categories(
            _distance_simulations,
            _is_single_cell=False,
            _is_heterogeneity=False,
            _is_low_connectivity=False,
            _is_causality=False,
            _is_dominant_passive=False
        )
        print('\tTotal cell pairs no heterogeneity (no causality, no dominant):', len(_no_heterogeneity_simulations))

        # heterogeneity without causality
        _filtered_simulations = filtering.by_categories(
            _distance_simulations,
            _is_single_cell=False,
            _is_heterogeneity=True,
            _is_low_connectivity=False,
            _is_causality=False,
            _is_dominant_passive=False
        )
        print('\tTotal cell pairs heterogeneity (no causality, no dominant):', len(_filtered_simulations))

        for _std in STDS:
            _std_simulations = filtering.by_heterogeneity(_filtered_simulations, _std=_std)
            print('\t\tTotal with std. ' + str(_std) + ':', len(_std_simulations))

        # heterogeneity with causality
        _filtered_simulations = filtering.by_categories(
            _distance_simulations,
            _is_single_cell=False,
            _is_heterogeneity=True,
            _is_low_connectivity=False,
            _is_causality=True,
            _is_dominant_passive=False
        )
        print('\tTotal cell pairs heterogeneity (causality):', len(_filtered_simulations))

        for _alpha, _beta in product(ALPHAS, BETAS):
            _causality_simulations = filtering.by_causality(_filtered_simulations, _alpha=_alpha, _beta=_beta)
            print('\t\tTotal with alpha = ' + str(_alpha) + ', beta = ' + str(_beta) + ':', len(_causality_simulations))


if __name__ == '__main__':
    main()

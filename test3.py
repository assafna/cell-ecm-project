from libs.simulations import load, filtering

if __name__ == '__main__':
    _simulations = load.structured()
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=True,
        _is_low_connectivity=True,
        _is_causality=False,
        _is_dominant_passive=False
    )
    print(len(_simulations))
    print(_simulations)

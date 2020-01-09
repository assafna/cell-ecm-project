from libs.experiments import load, config, filtering


def main():
    _experiments = load.fibers_density_dictionary(config.PAIRS)
    _filtered = filtering.by_distance(_experiments, _min_distance=5, _max_distance=7)
    _filtered_by_time = filtering.by_time_points_amount(_filtered, _time_points=10)
    print('hi')


if __name__ == '__main__':
    main()

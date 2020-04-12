import fibers_density


def figure_1():
    print('Panel C-ii')
    fibers_density.fibers_vs_distance_single_cells.main()

    print('Panel D-ii')
    fibers_density.fibers_vs_distance_pairs.main()

    print('Panel D-iii')
    fibers_density.fibers_vs_distance_differences.main()

    print('Panel E-i')
    fibers_density.simulations.fibers_vs_distance_multiple_cell_distances.main()

    print('Panel E-ii')
    fibers_density.experiments.fibers_vs_distance_multiple_cell_distances.main()


def figure_si_1_correlation_distance_density_offset():
    print('Panel B')
    fibers_density.experiments.fibers_and_pairs_cells_distance_correlations_vs_distance.main()

    print('Panel D')
    fibers_density.experiments.fibers_vs_pairs_cells_distance_in_offset.main()


def figure_2():
    print('Panel A-ii')
    fibers_density.fibers_in_time_single_cells.main()

    print('Panel B-ii')
    fibers_density.fibers_in_time_pairs.main()

    # print('Panel B-iii')
    # fibers_density.reaching_std_multiple_cell_distances.main()

    print('Panel C-i')
    fibers_density.simulations.fibers_vs_change.main()

    print('Panel C-ii')
    fibers_density.experiments.fibers_vs_change.main()


def figure_3():
    print('Panel B')
    fibers_density.simulations.correlations_by_derivatives_pairs.main()

    print('Panels C & D')
    fibers_density.simulations.insides_vs_outsides_stds.main()


def figure_si_3_stationary_and_detrending():
    print('Panel A')
    fibers_density.simulations.fibers_vs_time_derivatives.main()

    print('Panel B')
    fibers_density.simulations.stationary_vs_fibers_derivatives_pairs.main()


def figure_si_3_single_cell_detrending():
    print('Panel A-i')
    fibers_density.simulations.stationary_vs_fibers_derivatives_single_cells.main()

    print('Panel A-ii')
    fibers_density.simulations.correlations_by_derivative_single_cells.main()


def figure_4():
    fibers_density.simulations.communicated_vs_non_communicated.main()
    fibers_density.experiments.communicated_vs_non_communicated.main()
    fibers_density.simulations.communicated_vs_non_communicated_cell_distances.main()


if __name__ == '__main__':
    print('Figure 1')
    figure_1()

    print('Figure SI 1 - Correlation Distance Density Offset')
    figure_si_1_correlation_distance_density_offset()

    print('Figure 2')
    figure_2()

    print('Figure 3')
    figure_3()

    print('Figure SI 3 - Stationary & Detrending')
    figure_si_3_stationary_and_detrending()

    print('Figure SI 3 - Single Cell Detrending')
    figure_si_3_single_cell_detrending()

    print('Figure 4')
    figure_4()

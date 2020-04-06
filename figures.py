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


def figure_1_supp_correlation_distance_density_offset():
    print('Panel B')
    fibers_density.experiments.fibers_and_pairs_cells_distance_correlations_vs_distance.main()

    print('Panel D-i')
    fibers_density.experiments.fibers_vs_pairs_cells_distance_in_offset.main(_offset_x=1)

    print('Panel D-ii')
    fibers_density.experiments.fibers_vs_pairs_cells_distance_in_offset.main(_offset_x=2.5)


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
    print('?')
    fibers_density.simulations.stationary_vs_fibers_derivatives.main()

    print('?')
    fibers_density.simulations.insides_vs_outsides_derivatives.main()

    print('?')
    fibers_density.experiments.insides_vs_outsides_derivatives.main()

    # fibers_density.simulations.fibers_vs_derivatives_single_cell.main()
    # fibers_density.simulations.insides_vs_outsides_stds.main()
    # fibers_density.simulations.insides_vs_outsides_cell_distances.main()


def figure_4():
    fibers_density.simulations.communicated_vs_non_communicated.main()
    fibers_density.experiments.communicated_vs_non_communicated.main()
    fibers_density.simulations.communicated_vs_non_communicated_cell_distances.main()


if __name__ == '__main__':
    print('Figure 1')
    figure_1()

    print('Figure 1 Supporting Information - Correlation Distance Density Offset')
    figure_1_supp_correlation_distance_density_offset()

    print('Figure 2')
    figure_2()

    print('Figure 3')
    figure_3()

    print('Figure 4')
    figure_4()

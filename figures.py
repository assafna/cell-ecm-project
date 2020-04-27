import fibers_density


def figure_1():
    print('Panel C-ii')
    fibers_density.fibers_vs_distance_single_cells.main(_low_connectivity=False)

    print('Panel D-ii')
    fibers_density.fibers_vs_distance_pairs.main(_low_connectivity=False)

    print('Panel D-iii')
    fibers_density.fibers_vs_distance_differences.main(_low_connectivity=False)

    print('Panel E-i')
    fibers_density.simulations.fibers_vs_distance_multiple_cell_distances.main(_low_connectivity=False)

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

    print('Panel C-i')
    fibers_density.simulations.fibers_vs_change.main()

    print('Panel C-ii')
    fibers_density.experiments.fibers_vs_change.main(_early_time_points=True)


def figure_si_2_density_vs_change_in_density_in_late_time_points():
    print('Main panel')
    fibers_density.experiments.fibers_vs_change.main(_early_time_points=False)


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
    print('Panel A & B')
    fibers_density.simulations.stationary_vs_fibers_derivatives_single_cells.main()

    print('Panel C')
    fibers_density.simulations.correlations_by_derivatives_single_cells.main()


def figure_si_3_insides_vs_outsides_distances():
    print('Panels A & B')
    fibers_density.simulations.insides_vs_outsides_distances.main()


def figure_si_3_insides_vs_outsides_offsets():
    print('Panels A & B')
    fibers_density.simulations.insides_vs_outsides_offsets.main()


def figure_4():
    print('Panel B-i')
    fibers_density.simulations.same_vs_different.main(_low_connectivity=False)

    print('Panel B-ii')
    fibers_density.simulations.same_vs_different_distances.main(_low_connectivity=False)

    print('Panel C-i')
    fibers_density.experiments.same_vs_different.main(_band=True, _high_time_resolution=False)

    print('Panel C-ii')
    fibers_density.experiments.same_vs_different_distances.main(_band=True, _high_time_resolution=False)

    print('Panel D-ii')
    fibers_density.experiments.same_vs_different_triplets.main()

    print('Panel D-iii')
    fibers_density.experiments.same_vs_different_triplet.main()

    print('Panel F')
    fibers_density.experiments.same_vs_different.main(_band=False, _high_time_resolution=False)

    print('?')
    fibers_density.experiments.fibers_vs_time.compute_tuples([('SN16', 3, 'cells_0_1')])

    print('?')
    fibers_density.experiments.same_vs_different_compare.main(_high_time_resolution=False)


def figure_si_4_low_connectivity():
    print('?')


def figure_si_4_single_cell_detrending():
    print('Panel A & B')
    fibers_density.experiments.stationary_vs_fibers_derivatives_single_cells.main()

    print('Panel C')
    fibers_density.experiments.correlations_by_derivatives_single_cells.main()


def figure_si_4_temporal_assessment():
    print('?')
    fibers_density.experiments.same_vs_different_temporal_reduction.main()


def figure_5():
    print('Panel B-i')
    fibers_density.experiments.fibers_vs_offsets_in_axes.main()


def figure_si_5():
    print('?')
    fibers_density.experiments.fibers_vs_same_vs_different_offsets_in_axes.main()


if __name__ == '__main__':
    print('Figure 1')
    figure_1()

    print('Figure SI 1 - Correlation Distance Density Offset')
    figure_si_1_correlation_distance_density_offset()

    print('Figure 2')
    figure_2()

    print('Figure SI 2 - Density vs. Change in Density in Late Time-points')
    figure_si_2_density_vs_change_in_density_in_late_time_points()

    print('Figure 3')
    figure_3()

    print('Figure SI 3 - Stationary & Detrending')
    figure_si_3_stationary_and_detrending()

    print('Figure SI 3 - Single Cell Detrending')
    figure_si_3_single_cell_detrending()

    print('Figure SI 3 - Insides vs. Outsides Distances')
    figure_si_3_insides_vs_outsides_distances()

    print('Figure SI 3 - Insides vs. Outsides Offsets')
    figure_si_3_insides_vs_outsides_offsets()

    print('Figure 4')
    figure_4()

    print('Figure SI 4 - Low Connectivity')
    figure_si_4_low_connectivity()

    print('Figure SI 4 - Single Cell Detrending')
    figure_si_4_single_cell_detrending()

    print('Figure 5')
    figure_5()

    print('Figure SI 5 - ?')
    figure_si_5()

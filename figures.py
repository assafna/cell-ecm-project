import fibers_density


def figure_1_main():
    print('\nPanel D')
    fibers_density.fibers_vs_distance_single_cells.main(_low_connectivity=False)

    print('\nPanel F')
    fibers_density.fibers_vs_distance_pairs.main(_low_connectivity=False)

    print('\nPanel G')
    fibers_density.fibers_vs_distance_differences.main(_low_connectivity=False)


def figure_si_1_distances():
    print('\nPanel A')
    fibers_density.simulations.fibers_vs_distance_multiple_cell_distances.main(_low_connectivity=False)

    print('\nPanel B')
    fibers_density.experiments.fibers_vs_distance_multiple_cell_distances.main()


def figure_si_1_correlation_distance_density_offset():
    print('\nPanels B (main panel) & C')
    fibers_density.experiments.fibers_and_pairs_cells_distance_correlations_vs_distance.main()

    print('\nPanel B (secondary panels)')
    fibers_density.experiments.fibers_vs_pairs_cells_distance_in_offset.main()


def figure_1_all():
    print('\nFigure 1 - Main Figure')
    figure_1_main()

    print('\nFigure SI 1')
    print('\nFigure SI 1 - Distances')
    figure_si_1_distances()

    print('\nFigure SI 1 - Correlation Distance Density Offset')
    figure_si_1_correlation_distance_density_offset()


def figure_2_main():
    print('\nPanel B')
    fibers_density.fibers_in_time_single_cells.main()

    print('\nPanel D')
    fibers_density.fibers_in_time_pairs.main()

    print('\nPanel E')
    fibers_density.simulations.fibers_vs_change.main()

    print('\nPanel F')
    fibers_density.experiments.fibers_vs_change.main(_early_time_points=True)


def figure_si_2_density_vs_change_in_density_in_late_time_points():
    fibers_density.experiments.fibers_vs_change.main(_early_time_points=False)


def figure_2_all():
    print('\nFigure 2 - Main Figure')
    figure_2_main()

    print('\nFigure SI 2')
    print('\nFigure SI 2 - Density vs. Change in Density in Late Time-points')
    figure_si_2_density_vs_change_in_density_in_late_time_points()


def figure_3_main():
    print('\nPanels B & C')
    fibers_density.simulations.correlations_by_derivatives_pairs.main()

    print('\nPanels D & E')
    fibers_density.simulations.insides_vs_outsides_stds.main()


def figure_si_3_detrending_simulations():
    print('\nPanels A, B & C')
    fibers_density.simulations.fibers_vs_time_derivatives.main()

    print('\nPanel D')
    fibers_density.simulations.correlations_by_derivatives_pairs.main(_directions=['inside'])

    print('\nPanels E & F')
    fibers_density.simulations.stationary_vs_fibers_derivatives_pairs.main()

    print('\nPanel G')
    fibers_density.simulations.correlations_by_derivatives_single_cells.main()

    print('\nPanels H & I')
    fibers_density.simulations.stationary_vs_fibers_derivatives_single_cells.main()


def figure_si_3_insides_vs_outsides_distances():
    print('\nPanels A & B')
    fibers_density.simulations.insides_vs_outsides_distances.main()


def figure_si_3_insides_vs_outsides_offsets():
    print('\nPanels A & B')
    fibers_density.simulations.insides_vs_outsides_offsets.main()


def figure_3_all():
    print('\nFigure 3 - Main Figure')
    figure_3_main()

    print('\nFigure SI 3')
    print('\nFigure SI 3 - Detrending Simulations')
    figure_si_3_detrending_simulations()

    print('\nFigure SI 3 - Insides vs. Outsides Distances')
    figure_si_3_insides_vs_outsides_distances()

    print('\nFigure SI 3 - Insides vs. Outsides Offsets')
    figure_si_3_insides_vs_outsides_offsets()


def figure_4_main():
    print('\nPanel B')
    fibers_density.simulations.same_vs_different.main(_low_connectivity=False)

    print('\nPanel C')
    fibers_density.simulations.same_vs_different_distances.main(_low_connectivity=False)

    print('\nPanel D')
    fibers_density.experiments.same_vs_different.main(_band=True, _high_time_resolution=False)

    print('\nPanel E')
    fibers_density.experiments.same_vs_different_distances.main(_band=True, _high_time_resolution=False)

    print('\nPanel G')
    fibers_density.experiments.same_vs_different_triplet.main()

    print('\nPanel H')
    fibers_density.experiments.same_vs_different_triplets.main()


def figure_si_4_same_vs_diff_offsets():
    print('\nPanel A')
    fibers_density.experiments.same_vs_different.main(_offset_y=0, _offset_z=0)

    print('\nPanel B')
    fibers_density.experiments.same_vs_different.main(_offset_y=-0.5, _offset_z=0)

    print('\nPanel C')
    fibers_density.experiments.same_vs_different.main(_offset_y=0.5, _offset_z=0)

    print('\nPanel D')
    fibers_density.experiments.same_vs_different.main(_offset_y=0, _offset_z=-0.5)

    print('\nPanel E')
    fibers_density.experiments.same_vs_different.main(_offset_y=0, _offset_z=0.5)


def figure_si_4_detrending_experiments():
    print('\nPanel A')
    fibers_density.experiments.correlations_by_derivatives_pairs.main(_directions=['inside'])

    print('\nPanels B & C')
    fibers_density.experiments.stationary_vs_fibers_derivatives_pairs.main()

    print('\nPanel D')
    fibers_density.experiments.correlations_by_derivatives_single_cells.main()

    print('\nPanels E & F')
    fibers_density.experiments.stationary_vs_fibers_derivatives_single_cells.main()


def figure_si_4_same_vs_diff_distances():
    print('\nPanel A')
    fibers_density.experiments.same_vs_different.main(
        _band=True, _high_time_resolution=False, _cells_distance_range=[4, 6])

    print('\nPanel B')
    fibers_density.experiments.same_vs_different.main(
        _band=True, _high_time_resolution=False, _cells_distance_range=[6, 8])

    print('\nPanel C')
    fibers_density.experiments.same_vs_different.main(
        _band=True, _high_time_resolution=False, _cells_distance_range=[8, 10])


def figure_si_4_same_vs_diff_controls():
    print('\nPanel B')
    fibers_density.experiments.same_vs_different.main(
        _real_cells=False, _static=False, _band=True, _high_time_resolution=False)

    print('\nPanel D')
    fibers_density.experiments.same_vs_different.main(
        _real_cells=False, _static=True, _band=False, _high_time_resolution=False)


def figure_si_4_same_vs_diff_high_time_res():
    print('\nPanel A')
    fibers_density.experiments.same_vs_different.main(_band=True, _high_time_resolution=True)

    print('\nPanel B')
    fibers_density.experiments.same_vs_different.main(_band=False, _high_time_resolution=True)


def figure_si_4_temporal_assessment():
    print('\nPanel A')
    fibers_density.experiments.same_vs_different_temporal_assessment.main(_high_time_resolution=False)

    print('\nPanel B')
    fibers_density.experiments.same_vs_different_temporal_assessment.main(_high_time_resolution=True)


def figure_4_all():
    print('\nFigure 4 - Main Figure')
    figure_4_main()

    print('\nFigure SI 4')
    print('\nFigure SI 4 - Same vs. Different - Offsets')
    figure_si_4_same_vs_diff_offsets()

    print('\nFigure SI 4 - Detrending Experiments')
    figure_si_4_detrending_experiments()

    print('\nFigure SI 4 - Same vs. Different - Distances')
    figure_si_4_same_vs_diff_distances()

    print('\nFigure SI 4 - Same vs. Different - Controls')
    figure_si_4_same_vs_diff_controls()

    print('\nFigure SI 4 - Real vs. Fake')
    figure_si_4_real_vs_fake()

    print('\nFigure SI 4 - Same vs. Different - High Temporal Resolution')
    figure_si_4_same_vs_diff_high_time_res()

    print('\nFigure SI 4 - Temporal Assessment')
    figure_si_4_temporal_assessment()


def figure_5_main():
    print('\nPanel A')
    fibers_density.experiments.fibers_vs_offsets_in_axes.main(_high_time_resolution=False)

    print('\nPanel B')
    fibers_density.experiments.same_vs_different_offsets_in_axes.main(_high_time_resolution=False)

    print('\nPanel C')
    fibers_density.experiments.density_vs_communication_vs_z_offset.main()

    print('\nPanel F')
    fibers_density.experiments.fibers_vs_time.compute_tuples([('SN16', 3, 'cells_0_1')])

    print('\nPanel G')
    fibers_density.experiments.same_vs_different.main(_band=False, _high_time_resolution=False)

    print('\nPanel H')
    fibers_density.experiments.same_vs_different_band_vs_no_band.main(_high_time_resolution=False)


def figure_si_5_high_time_res():
    print('\nPanel A')
    fibers_density.experiments.fibers_vs_offsets_in_axes.main(_high_time_resolution=True)

    print('\nPanel B')
    fibers_density.experiments.same_vs_different_offsets_in_axes.main(_high_time_resolution=True)


def figure_si_5_simulations():
    print('\nPanel A')
    fibers_density.simulations.fibers_vs_offsets_in_axes.main(_low_connectivity=False)

    print('\nPanel B')
    fibers_density.simulations.same_vs_different_offsets_in_axes.main(_low_connectivity=False)


def figure_si_5_no_band():
    print('\nPanel A')
    fibers_density.experiments.fibers_vs_offsets_in_axes.main(_band=False, _high_time_resolution=False)

    print('\nPanel B')
    fibers_density.experiments.same_vs_different_offsets_in_axes.main(_band=False, _high_time_resolution=False)


def figure_si_4_real_vs_fake():
    print('\nPanels A & B')
    fibers_density.experiments.same_vs_different_real_vs_fake.main(_offset_y=0, _high_time_resolution=False)

    print('\nPanels C & D')
    fibers_density.experiments.same_vs_different_real_vs_fake.main(_offset_y=0.5, _high_time_resolution=False)

    print('\nPanels E & F')
    fibers_density.experiments.same_vs_different_real_vs_fake.main(_offset_y=0, _high_time_resolution=True)

    print('\nPanels G & H')
    fibers_density.experiments.same_vs_different_real_vs_fake.main(_offset_y=0.5, _high_time_resolution=True)


def figure_si_5_communication_vs_density():
    fibers_density.experiments.fibers_vs_same_vs_different_offsets_in_axes.main()


def figure_5_all():
    print('\nFigure 5 - Main Figure')
    figure_5_main()

    print('\nFigure SI 5')
    print('\nFigure SI 5 - High Temporal Resolution')
    figure_si_5_high_time_res()

    # print('\nFigure SI 5 - Simulations')
    # figure_si_5_simulations()

    print('\nFigure SI 5 - No Band')
    figure_si_5_no_band()

    print('\nFigure SI 5 - Communication vs. Density')
    figure_si_5_communication_vs_density()


def figure_6_main():
    print('\nPanels B & C')
    fibers_density.experiments.matchmaking_by_correlation.main(
        _real_cells=True, _offset_y=0.5, _high_time_resolution=False)

    print('\nPanel D')
    fibers_density.experiments.matchmaking_by_correlation.main(
        _real_cells=True, _offset_y=0, _high_time_resolution=False)

    print('\nPanel E')
    fibers_density.experiments.matchmaking_by_correlation.main(
        _real_cells=False, _offset_y=0.5, _high_time_resolution=False)

    print('\nPanel F')
    fibers_density.experiments.matchmaking_by_correlation.main(
        _real_cells=True, _offset_y=0.5, _high_time_resolution=True)


def figure_6_all():
    print('\nFigure 6 - Main Figure')
    figure_6_main()


def figure_7_main():
    print('\nPanel B')
    fibers_density.simulations.same_vs_different_cross_correlation.main(
        _alpha=1, _beta=1, _low_connectivity=False, _plots=['same'], _plot_types=['box'])

    print('\nPanel C')
    fibers_density.simulations.same_vs_different_cross_correlation_alphas_betas.main(
        _type='alpha', _low_connectivity=False, _plots=['same'], _plot_types=['box'])

    print('\nPanel D')
    fibers_density.simulations.same_vs_different_cross_correlation_alphas_betas.main(
        _type='beta', _low_connectivity=False, _plots=['same'], _plot_types=['box'])

    print('\nPanel E')
    fibers_density.experiments.inner_regions_granger_causality.main(
        _band=None, _high_time_resolution=True, _plots=['granger'])


def figure_si_7_simulate_time_lag_distribution():
    print('\nPanel A')
    fibers_density.simulations.same_vs_different_cross_correlation.main(
        _alpha=1, _beta=1, _low_connectivity=False, _plots=['same'], _plot_types=['bar'])

    print('\nPanel B')
    fibers_density.simulations.same_vs_different_cross_correlation_alphas_betas.main(
        _type='alpha', _low_connectivity=False, _plots=['same'], _plot_types=['bar'])

    print('\nPanel C')
    fibers_density.simulations.same_vs_different_cross_correlation_alphas_betas.main(
        _type='beta', _low_connectivity=False, _plots=['same'], _plot_types=['bar'])


def figure_si_7_leader_follower_cross_correlation_experiments():
    fibers_density.experiments.same_vs_different_cross_correlation.main(
        _plots=['same', 'different'], _plot_types=['box', 'bar'])


def figure_si_7_different_network():
    print('\nPanels A & B')
    fibers_density.simulations.same_vs_different_cross_correlation.main(
        _alpha=1, _beta=1, _low_connectivity=False, _plots=['different'])

    print('\nPanels C & D')
    fibers_density.simulations.same_vs_different_cross_correlation_alphas_betas.main(
        _type='alpha', _low_connectivity=False, _plots=['different'])

    print('\nPanels E & F')
    fibers_density.simulations.same_vs_different_cross_correlation_alphas_betas.main(
        _type='beta', _low_connectivity=False, _plots=['different'])


def figure_si_7_granger_causality():
    print('\nPanels B & C')
    fibers_density.experiments.inner_regions_granger_causality.main(
        _band=None, _high_time_resolution=True, _tuples_to_plot=[('SN45', 2, 'cells_2_3')])

    print('\nPanels E & F')
    fibers_density.experiments.inner_regions_granger_causality.main(
        _band=None, _high_time_resolution=True, _tuples_to_plot=[('SN45', 1, 'cells_1_3')])


def figure_si_7_whiteness():
    fibers_density.experiments.inner_regions_granger_causality.main(
        _band=None, _high_time_resolution=True, _plots=['whiteness'])


def figure_7_all():
    print('\nFigure 7 - Main Figure')
    figure_7_main()

    print('\nFigure SI 7')
    print('\nFigure SI 7 - Simulate Time Lag Distribution')
    figure_si_7_simulate_time_lag_distribution()

    print('\nFigure SI 7 - Leader Follower Cross Correlation Experiments')
    figure_si_7_leader_follower_cross_correlation_experiments()

    print('\nFigure SI 7 - Different Network')
    figure_si_7_different_network()

    print('\nFigure SI 7 - Granger Causality')
    figure_si_7_granger_causality()

    print('\nFigure 7 SI - Whiteness')
    figure_si_7_whiteness()


def figure_si_discussion_saturation():
    print('\nPanel A')
    fibers_density.experiments.inner_regions_saturation.main(_offset_y=0)

    print('\nPanel B')
    fibers_density.experiments.inner_regions_saturation.main(_offset_y=0.5)


def figure_discussion_all():
    print('\nFigure SI Discussion - Saturation')
    figure_si_discussion_saturation()


if __name__ == '__main__':
    print('\nFigures')
    figure_1_all()
    figure_2_all()
    figure_3_all()
    figure_4_all()
    figure_5_all()
    figure_6_all()
    figure_7_all()
    figure_discussion_all()

import fiber_density


def figure_1_main():
    print('\nPanel D')
    fiber_density.inner_density_vs_window_distance_single_cells.main(_low_connectivity=False)

    print('\nPanel F')
    fiber_density.inner_density_vs_window_distance_cell_pairs.main(_low_connectivity=False)

    print('\nPanel G')
    fiber_density.inner_density_vs_window_distance_differences.main(_low_connectivity=False)


def figure_si_1_distances():
    print('\nPanel A')
    fiber_density.simulations.inner_density_vs_window_distance.main(_low_connectivity=False)

    print('\nPanel B')
    fiber_density.experiments.inner_density_vs_window_distance.main()


def figure_si_1_correlation_distance_density_offset():
    print('\nPanels B (main panel) & C')
    fiber_density.experiments.inner_density_vs_pair_distance_correlations_vs_window_distance.main()

    print('\nPanel B (secondary panels)')
    fiber_density.experiments.inner_density_vs_pair_distance.main()


def figure_1_all():
    print('\nFigure 1 - Main Figure')
    figure_1_main()

    print('\nFigure SI 1')
    print('\nFigure SI 1 - Distances (S2)')
    figure_si_1_distances()

    print('\nFigure SI 1 - Correlation Distance Density Offset (S3)')
    figure_si_1_correlation_distance_density_offset()


def figure_2_main():
    print('\nPanel B')
    fiber_density.inner_density_vs_time_single_cells.main()

    print('\nPanel D')
    fiber_density.inner_density_vs_time_cell_pairs.main()

    print('\nPanel E')
    fiber_density.simulations.inner_density_vs_change_in_inner_density.main()

    print('\nPanel F')
    fiber_density.experiments.inner_density_vs_change_in_inner_density.main(_early_time_frames=True)


def figure_si_2_density_vs_change_in_density_in_late_time_points():
    fiber_density.experiments.inner_density_vs_change_in_inner_density.main(_early_time_frames=False)


def figure_2_all():
    print('\nFigure 2 - Main Figure')
    figure_2_main()

    print('\nFigure SI 2')
    print('\nFigure SI 2 - Density vs. Change in Density in Late Time-points (S4)')
    figure_si_2_density_vs_change_in_density_in_late_time_points()


def figure_3_main():
    print('\nPanels B & C')
    fiber_density.simulations.inner_correlation_by_derivatives_cell_pairs.main()

    print('\nPanels D & E')
    fiber_density.simulations.inner_correlation_vs_outer_correlation_heterogeneity.main()


def figure_si_3_detrending_simulations():
    print('\nPanels A, B & C')
    fiber_density.simulations.inner_density_vs_time_derivatives.main()

    print('\nPanel D')
    fiber_density.simulations.inner_correlation_by_derivatives_cell_pairs.main(_directions=['inside'])

    print('\nPanels E & F')
    fiber_density.simulations.stationary_vs_inner_dynamics_derivatives_cell_pairs.main()

    print('\nPanel G')
    fiber_density.simulations.inner_correlation_by_derivatives_single_cells.main()

    print('\nPanels H & I')
    fiber_density.simulations.stationary_vs_inner_dynamics_derivatives_single_cells.main()


def figure_si_3_insides_vs_outsides_distances():
    print('\nPanels A & B')
    fiber_density.simulations.inner_correlation_vs_outer_correlation_pair_distance.main()


def figure_si_3_insides_vs_outsides_offsets():
    print('\nPanels A & B')
    fiber_density.simulations.inner_correlation_vs_outer_correlation_offsets.main()


def figure_3_all():
    print('\nFigure 3 - Main Figure')
    figure_3_main()

    print('\nFigure SI 3')
    print('\nFigure SI 3 - Detrending Simulations (S5)')
    figure_si_3_detrending_simulations()

    print('\nFigure SI 3 - Insides vs. Outsides Distances (S6)')
    figure_si_3_insides_vs_outsides_distances()

    print('\nFigure SI 3 - Insides vs. Outsides Offsets (S7)')
    figure_si_3_insides_vs_outsides_offsets()


def figure_4_main():
    print('\nPanel B')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation.main(_low_connectivity=False)

    print('\nPanel C')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_pair_distance.main(
        _low_connectivity=False)

    print('\nPanel D')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_band=True,
                                                                                         _high_time_resolution=False)

    print('\nPanel E')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_vs_pair_distance.main(
        _band=True, _high_time_resolution=False)

    print('\nPanel G')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_triplet.main()

    print('\nPanel H')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_triplets.main()


def figure_si_4_same_vs_diff_offsets():
    print('\nPanel A')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_offset_y=0, _offset_z=0)

    print('\nPanel B')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_offset_y=-0.5, _offset_z=0)

    print('\nPanel C')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_offset_y=0.5, _offset_z=0)

    print('\nPanel D')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_offset_y=0, _offset_z=-0.5)

    print('\nPanel E')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_offset_y=0, _offset_z=0.5)


def figure_si_4_detrending_experiments():
    print('\nPanel A')
    fiber_density.experiments.inner_correlation_by_derivatives_cell_pairs.main(_directions=['inside'])

    print('\nPanels B & C')
    fiber_density.experiments.stationary_vs_inner_dynamics_derivatives_cell_pairs.main()

    print('\nPanel D')
    fiber_density.experiments.inner_correlation_by_derivatives_single_cells.main()

    print('\nPanels E & F')
    fiber_density.experiments.stationary_vs_inner_dynamics_derivatives_single_cells.main()


def figure_si_4_same_vs_diff_distances():
    print('\nPanel A')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _band=True, _high_time_resolution=False, _pair_distance_range=[4, 6])

    print('\nPanel B')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _band=True, _high_time_resolution=False, _pair_distance_range=[6, 8])

    print('\nPanel C')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _band=True, _high_time_resolution=False, _pair_distance_range=[8, 10])


def figure_si_4_same_vs_diff_controls():
    print('\nPanel B')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _real_cells=False, _static=False, _band=True, _high_time_resolution=False)

    print('\nPanel D')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _real_cells=False, _static=True, _band=False, _high_time_resolution=False)


def figure_si_4_real_vs_fake():
    print('\nPanels A & B')
    fiber_density.experiments.same_real_inner_correlation_vs_different_fake_inner_correlation.main(
        _offset_y=0, _high_time_resolution=False)

    print('\nPanels C & D')
    fiber_density.experiments.same_real_inner_correlation_vs_different_fake_inner_correlation.main(
        _offset_y=0.5, _high_time_resolution=False)

    print('\nPanels E & F')
    fiber_density.experiments.same_real_inner_correlation_vs_different_fake_inner_correlation.main(
        _offset_y=0, _high_time_resolution=True)

    print('\nPanels G & H')
    fiber_density.experiments.same_real_inner_correlation_vs_different_fake_inner_correlation.main(
        _offset_y=0.5, _high_time_resolution=True)


def figure_si_4_same_vs_diff_high_time_res():
    print('\nPanel A')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _band=True, _high_time_resolution=True)

    print('\nPanel B')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _band=False, _high_time_resolution=True)


def figure_si_4_temporal_assessment():
    print('\nPanel A')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_temporal_assessment.main(
        _high_time_resolution=False)

    print('\nPanel B')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_temporal_assessment.main(
        _high_time_resolution=True)


def figure_4_all():
    print('\nFigure 4 - Main Figure')
    figure_4_main()

    print('\nFigure SI 4')
    print('\nFigure SI 4 - Detrending Experiments (S8)')
    figure_si_4_detrending_experiments()

    print('\nFigure SI 4 - Same vs. Different - Offsets (S9)')
    figure_si_4_same_vs_diff_offsets()

    print('\nFigure SI 4 - Same vs. Different - Distances (S10)')
    figure_si_4_same_vs_diff_distances()

    print('\nFigure SI 4 - Same vs. Different - High Temporal Resolution (S11)')
    figure_si_4_same_vs_diff_high_time_res()

    print('\nFigure SI 4 - Temporal Assessment (S12)')
    figure_si_4_temporal_assessment()

    print('\nFigure SI 4 - Same vs. Different - Controls (S13)')
    figure_si_4_same_vs_diff_controls()

    print('\nFigure SI 4 - Real vs. Fake (S14)')
    figure_si_4_real_vs_fake()


def figure_5_main():
    print('\nPanel A')
    fiber_density.experiments.inner_density_vs_offsets_in_axes.main(_high_time_resolution=False)

    print('\nPanel B')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_offsets_in_axes.main(
        _high_time_resolution=False)

    print('\nPanel C')
    fiber_density.experiments.inner_density_vs_same_inner_correlation_vs_different_inner_correlation_vs_z_offset.main()

    print('\nPanel F')
    fiber_density.experiments.inner_density_vs_time.compute_tuples([('SN16', 3, 'cells_0_1')])

    print('\nPanel G')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_band=False,
                                                                                         _high_time_resolution=False)

    print('\nPanel H')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_band_vs_no_band.main(
        _high_time_resolution=False)


def figure_si_5_saturation():
    print('\nPanel A')
    fiber_density.experiments.inner_correlation_vs_saturation.main(_offset_y=0)

    print('\nPanel B')
    fiber_density.experiments.inner_correlation_vs_saturation.main(_offset_y=0.5)


def figure_si_5_high_time_res():
    print('\nPanel A')
    fiber_density.experiments.inner_density_vs_offsets_in_axes.main(_high_time_resolution=True)

    print('\nPanel B')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_offsets_in_axes.main(
        _high_time_resolution=True)


def figure_si_5_simulations():
    print('\nPanel A')
    fiber_density.simulations.inner_density_vs_offsets_in_axes.main(_low_connectivity=False)

    print('\nPanel B')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_offsets_in_axes.main(
        _low_connectivity=False)


def figure_si_5_no_band():
    print('\nPanel A')
    fiber_density.experiments.inner_density_vs_offsets_in_axes.main(_band=False, _high_time_resolution=False)

    print('\nPanel B')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_offsets_in_axes.main(
        _band=False, _high_time_resolution=False)


def figure_5_all():
    print('\nFigure 5 - Main Figure')
    figure_5_main()

    print('\nFigure SI 5')
    print('\nFigure SI 5 - Saturation (S16)')
    figure_si_5_saturation()

    print('\nFigure SI 5 - High Temporal Resolution (S17)')
    figure_si_5_high_time_res()

    # print('\nFigure SI 5 - Simulations')
    # figure_si_5_simulations()

    print('\nFigure SI 5 - No Band (S18)')
    figure_si_5_no_band()


def figure_6_main():
    print('\nPanels B & C')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(
        _real_cells=True, _offset_y=0.5, _high_time_resolution=False)

    print('\nPanel D')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(
        _real_cells=True, _offset_y=0, _high_time_resolution=False)

    print('\nPanel E')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(
        _real_cells=False, _offset_y=0.5, _high_time_resolution=False)

    print('\nPanel F')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(
        _real_cells=True, _offset_y=0.5, _high_time_resolution=True)

    print('\nPanel G')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(
        _real_cells=True, _offset_y=0.5, _high_time_resolution=True, _band=False)


def figure_6_all():
    print('\nFigure 6 - Main Figure')
    figure_6_main()


def figure_7_main():
    print('\nPanel B')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation.main(
        _alpha=1, _beta=1, _low_connectivity=False, _plots=['same'], _plot_types=['box'])

    print('\nPanel C')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _type='alpha', _low_connectivity=False, _plots=['same'], _plot_types=['box'])

    print('\nPanel D')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _type='alpha', _low_connectivity=False, _plots=['same'], _plot_types=['stacked_bar']
    )

    print('\nPanel E')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _type='beta', _low_connectivity=False, _plots=['same'], _plot_types=['box'])

    print('\nPanel F')
    fiber_density.experiments.inner_dynamics_granger_causality.main(
        _band=None, _high_time_resolution=True, _tuples_to_mark=[('SN45', 2, 'cells_2_3'), ('SN45', 1, 'cells_1_3')],
        _plots=['granger'])


def figure_si_7_simulate_time_lag_distribution():
    print('\nPanel A')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation.main(
        _alpha=1, _beta=1, _low_connectivity=False, _plots=['same'], _plot_types=['bar'])

    print('\nPanel B')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _type='alpha', _low_connectivity=False, _plots=['same'], _plot_types=['bar'])

    print('\nPanel C')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _type='beta', _low_connectivity=False, _plots=['same'], _plot_types=['bar'])


def figure_si_7_leader_follower_cross_correlation_experiments():
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_cross_correlation.main(
        _plots=['same', 'different'], _plot_types=['box', 'bar'])


def figure_si_7_different_network():
    print('\nPanels A & B')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation.main(
        _alpha=1, _beta=1, _low_connectivity=False, _plots=['different'])

    print('\nPanels C & D')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _type='alpha', _low_connectivity=False, _plots=['different'])

    print('\nPanels E & F')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _type='beta', _low_connectivity=False, _plots=['different'])


def figure_si_7_granger_causality():
    print('\nPanels B & C')
    fiber_density.experiments.inner_dynamics_granger_causality.main(
        _band=None, _high_time_resolution=True, _tuples_to_plot=[('SN45', 2, 'cells_2_3')])

    print('\nPanels E & F')
    fiber_density.experiments.inner_dynamics_granger_causality.main(
        _band=None, _high_time_resolution=True, _tuples_to_plot=[('SN45', 1, 'cells_1_3')])


def figure_si_7_whiteness():
    fiber_density.experiments.inner_dynamics_granger_causality.main(
        _band=None, _high_time_resolution=True, _tuples_to_mark=[('SN45', 2, 'cells_2_3'), ('SN45', 1, 'cells_1_3')],
        _plots=['whiteness'])


def figure_7_all():
    print('\nFigure 7 - Main Figure')
    figure_7_main()

    print('\nFigure SI 7')
    print('\nFigure SI 7 - Simulate Time Lag Distribution (S19)')
    figure_si_7_simulate_time_lag_distribution()

    print('\nFigure SI 7 - Different Network (S20)')
    figure_si_7_different_network()

    print('\nFigure SI 7 - Leader Follower Cross Correlation Experiments (S21)')
    figure_si_7_leader_follower_cross_correlation_experiments()

    print('\nFigure SI 7 - Granger Causality (S22)')
    figure_si_7_granger_causality()

    print('\nFigure 7 SI - Whiteness (S23)')
    figure_si_7_whiteness()


def figure_si_discussion_plasticity():
    fiber_density.experiments.inner_density_vs_same_inner_correlation_vs_different_inner_correlation_offsets_in_axes.\
        main()


def figure_discussion_all():
    print('\nFigure SI Discussion - Plasticity (S24)')
    figure_si_discussion_plasticity()


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

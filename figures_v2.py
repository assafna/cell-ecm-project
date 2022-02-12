import fiber_density


def main():
    # figure 1
    print('\nPanel C')
    fiber_density.simulations.inner_density_vs_window_distance_with_single_cells.main(_low_connectivity=False)

    # figure 2
    print('\nPanel A')
    fiber_density.simulations.inner_density_vs_time_cell_pairs.main()

    print('\nPanel C')
    fiber_density.simulations.communicating_inner_correlation_vs_non_communicating_inner_correlation.main()

    # figure s1
    print('\nPanel A')
    fiber_density.simulations.inner_correlation_communicating_vs_non_communicating_vs_derivatives.main()

    print('\nPanel B')
    fiber_density.simulations.inner_density_vs_change_in_inner_density.main()

    print('\nPanel C')
    fiber_density.simulations.inner_correlation_communicating_vs_non_communicating_vs_pair_distance.main()

    print('\nPanel D')
    fiber_density.simulations.inner_correlation_communicating_vs_non_communicating_vs_window_distance.main()

    # figure 3
    print('\nPanel C')
    fiber_density.experiments.inner_density_vs_window_distance_with_single_cells.main()

    print('\nPanel E')
    fiber_density.experiments.inner_density_vs_time_with_single_cells.main()

    # figure s3
    print('\nPanels B (main panel) & C')
    fiber_density.experiments.inner_density_vs_pair_distance_correlations_vs_window_distance.main()

    print('\nPanel B (secondary panels)')
    fiber_density.experiments.inner_density_vs_pair_distance.main()

    # figure 4
    print('\nPanel B - left')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation.main()

    print('\nPanel B - right')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_pair_distance.main()

    print('\nPanel D - left')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _band=True, _high_temporal_resolution=False)

    print('\nPanel D - right')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_vs_pair_distance.main(
        _band=True, _high_temporal_resolution=False)

    print('\nPanel E - middle')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_triplet.main()

    print('\nPanel E - right')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_triplets.main()

    # figure s5
    print('\nPanel A')
    fiber_density.experiments.inner_correlation_by_derivatives_cell_pairs.main(_directions=['inside'])

    print('\nPanels B & C')
    fiber_density.experiments.stationary_vs_inner_dynamics_derivatives_cell_pairs.main()

    print('\nPanel D')
    fiber_density.experiments.inner_correlation_by_derivatives_single_cells.main()

    print('\nPanels E & F')
    fiber_density.experiments.stationary_vs_inner_dynamics_derivatives_single_cells.main()

    # figure s6
    print('\nTop panel')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_offset_y=0.5, _offset_z=0)

    print('\nLeft panel')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_offset_y=0, _offset_z=-0.5)

    print('\nMiddle panel')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_offset_y=0, _offset_z=0)

    print('\nRight panel')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_offset_y=0, _offset_z=0.5)

    print('\nBottom panel')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_offset_y=-0.5, _offset_z=0)

    # figure s7
    print('\nPanel A')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _band=True, _high_temporal_resolution=False, _pair_distance_range=[4, 6])

    print('\nPanel B')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _band=True, _high_temporal_resolution=False, _pair_distance_range=[6, 8])

    print('\nPanel C')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _band=True, _high_temporal_resolution=False, _pair_distance_range=[8, 10])

    # figure 5
    print('\nPanel B')
    fiber_density.experiments.real_real_inner_correlation_vs_real_fake_inner_correlation.main()

    print('\nPanels D & E')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(
        _real_cells=True, _offset_y=0.5, _high_temporal_resolution=False)

    print('\nPanel F')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(
        _real_cells=True, _offset_y=0, _high_temporal_resolution=False)

    print('\nPanel G')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(
        _real_cells=False, _offset_y=0.5, _high_temporal_resolution=False)

    # figure s8
    print('\nPanel B')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _real_cells=False, _static=False, _band=True, _high_temporal_resolution=False)

    print('\nPanels C')
    fiber_density.experiments.same_real_inner_correlation_vs_different_fake_inner_correlation.main(
        _offset_y=0, _high_temporal_resolution=False)

    print('\nPanels D')
    fiber_density.experiments.same_real_inner_correlation_vs_different_fake_inner_correlation.main(
        _offset_y=0.5, _high_temporal_resolution=False)

    print('\nPanels E')
    fiber_density.experiments.same_real_inner_correlation_vs_different_fake_inner_correlation.main(
        _offset_y=0, _high_temporal_resolution=True)

    print('\nPanels F')
    fiber_density.experiments.same_real_inner_correlation_vs_different_fake_inner_correlation.main(
        _offset_y=0.5, _high_temporal_resolution=True)

    # figure s9
    # TODO: add new figure

    # figure s10
    print('\nPanel A')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_high_temporal_resolution=True)

    print('\nPanel B - left')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_temporal_assessment.main(
        _high_temporal_resolution=False)

    print('\nPanel B - right')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_temporal_assessment.main()

    print('\nPanel C')
    fiber_density.experiments.real_real_inner_correlation_vs_real_fake_inner_correlation.main(
        _high_temporal_resolution=True)

    print('\nPanel D')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(_high_temporal_resolution=True)

    # figure 6
    print('\nPanel A - left')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main()

    print('\nPanel A - middle')
    fiber_density.experiments.real_real_inner_correlation_vs_real_fake_inner_correlation.main()

    print('\nPanel A - right')
    fiber_density.experiments.matchmaking_by_inner_correlation.main()

    print('\nPanel B - left')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_bead=True, _band=False)

    print('\nPanel B - middle')
    fiber_density.experiments.real_real_inner_correlation_vs_real_fake_inner_correlation.main(_bead=True, _band=False)

    print('\nPanel B - right')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(_bead=True, _band=False)

    print('\nPanel C - left')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _dead_dead=True, _dead=True, _band=False)

    print('\nPanel C - middle')
    fiber_density.experiments.real_real_inner_correlation_vs_real_fake_inner_correlation.main(
        _dead_dead=True, _dead=True, _band=False)

    print('\nPanel C - right')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(
        _dead_dead=True, _dead=True, _band=False)

    print('\nPanel D - left')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _live_dead=True, _dead=True, _band=False)

    print('\nPanel D - middle')
    fiber_density.experiments.real_real_inner_correlation_vs_real_fake_inner_correlation.main(
        _live_dead=True, _dead=True, _band=False)

    print('\nPanel D - right')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(_live_dead=True, _dead=True, _band=False)

    print('\nPanel E - left')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _live_dead=True, _dead=True, _live=True, _band=False)

    print('\nPanel E - middle')
    fiber_density.experiments.real_real_inner_correlation_vs_real_fake_inner_correlation.main(
        _live_dead=True, _dead=True, _live=True, _band=False)

    print('\nPanel E - right')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(
        _live_dead=True, _dead=True, _live=True, _band=False)

    # figure 7
    print('\nPanel B')
    fiber_density.experiments.inner_density_vs_offsets_in_axes.main()

    print('\nPanel C')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_offsets_in_axes.main()

    print('\nPanel D')
    fiber_density.experiments.inner_density_vs_same_inner_correlation_vs_different_inner_correlation_vs_z_offset.main()

    print('\nPanel G')
    fiber_density.experiments.inner_density_vs_time_tuple.compute_tuples([('SN16', 3, 'cells_0_1')])

    print('\nPanel H - left')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(_band=False)

    print('\nPanel H - middle')
    fiber_density.experiments.real_real_inner_correlation_vs_real_fake_inner_correlation.main(_band=False)

    print('\nPanel H - right')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(_band=False)

    print('\nPanel I - left')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _band=False, _high_temporal_resolution=True)

    print('\nPanel I - middle')
    fiber_density.experiments.real_real_inner_correlation_vs_real_fake_inner_correlation.main(
        _band=False, _high_temporal_resolution=True)

    print('\nPanel I - right')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(
        _band=False, _high_temporal_resolution=True)

    # figure s11
    print('\nPanel A')
    fiber_density.experiments.inner_correlation_vs_saturation.main(_offset_y=0)

    print('\nPanel B')
    fiber_density.experiments.inner_correlation_vs_saturation.main()

    # figure s12
    print('\nPanel A')
    fiber_density.experiments.inner_density_vs_offsets_in_axes.main(_high_temporal_resolution=True)

    print('\nPanel B')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_offsets_in_axes.main(
        _high_temporal_resolution=True)

    # figure s13
    print('\nPanel A')
    fiber_density.experiments.inner_density_vs_offsets_in_axes.main(_band=False)

    print('\nPanel B')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_offsets_in_axes.main(_band=False)

    # figure s15
    print('\nPanel A - left')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _band=False, _bleb=True, _bleb_amount_um=85)

    print('\nPanel A - middle')
    fiber_density.experiments.real_real_inner_correlation_vs_real_fake_inner_correlation.main(
        _band=False, _bleb=True, _bleb_amount_um=85)

    print('\nPanel A - right')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(_band=False, _bleb=True, _bleb_amount_um=85)

    print('\nPanel B - left')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _band=False, _bleb=True, _bleb_amount_um=150)

    print('\nPanel B - middle')
    fiber_density.experiments.real_real_inner_correlation_vs_real_fake_inner_correlation.main(
        _band=False, _bleb=True, _bleb_amount_um=150)

    print('\nPanel B - right')
    fiber_density.experiments.matchmaking_by_inner_correlation.main(_band=False, _bleb=True, _bleb_amount_um=150)

    # figure s14
    print('\nPanel A')
    fiber_density.experiments.inner_density_vs_same_inner_correlation_vs_different_inner_correlation_offsets_in_axes.\
        main()

    # TODO: fix red line for panel B graphs
    print('\nPanel B - top left')
    fiber_density.experiments.bleb_before_vs_after_inner_density.main()

    print('\nPanel B - top right')
    fiber_density.experiments.bleb_before_vs_after_inner_density.main(_offset_y=0.5)

    print('\nPanel B - bottom left')
    fiber_density.experiments.bleb_before_vs_after_inner_density_change.main()

    print('\nPanel B - bottom right')
    fiber_density.experiments.bleb_before_vs_after_inner_density_change.main(_offset_y=0.5)

    print('\nPanel C - top left')
    fiber_density.experiments.inner_density_vs_offsets_in_axes.main()

    print('\nPanel C - top middle')
    fiber_density.experiments.inner_density_change_vs_offsets_in_axes.main()

    print('\nPanel C - top right')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_offsets_in_axes.main()

    print('\nPanel C - bottom left')
    fiber_density.experiments.inner_density_vs_offsets_in_axes.main(_high_temporal_resolution=True)

    print('\nPanel C - bottom middle')
    fiber_density.experiments.inner_density_change_vs_offsets_in_axes.main(_high_temporal_resolution=True)

    print('\nPanel C - bottom right')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_offsets_in_axes.main(
        _high_temporal_resolution=True)

    # figure s16
    print('\nPanel B')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation.main(
        _plots=['same'], _plot_types=['box'])

    print('\nPanel C')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _plots=['same'], _plot_types=['box'])

    print('\nPanel D')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _plots=['same'], _plot_types=['stacked_bar']
    )

    print('\nPanel E')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _type='beta', _plots=['same'], _plot_types=['box'])

    # figure s17
    print('\nPanel A')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation.main(
        _plots=['same'], _plot_types=['bar'])

    print('\nPanel B')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _plots=['same'], _plot_types=['bar'])

    print('\nPanel C')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _type='beta', _plots=['same'], _plot_types=['bar'])

    # figure s18
    print('\nPanels A & B')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation.main(
        _plots=['different'])

    print('\nPanels C & D')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _plots=['different'])

    print('\nPanels E & F')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_cross_correlation_alphas_betas.main(
        _type='beta', _plots=['different'])

    # figure s19
    print('\nAll panels')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_cross_correlation.main(
        _plots=['same', 'different'], _plot_types=['box', 'bar'])


if __name__ == '__main__':
    main()

import fiber_density


def main():
    # figure 1
    print('\nPanel C')
    fiber_density.simulations.inner_density_vs_window_distance_with_single_cells.main(_low_connectivity=False)

    # figure 2
    print('\nPanel A')
    fiber_density.simulations.inner_density_vs_time_cell_pairs.main()

    print('\nPanel F')
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
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation.main(_low_connectivity=False)

    print('\nPanel B - right')
    fiber_density.simulations.same_inner_correlation_vs_different_inner_correlation_pair_distance.main(
        _low_connectivity=False)

    print('\nPanel D - top')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation.main(
        _band=True, _high_temporal_resolution=False)

    print('\nPanel D - bottom')
    fiber_density.experiments.same_inner_correlation_vs_different_inner_correlation_vs_pair_distance.main(
        _band=True, _high_temporal_resolution=False)

    print('\nPanel E - left')
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


if __name__ == '__main__':
    main()

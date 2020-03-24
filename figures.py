import fibers_density


def figure_1():
    fibers_density.fibers_vs_distance_single_cells.main()
    fibers_density.fibers_vs_distance_pairs.main()
    fibers_density.fibers_vs_distance_differences.main()
    fibers_density.simulations.fibers_vs_distance_multiple_cell_distances.main()
    fibers_density.experiments.fibers_vs_distance_multiple_cell_distances.main()


def figure_2():
    fibers_density.fibers_in_time_single_cells.main()
    fibers_density.fibers_in_time_pairs.main()


def figure_3():
    fibers_density.simulations.fibers_vs_derivatives_single_cell.main()
    fibers_density.simulations.insides_vs_outsides_stds.main()
    # TODO: run when simulations are available
    fibers_density.simulations.insides_vs_outsides_cell_distances.main()


def figure_4():
    # TODO: create a single graph for all derivatives
    fibers_density.simulations.communicated_vs_non_communicated.main()
    fibers_density.experiments.communicated_vs_non_communicated.main()
    fibers_density.simulations.communicated_vs_non_communicated_cell_distances.main()


if __name__ == '__main__':
    figure_1()
    figure_2()
    figure_3()
    figure_4()

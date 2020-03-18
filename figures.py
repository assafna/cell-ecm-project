import fibers_density.simulations as simulations
import fibers_density.experiments as experiments
import fibers_density


def figure_1():
    fibers_density.fibers_vs_distance_single_cells.main()
    fibers_density.fibers_vs_distance_pairs.main()
    fibers_density.fibers_vs_distance_differences.main()
    simulations.fibers_vs_distance_multiple_cell_distances.main()
    experiments.fibers_vs_distance_multiple_cell_distances.main()

    # sup
    fibers_density.fibers_vs_distance_single_cells.main(_low_connectivity=True)
    fibers_density.fibers_vs_distance_pairs.main(_low_connectivity=True)
    fibers_density.fibers_vs_distance_differences.main(_low_connectivity=True)
    simulations.fibers_vs_distance_multiple_cell_distances.main(_low_connectivity=True)


def figure_2():
    fibers_density.fibers_in_time_single_cells.main()
    fibers_density.fibers_in_time_pairs.main()


def figure_3():
    simulations.fibers_vs_derivatives_single_cell.main()
    simulations.insides_vs_outsides_stds.main()
    # TODO: run when simulations are available
    simulations.insides_vs_outsides_cell_distances.main()


def figure_4():
    # TODO: create a single graph for all derivatives
    simulations.master_vs_slave.main()
    experiments.master_vs_slave.main()
    simulations.master_vs_slave_cell_distances.main()


if __name__ == '__main__':
    figure_1()
    figure_2()
    figure_3()
    figure_4()

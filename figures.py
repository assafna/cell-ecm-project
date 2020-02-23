import fibers_density.simulations as simulations
import fibers_density.experiments as experiments


def figure_1():
    simulations.fibers_vs_distance.main()
    experiments.fibers_vs_distance.main()


def figure_2():
    simulations.fibers_in_time.main()


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

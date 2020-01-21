from fibers_density.simulations import fibers_vs_distance, fibers_in_time, different_single_cell_correlations, \
    insides_vs_outsides_stds, insides_vs_outsides_cell_distances, master_vs_slave, master_vs_slave_cell_distances


def figure_1():
    fibers_vs_distance.main()


def figure_2():
    fibers_in_time.main()


def figure_3():
    different_single_cell_correlations.main()
    insides_vs_outsides_stds.main()
    # TODO: run when simulations are available
    insides_vs_outsides_cell_distances.main()


def figure_4():
    master_vs_slave.main()
    master_vs_slave_cell_distances.main()


if __name__ == '__main__':
    figure_1()
    figure_2()
    figure_3()
    figure_4()

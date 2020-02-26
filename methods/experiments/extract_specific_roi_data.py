from libs import compute_lib
from libs.experiments import load, compute
from methods.experiments import export_video


def main():
    _experiment = 'SN16'
    _series_id = 1
    _cells = '2_3'
    _group_real = 'cells_' + _cells
    _group_fake = 'fake_' + _cells
    _roi_length = 1
    _roi_height = 1
    _roi_width = 1
    _offset_x = 0
    _offset_y = 0
    _offset_z = 0
    _direction = 'inside'
    _start_time_point = 0
    _end_time_point = 500
    _derivative = 0
    _arguments = {
        'experiment': _experiment,
        'series_id': _series_id,
        'group': _group_real,
        'length_x': _roi_length,
        'length_y': _roi_height,
        'length_z': _roi_width,
        'offset_x': _offset_x,
        'offset_y': _offset_y,
        'offset_z': _offset_z,
        'direction': _direction,
        'time_points': _end_time_point
    }
    _fibers_densities_real = compute.roi_fibers_density_by_time_pairs(_arguments)
    # _fibers_densities_fake, _out_of_boundaries = compute.roi_fibers_density_by_time_pairs(
    #     _experiment, _series_id, _group_fake, _roi_length, _roi_height, _roi_width,
    #     _offset_x, _offset_y, _offset_z_, _direction, _end_time_point
    # )
    # for _end in range(_start_time_point + 5, _end_time_point):
    #     _correlation = compute_lib.correlation(
    #         compute_lib.derivative(_fibers_densities_real['left_cell'][_start_time_point:_end], _n=_derivative),
    #         compute_lib.derivative(_fibers_densities_real['right_cell'][_start_time_point:_end], _n=_derivative)
    #     )
    #     print(_correlation)
    # export_video.process_group(_experiment, _series_id, _group_fake)
    for _left, _right in zip(_fibers_densities_real['left_cell'], _fibers_densities_real['right_cell']):
        print(_left[0], _right[0], sep='\t')


if __name__ == '__main__':
    main()

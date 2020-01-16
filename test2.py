import os

from libs import compute_lib

if __name__ == '__main__':
    _path = 'G:\My Drive\BGU\Thesis\Cell-ECM & Cell-ECM-Cell Project\Old\Data\Experiments\Manipulations\Fibers Density\SN16_CZI\Series 1\cells_3_4\cells_1_2'
    _list = os.listdir(_path)
    _left_array = list([None] * 15)
    _right_array = list([None] * 15)
    for _file in _list:
        _time_point = int(str(_file.split('_')[1]).split('.txt')[0])
        if _file.endswith('.txt'):
            with open(os.path.join(_path, _file), 'r') as _f:
                _lines = _f.readlines()
                _left_array[_time_point - 1] = float(str(_lines[5].split('\t')[1]).split('\n')[0])
                _right_array[_time_point - 1] = float(str(_lines[-5].split('\t')[1]).split('\n')[0])

    _cor = compute_lib.correlation(_left_array, _right_array)
    _cor1 = compute_lib.correlation(compute_lib.derivative(_left_array, 1), compute_lib.derivative(_right_array, 1))
    _cor2 = compute_lib.correlation(compute_lib.derivative(_left_array, 2), compute_lib.derivative(_right_array, 2))

    print(_cor, _cor1, _cor2)

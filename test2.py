import os

from libs.load_lib import from_pickle
from libs.save_lib import to_pickle
from libs.simulations import paths


def main():
    _raw = os.listdir(paths.RAW)
    for _simulation in _raw:

        if '_a_' in _simulation:
            _part_1, _part_2 = _simulation.split('_a_')

            if _part_2 == '00':
                _new_name = _part_1 + '_alpha_000'
            elif _part_2 == '025':
                _new_name = _part_1 + '_alpha_025_beta_100'
            elif _part_2 == '050':
                _new_name = _part_1 + '_alpha_050_beta_100'
            elif _part_2 == '075':
                _new_name = _part_1 + '_alpha_075_beta_100'
            else:
                _new_name = _simulation

            _old_path = paths.raw(_simulation)
            _new_path = paths.raw(_new_name)
            os.rename(_old_path, _new_path)

        elif '_beta_' in _simulation:
            _part_1, _part_2 = _simulation.split('_beta_')

            if _part_2 == '105':
                _new_name = _part_1 + '_alpha_100_beta_105'
            elif _part_2 == '110':
                _new_name = _part_1 + '_alpha_100_beta_110'
            elif _part_2 == '120':
                _new_name = _part_1 + '_alpha_100_beta_120'
            else:
                _new_name = _simulation

            _old_path = paths.raw(_simulation)
            _new_path = paths.raw(_new_name)
            os.rename(_old_path, _new_path)


if __name__ == '__main__':
    main()

import gzip
import os
import pickle


def to_pickle(_object, _pickle_path):
    os.makedirs(os.path.dirname(_pickle_path)) if not os.path.isdir(os.path.dirname(_pickle_path)) else None
    try:
        with gzip.open(_pickle_path, 'wb') as _pickle:
            pickle.dump(_object, _pickle, -1)
    finally:
        _pickle.close()

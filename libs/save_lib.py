import gzip
import json
import os
import pickle


def to_pickle(_object, _pickle_path):
    os.makedirs(os.path.dirname(_pickle_path)) if not os.path.isdir(os.path.dirname(_pickle_path)) else None
    with gzip.open(_pickle_path, 'wb') as _pickle:
        pickle.dump(_object, _pickle, -1)


def to_json(_json_data, _json_path):
    os.makedirs(os.path.dirname(_json_path)) if not os.path.isdir(os.path.dirname(_json_path)) else None
    with open(_json_path, 'w') as _json:
        json.dump(_json_data, _json, indent=4)

import gzip
import json
import pickle


def from_pickle(_path):
    try:
        with gzip.open(_path, 'rb') as _pickle:
            return pickle.load(_pickle)
    except (FileNotFoundError, EOFError):
        return {}


def from_json(_path):
    try:
        with open(_path) as _json:
            return json.load(_json)
    except (FileNotFoundError, EOFError):
        return {}

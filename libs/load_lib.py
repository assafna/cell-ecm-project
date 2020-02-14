import gzip
import json
import pickle


def from_pickle(_path):
    try:
        with gzip.open(_path, 'rb') as _pickle:
            return pickle.load(_pickle)
    finally:
        _pickle.close()


def from_json(_path):
    try:
        with open(_path) as _json:
            return json.load(_json)
    finally:
        _json.close()

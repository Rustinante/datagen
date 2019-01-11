import numpy as np


def map_counts_to_vec(a, g, c, t, x):
    # A G C T X
    return np.array([a, g, c, t, x])


def map_counts_to_revcomp_vec(a, g, c, t, x):
    # A G C T X
    return np.array([t, c, g, a, x])


mapping = {
    'a': np.array([1, 0, 0, 0, 0]),
    'A': np.array([1, 0, 0, 0, 0]),

    'g': np.array([0, 1, 0, 0, 0]),
    'G': np.array([0, 1, 0, 0, 0]),

    'c': np.array([0, 0, 1, 0, 0]),
    'C': np.array([0, 0, 1, 0, 0]),

    't': np.array([0, 0, 0, 1, 0]),
    'T': np.array([0, 0, 0, 1, 0]),

    'X': np.array([0, 0, 0, 0, 1]),

    'N': np.array([0, 0, 0, 0, 0]),
    'n': np.array([0, 0, 0, 0, 0])
}

complement_mapping = {
    'a': mapping['t'],
    'A': mapping['T'],

    'g': mapping['c'],
    'G': mapping['C'],

    'c': mapping['g'],
    'C': mapping['G'],

    't': mapping['a'],
    'T': mapping['A'],

    'X': mapping['X'],

    'N': mapping['N'],
    'n': mapping['n']
}

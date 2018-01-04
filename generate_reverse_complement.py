import pickle
import numpy as np

def generate_reverse_complement(filename, output_filename):
    print('output filename: {}'.format(output_filename))
    complement_map = {
        'a': 't',
        't': 'a',
        'c': 'g',
        'g': 'c',
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'N': 'N',
        'X': 'X'
    }
    reverse_complement_list = []
    with open(filename, 'rb') as file:
        data = pickle.load(file)

    for matrix in data:
        reverse_complement_matrix = np.empty_like(matrix)
        # each matrix is 100 x 1000 where each row corresponds to one species
        for row in range(100):
            strand = matrix[row]
            reverse_complement_strand = reverse_complement_matrix[row]
            for column in range(1000):
                reverse_complement_strand[999 - column] = complement_map[strand[column]]

        reverse_complement_list.append(reverse_complement_matrix)

    print('=> pickling...')
    with open(output_filename, 'wb') as output_file:
        pickle.dump(reverse_complement_list, output_file, protocol=pickle.HIGHEST_PROTOCOL)
    print('=> finished pickling')

import argparse
import h5py
import numpy as np
import os
import time

from binary_search import search, get_location_from_line

mapping = {
    'a': np.array([1, 0, 0, 0]),
    'A': np.array([1, 0, 0, 0]),

    'g': np.array([0, 1, 0, 0]),
    'G': np.array([0, 1, 0, 0]),

    'c': np.array([0, 0, 1, 0]),
    'C': np.array([0, 0, 1, 0]),

    't': np.array([0, 0, 0, 1]),
    'T': np.array([0, 0, 0, 1]),

    'X': np.array([0, 0, 0, 0]),

    'N': np.array([0, 0, 0, 0])
}


def scan_through_line_for_number(alignment_file, start_line_hint, number):
    alignment_file.seek(start_line_hint)

    for line in alignment_file:
        location = get_location_from_line(line)

        if not location:
            raise ValueError('There still lines in the alignment file but cannot obtain coordinate location')

        if location == number:
            return line, start_line_hint

        # The lines are sorted so if the current location is already greater than the one
        # we're searching for we know what we search does not exist.
        elif location > number:
            return None

        start_line_hint += len(bytes(line, 'ascii'))

    return None


def extend_dataset(chr, purpose):
    array_list = []

    coordinate_filename = os.path.join('data', '{}_{}'.format(chr, purpose))
    alignment_filename = '{}_maf_sequence.csv'.format(chr)
    pickle_filename = '{}_{}.align.hdf5'.format(chr, purpose)

    print('=> coordinate_filename: {}'.format(coordinate_filename))
    print('=> alignment_filename: {}'.format(alignment_filename))
    print('=> target pickle_filename: {}'.format(pickle_filename))

    with open(coordinate_filename, 'r') as file, open(alignment_filename, 'r') as alignment_file:
        header = alignment_file.readline().strip().split(',')
        human_index = header.index('hg19')

        processed_line_count = 0
        start_time = time.time()
        flanking_number = 400

        for line in file:
            processed_line_count += 1
            (start_coordinate, sequence) = line.strip().split(',')
            start_coordinate = int(start_coordinate)

            alignment_matrix = np.zeros((1000, 100, 4), dtype='uint8')

            for index, hg_letter in enumerate(sequence):
                alignment_matrix[index, 0, :] = mapping[hg_letter]

            start_line_hint = None
            for index, coordinate in enumerate(range(start_coordinate - flanking_number, start_coordinate + 200 + flanking_number)):
                if not start_line_hint:
                    result = search(alignment_file, coordinate, alignment_filename)
                else:
                    result = scan_through_line_for_number(alignment_file=alignment_file, start_line_hint=start_line_hint, number=coordinate)

                if result:
                    start_line_hint = result[1]
                    tokens = result[0].strip().split(',')
                    # Important: pop the human_index first before removing the start index.
                    # So the human_index will be the correct index.
                    tokens.pop(human_index)
                    # aligned_letters = tokens[1:]
                    # aligned_letters is 100x4
                    aligned_letters = alignment_matrix[index]
                    for species_index, letter in enumerate(tokens[1:]):
                        # +1 because we want to put hg19 in the first row
                        aligned_letters[species_index + 1, :] = mapping[letter]

            array_list.append(alignment_matrix.transpose((1, 0, 2)))

            if processed_line_count % 100 == 0:
                elapsed_time = time.time() - start_time
                time_per_line = elapsed_time / processed_line_count
                print('Processed {} lines in {:5f}s, averaging: {:5f}s per line'
                      .format(processed_line_count, elapsed_time, time_per_line))

    stamp = time.time()
    print('=> Serializing...')
    with h5py.File(pickle_filename, 'w') as file:
        feature_group = file.create_group('feature')

        feature_data = feature_group.create_dataset('data', (len(array_list), 100, 1000, 4), dtype='uint8')
        for index, matrix in enumerate(array_list):
            feature_data[index] = matrix
            # for row in range(100):
            #     species_sequence = matrix[row]
            #     for column in range(1000):
            #         feature_data[index][row][column] = mapping[species_sequence[column]]()

    print('=> Finished serializeing in {:.5f}s'.format(time.time() - stamp))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('chr')
    parser.add_argument('purpose')
    args = parser.parse_args()
    extend_dataset(args.chr, args.purpose)

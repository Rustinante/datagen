import argparse
import h5py
import numpy as np
import os
import time

from binary_search import search, get_location_from_line


class LineCache:
    def __init__(self):
        self.capacity = 1000
        self.count_per_eviction = 200
        self.cache = {}
        self.key_list = []
        
    def evict(self):
        for key in self.key_list[:self.count_per_eviction]:
            del self.cache[key]
        self.key_list = self.key_list[self.count_per_eviction:]
    
    def __contains__(self, item):
        return item in self.cache
    
    def __getitem__(self, item):
        return self.cache[item]
    
    def __setitem__(self, key, value):
        if len(self.cache) >= self.capacity:
            self.evict()
        self.cache[key] = value
        if key not in self.key_list:
            self.key_list.append(key)
        
    def __len__(self):
        return len(self.cache)
    
    def keys(self):
        return self.cache.keys()
    
    def values(self):
        return self.cache.values()
    
    def items(self):
        return self.cache.items()
        
        
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
    cache = LineCache()
    
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
            # c1=c2=c3=0
            processed_line_count += 1
            (start_coordinate, sequence) = line.strip().split(',')
            start_coordinate = int(start_coordinate)

            alignment_matrix = np.zeros((1000, 100, 4), dtype='uint8')

            for letter_index, hg_letter in enumerate(sequence):
                alignment_matrix[letter_index, 0, :] = mapping[hg_letter]

            start_line_hint = None
            for letter_index, coordinate in enumerate(range(start_coordinate - flanking_number, start_coordinate + 200 + flanking_number)):
                if coordinate in cache:
                    # c1+=1
                    cached_result = cache[coordinate]
                    alignment_matrix[letter_index][1:, :] = cached_result[0]
                    start_line_hint = cached_result[1]
                    continue
                    
                elif not start_line_hint:
                    # c2+=1
                    result = search(alignment_file, coordinate, alignment_filename)
                else:
                    # c3+=1
                    result = scan_through_line_for_number(alignment_file=alignment_file, start_line_hint=start_line_hint, number=coordinate)

                if result:
                    start_line_hint = result[1]
                    tokens = result[0].strip().split(',')
                    # Important: pop the human_index first before removing the start index,
                    # so that the human_index will be the correct index.
                    del tokens[human_index]
                    del tokens[0]
                    
                    # aligned_letters is 100x4
                    aligned_letters = alignment_matrix[letter_index]
                    for species_index, letter in enumerate(tokens):
                        # aligned_letters is 100x4
                        # +1 because we put hg19 in the first row
                        aligned_letters[species_index + 1, :] = mapping[letter]

                    cache[coordinate] = (aligned_letters[1:, :], start_line_hint)

            # print("cache: {} binary search: {} line search {}".format(c1,c2,c3))
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
        
        array_list_length = len(array_list)

        feature_data = feature_group.create_dataset('data', (array_list_length, 100, 1000, 4), dtype='uint8')
        for index, matrix in enumerate(array_list):
            feature_data[index] = matrix
            
            if index % 1000 == 0:
                print("{}/{} in {:5f}s".format(index, array_list_length, time.time() - stamp))

    print('=> Finished serializeing in {:.5f}s'.format(time.time() - stamp))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('chr')
    parser.add_argument('purpose')
    args = parser.parse_args()
    extend_dataset(args.chr, args.purpose)

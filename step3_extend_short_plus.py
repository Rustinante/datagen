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


a_array = np.array([1, 0, 0, 0, 0], dtype='uint8')
g_array = np.array([0, 1, 0, 0, 0], dtype='uint8')
c_array = np.array([0, 0, 1, 0, 0], dtype='uint8')
t_array = np.array([0, 0, 0, 1, 0], dtype='uint8')
x_array = np.array([0, 0, 0, 0, 1], dtype='uint8')
zero_array = np.array([0, 0, 0, 0, 0], dtype='uint8')

mapping = {
    'a': a_array,
    'A': a_array,
    
    'g': g_array,
    'G': g_array,
    
    'c': c_array,
    'C': c_array,
    
    't': t_array,
    'T': t_array,
    
    'X': x_array,
    
    'N': zero_array,
    'n': zero_array
}

complement_mapping = {
    'a': t_array,
    'A': t_array,
    
    'g': c_array,
    'G': c_array,
    
    'c': g_array,
    'C': g_array,
    
    't': a_array,
    'T': a_array,
    
    'X': x_array,
    
    'N': zero_array,
    'n': zero_array
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


def extend_dataset(chrom, purpose):
    checkpoint_time_str = time.strftime('%a %b %d %Y %H:%M:%S UTC%z', time.localtime(time.time()))
    print('Current time: {}'.format(checkpoint_time_str))
    
    cache = LineCache()
    
    array_list = []
    reverse_complement_array_list = []
    
    human_seq_list = []
    human_revcomp_seq_list = []
    
    coordinate_filename = os.path.join('data', '{}_{}'.format(chrom, purpose))
    alignment_filename = '{}_maf_sequence.csv'.format(chrom)
    hdf5_filename = '{}_{}.short.hdf5'.format(chrom, purpose)
    hdf5_revcomp_filename = '{}_{}.revcomp.short.hdf5'.format(chrom, purpose)
    species_indices = [42, 74, 39, 21, 78, 69, 83, 94, 81, 96, 71, 17, 75, 12]
    number_of_species = len(species_indices)
    
    print('=> coordinate_filename: {}'.format(coordinate_filename))
    print('=> alignment_filename: {}'.format(alignment_filename))
    print('=> target hdf5_filename: {}'.format(hdf5_filename))
    print('=> target reverse complement hdf5_filename: {}'.format(hdf5_revcomp_filename))
    
    with open(coordinate_filename, 'r') as file, open(alignment_filename, 'r') as alignment_file:
        header = alignment_file.readline().strip().split(',')
        
        processed_line_count = 0
        start_time = time.time()
        
        seq_len = 200
        feature_dim = 5
        human_seq_len = 1000
        human_end_index = human_seq_len - 1
        
        for line in file:
            processed_line_count += 1
            start_coordinate = int(line.strip().split(',')[0])
            
            human_seq = line.strip().split(',')[1]
            human_matrix = np.zeros((human_seq_len, feature_dim), dtype='uint8')
            human_revcomp_matrix = np.zeros((human_seq_len, feature_dim), dtype='uint8')
            for index, letter in enumerate(human_seq):
                human_matrix[index] = mapping[letter]
                human_revcomp_matrix[human_end_index - index] = complement_mapping[letter]
            human_seq_list.append(human_matrix)
            human_revcomp_seq_list.append(human_revcomp_matrix)

            # revcomp_matrix is the reverse complement of the sequences
            alignment_matrix = np.zeros((seq_len, number_of_species, feature_dim), dtype='uint8')
            revcomp_matrix = np.zeros((seq_len, number_of_species, feature_dim), dtype='uint8')

            start_line_hint = None
            end_index = seq_len - 1
            for letter_index, coordinate in enumerate(range(start_coordinate, start_coordinate + seq_len)):
                if coordinate in cache:
                    cached_result = cache[coordinate]
                    alignment_matrix[letter_index] = cached_result[0]
                    revcomp_matrix[end_index - letter_index] = cached_result[1]
                    start_line_hint = cached_result[2]

                    continue

                elif not start_line_hint:
                    result = search(alignment_file, coordinate, alignment_filename)
                else:
                    result = scan_through_line_for_number(alignment_file=alignment_file,
                                                          start_line_hint=start_line_hint, number=coordinate)

                if result:
                    start_line_hint = result[1]
                    tokens = result[0].strip().split(',')
                    # the first token is pos
                    del tokens[0]

                    # aligned_letters is of shape (num_species , feature_dim)
                    # revcomp_letters is associated at the index end_index - letter_index
                    aligned_letters = alignment_matrix[letter_index]
                    revcomp_letters = revcomp_matrix[end_index - letter_index]

                    for row_number, species_index in enumerate(species_indices):
                        letter = tokens[species_index]
                        aligned_letters[row_number, :] = mapping[letter]
                        revcomp_letters[row_number, :] = complement_mapping[letter]

                    cache[coordinate] = (aligned_letters, revcomp_letters, start_line_hint)

            array_list.append(alignment_matrix.transpose((1, 0, 2)))
            reverse_complement_array_list.append(revcomp_matrix.transpose((1, 0, 2)))

            if processed_line_count % 1000 == 0:
                elapsed_time = time.time() - start_time
                time_per_line = elapsed_time / processed_line_count
                print('Processed {} lines in {:5f}s, averaging: {:5f}s per line'
                      .format(processed_line_count, elapsed_time, time_per_line), end='\r')

    stamp = time.time()
    print('\n=> Serializing...')
    with h5py.File(hdf5_filename, 'w') as file:
        feature_group = file.create_group('feature')

        array_list_length = len(array_list)

        feature_data = feature_group.create_dataset('data',
                                                    (array_list_length, number_of_species, seq_len, feature_dim),
                                                    dtype='uint8')
        for index, matrix in enumerate(array_list):
            feature_data[index] = matrix

            if index % 1000 == 0:
                print("{}/{} in {:5f}s".format(index, array_list_length, time.time() - stamp), end='\r')

    print('\n=> Serializing the reverse complement...')
    with h5py.File(hdf5_revcomp_filename, 'w') as file:
        feature_group = file.create_group('feature')

        array_list_length = len(reverse_complement_array_list)

        feature_data = feature_group.create_dataset('data',
                                                    (array_list_length, number_of_species, seq_len, feature_dim),
                                                    dtype='uint8')
        for index, matrix in enumerate(reverse_complement_array_list):
            feature_data[index] = matrix

            if index % 1000 == 0:
                print("{}/{} in {:5f}s".format(index, array_list_length, time.time() - stamp), end='\r')

    print('\n=> Finished serializeing in {:.5f}s'.format(time.time() - stamp))
    
    
    print('\n=> Serializing human seq...')
    with h5py.File(hdf5_filename, 'r+') as file:
        feature_group = file.create_group('human_seq')
        array_list_length = len(human_seq_list)
        
        feature_data = feature_group.create_dataset('data',
                                                    (array_list_length, human_seq_len, feature_dim),
                                                    dtype='uint8')
        for index, matrix in enumerate(human_seq_list):
            feature_data[index] = matrix

            if index % 1000 == 0:
                print("{}/{} in {:5f}s".format(index, array_list_length, time.time() - stamp), end='\r')
    
    print('\n=> Finished serializing the human seq')

    print('\n=> Serializing revcomp human seq...')
    with h5py.File(hdf5_revcomp_filename, 'r+') as file:
        feature_group = file.create_group('human_seq')
        array_list_length = len(human_revcomp_seq_list)
    
        feature_data = feature_group.create_dataset('data',
                                                    (array_list_length, human_seq_len, feature_dim),
                                                    dtype='uint8')
        for index, matrix in enumerate(human_revcomp_seq_list):
            feature_data[index] = matrix
        
            if index % 1000 == 0:
                print("{}/{} in {:5f}s".format(index, array_list_length, time.time() - stamp), end='\r')

    print('\n=> Finished serializing the human seq')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('chrom')
    parser.add_argument('purpose')
    args = parser.parse_args()
    extend_dataset(args.chrom, args.purpose)

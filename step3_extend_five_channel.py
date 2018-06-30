import argparse
import h5py
import numpy as np
import os
import time

from binary_search import search, get_location_from_line


def get_line_count(filename):
    count = 0
    with open(filename, 'r') as file:
        for _ in file:
            count += 1
    return count


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
    
    matrix_list = []
    
    coordinate_filename = os.path.join('data', '{}_{}'.format(chrom, purpose))
    alignment_filename = '{}_maf_sequence.csv'.format(chrom)
    hdf5_filename = '{}_{}.hundred.hdf5'.format(chrom, purpose)
    num_species = 100
    
    print('=> coordinate_filename: {}'.format(coordinate_filename))
    print('=> alignment_filename: {}'.format(alignment_filename))
    print('=> target hdf5_filename: {}'.format(hdf5_filename))
    file_byte_size = os.stat(alignment_filename).st_size
    
    line_count = get_line_count(coordinate_filename)
    
    flanking_number = 400
    seq_len = 200 + 2 * flanking_number
    feature_dim = 5
    
    serializing_index = 0
    
    with open(coordinate_filename, 'r') as coord_file, open(alignment_filename, 'r') as alignment_file, \
            h5py.File(hdf5_filename, 'w') as hdf5_file:
        
        feature_group = hdf5_file.create_group('feature')
        feature_data = feature_group.create_dataset('data', (line_count * 2, seq_len, num_species, feature_dim), dtype='uint8')
        
        header = alignment_file.readline().strip().split(',')
        human_index = header.index('hg19')
        
        processed_line_count = 0
        start_time = time.time()
        
        for line in coord_file:
            processed_line_count += 1
            (start_coordinate, sequence) = line.strip().split(',')
            start_coordinate = int(start_coordinate)
            
            # 1000 x 100 x 5
            alignment_matrix = np.zeros((seq_len, num_species, feature_dim), dtype='uint8')
            
            start_line_hint = None
            for bp_index, (hg_letter, coordinate) in enumerate(
                    zip(sequence, range(start_coordinate - flanking_number, start_coordinate + 200 + flanking_number))):
                
                # insert the human letter first
                alignment_matrix[bp_index, 0, :] = mapping[hg_letter]
                
                if coordinate in cache:
                    cached_result = cache[coordinate]
                    alignment_matrix[bp_index, 1:, :] = cached_result[0]
                    start_line_hint = cached_result[1]
                    continue
                
                elif not start_line_hint:
                    result = search(alignment_file, coordinate, file_byte_size)
                else:
                    result = scan_through_line_for_number(alignment_file=alignment_file,
                                                          start_line_hint=start_line_hint, number=coordinate)
                
                if result:
                    start_line_hint = result[1]
                    tokens = result[0].strip().split(',')
                    # Important: pop the human_index first before removing the start index,
                    # so that the human_index will be the correct index.
                    del tokens[human_index]
                    del tokens[0]
                    
                    # aligned_letters is 100x5
                    aligned_letters = alignment_matrix[bp_index]
                    
                    for row_number, letter in enumerate(tokens):
                        # +1 because we put hg19 in the first row
                        aligned_letters[row_number + 1, :] = mapping[letter]
                    
                    cache[coordinate] = (aligned_letters[1:, :], start_line_hint)
                
                elif hg_letter != 'N' and hg_letter != 'n':
                    # broadcasting along the species dimension
                    alignment_matrix[bp_index, 1:, :] = mapping['X']
            
            matrix_list.append(alignment_matrix)
            
            if processed_line_count % 100 == 1:
                for matrix in matrix_list:
                    feature_data[serializing_index] = matrix
                    # For the reverse complement strand
                    feature_data[serializing_index + line_count] = matrix[::-1]
                    serializing_index += 1
                
                matrix_list = []
                
                elapsed_time = time.time() - start_time
                time_per_line = elapsed_time / processed_line_count
                print(f'{processed_line_count}/{line_count} = {processed_line_count/line_count:.2%} in {elapsed_time:.4f}s'
                      f'averaging {time_per_line}s per line'
                      f' current serializing index: {serializing_index}',
                      end='\r')
        
        # Serializing the remaining data
        if matrix_list:
            for matrix in matrix_list:
                feature_data[serializing_index] = matrix
                feature_data[serializing_index + line_count] = matrix[::-1]
                serializing_index += 1
            print(f'{processed_line_count}/{line_count} = {processed_line_count/line_count:.2%} in {elapsed_time:.4f}s'
                  f'averaging {time_per_line}s per line'
                  f' current serializing index: {serializing_index}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('chr')
    parser.add_argument('purpose')
    args = parser.parse_args()
    extend_dataset(args.chr, args.purpose)

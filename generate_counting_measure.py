import argparse
from collections import defaultdict
import numpy as np
import time
import h5py
import os
from binary_search import search, scan_through_line_for_number

num_chromatin_states = 100


def open_maf_files():
    file_dict = {}
    # chromosome 1 to 22
    for i in range(1, 23):
        filename = f'chr{i}_maf_sequence.csv'
        file_dict[f'chr{i}'] = open(filename, 'r'), os.stat(filename).st_size
    file_dict['chrX'] = open('chrX_maf_sequence.csv', 'r'), os.stat('chrX_maf_sequence.csv').st_size
    file_dict['chrY'] = open('chrY_maf_sequence.csv', 'r'), os.stat('chrY_maf_sequence.csv').st_size
    
    return file_dict


def close_file_dict(file_dict):
    for file, _ in file_dict.values():
        file.close()


def get_line_count(filename):
    count = 0
    with open(filename, 'r') as file:
        for _ in file:
            count += 1
    return count


def get_counts(line):
    count_dict = defaultdict(int)
    tokens = line.strip().lower().split(',')[1:]
    # remove the hg19 base
    tokens.remove(42)
    for base in tokens:
        count_dict[base] += 1
    
    return np.array([count_dict['a'], count_dict['g'], count_dict['c'], count_dict['t']], dtype='uint8')


def generate_counting_measure(coord_filename, output_filename):
    maf_file_dict = open_maf_files()
    flanking_number = 400
    num_basepairs = 1000
    num_channels = 4
    states_list = []
    
    line_count = get_line_count(coord_filename)
    
    stamp = time.time()
    num_missing_states = 0
    serializing_index = 0
    
    zero_state = np.zeros(num_channels, dtype='uint8')
    
    with open(coord_filename, 'r') as coord_file, h5py.File(output_filename, 'w') as hdf5_file:
        feature_group = hdf5_file.create_group('state')
        # uint8 is enough because the count will not exceed the number of species
        feature_data = feature_group.create_dataset('data', (line_count * 2, num_channels, num_basepairs), dtype='uint8')
        
        for line_index, line in enumerate(coord_file):
            tokens = line.split()
            chrom, start, end_exclusive = tokens[0], int(tokens[1]), int(tokens[2])
            real_start = start - flanking_number
            real_end_exclusive = end_exclusive + flanking_number
            # assert num_basepairs == real_end_exclusive - real_start
            counting_states = np.zeros((num_channels, num_basepairs), dtype='uint8')
            coord_to_search = real_start
            start_byteoffset_hint = None
            file, file_bytesize = maf_file_dict[chrom]
            for i in range(num_basepairs):
                if start_byteoffset_hint:
                    result = scan_through_line_for_number(file, start_byteoffset_hint, coord_to_search)
                else:
                    result = search(file, coord_to_search, file_bytesize)
                
                if not result:
                    state_value = zero_state
                    num_missing_states += 1
                else:
                    line, start_byteoffset_hint = result
                    state_value = get_counts(line)
                
                counting_states[i] = state_value
                coord_to_search += 1
            
            states_list.append(counting_states)
            
            if line_index % 100 == 1:
                for vector in states_list:
                    # vector has shape num_channels x num_basepairs
                    feature_data[serializing_index] = vector
                    feature_data[serializing_index + line_count] = vector[:, ::-1]
                    serializing_index += 1
                
                states_list = []
                print(f'{line_index}/{line_count} = {line_index/line_count:.2%} in {time.time() - stamp:.4f}s'
                      f' current serializing index: {serializing_index}',
                      end='\r')
        
        # Serialze the remaining data
        if states_list:
            for vector in states_list:
                # vector has shape num_channels x num_basepairs
                feature_data[serializing_index] = vector
                feature_data[serializing_index + line_count] = vector[:, ::-1]
                serializing_index += 1
            
            print(f'{serializing_index}/{line_count} = {serializing_index/line_count:.2%} in {time.time() - stamp:.4f}s'
                  f' current serializing index: {serializing_index}',
                  end='\r')
    
    print(f'\n-> Number of missing states: {num_missing_states}')
    
    close_file_dict(maf_file_dict)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('coord_filename')
    parser.add_argument('output_filename')
    args = parser.parse_args()
    generate_counting_measure(args.coord_filename, args.output_filename)

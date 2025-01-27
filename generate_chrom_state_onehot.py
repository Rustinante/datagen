import argparse
import numpy as np
import time
import h5py
import os
from bed_binary_search import search, get_start_end_location_from_line, scan_through_line_for_number

num_chromatin_states = 100


def open_chrom_state_files():
    file_dict = {}
    # chromosome 1 to 22
    for i in range(1, 23):
        filename = f'chr{i}_segmentation.bed'
        file_dict[f'chr{i}'] = open(filename, 'r'), os.stat(filename).st_size
    file_dict['chrX'] = open('chrX_segmentation.bed', 'r'), os.stat('chrX_segmentation.bed').st_size
    file_dict['chrY'] = open('chrY_segmentation.bed', 'r'), os.stat('chrY_segmentation.bed').st_size
    
    return file_dict


def close_file_dict(file_dict):
    for file, _ in file_dict.values():
        file.close()


def get_chrom_state_onehot_mapping():
    mapping = {}
    for i in range(num_chromatin_states):
        encoding = np.zeros(num_chromatin_states, dtype='uint8')
        encoding[i] = 1
        mapping[f'U{i+1}'] = encoding
    return mapping


def get_line_count(filename):
    count = 0
    with open(filename, 'r') as file:
        for _ in file:
            count += 1
    return count


def generate_one_hot(coord_filename):
    chrom_state_file_dict = open_chrom_state_files()
    flanking_number = 400
    chrom_state_mapping = get_chrom_state_onehot_mapping()
    states_list = []
    num_basepairs = 1000
    
    line_count = get_line_count(coord_filename)
    
    stamp = time.time()
    num_missing_states = 0
    serializing_index = 0
    
    with open(coord_filename, 'r') as coord_file, h5py.File('chrom_states_onehot', 'w') as hdf5_file:
        feature_group = hdf5_file.create_group('state')
        # uint8 is enough because there are only 100 states
        feature_data = feature_group.create_dataset('data', (line_count, num_basepairs, num_chromatin_states),
                                                    dtype='uint8')
        
        for line_index, line in enumerate(coord_file):
            tokens = line.split()
            chrom, start, end_exclusive = tokens[0], int(tokens[1]), int(tokens[2])
            real_start = start - flanking_number
            real_end_exclusive = end_exclusive + flanking_number
            # assert num_basepairs == real_end_exclusive - real_start
            states = np.zeros((num_basepairs, num_chromatin_states), dtype='uint8')
            collected_basepair_count = 0
            coord_to_search = real_start
            start_byteoffset_hint = None
            file, file_bytesize = chrom_state_file_dict[chrom]
            while collected_basepair_count < num_basepairs:
                if start_byteoffset_hint:
                    result = scan_through_line_for_number(file, start_byteoffset_hint, coord_to_search)
                else:
                    result = search(file, coord_to_search, file_bytesize)
                
                if not result:
                    line = f'{chrom} {coord_to_search} {coord_to_search + 1} U96'
                    num_missing_states += 1
                else:
                    line, start_byteoffset_hint = result
                
                start, end_exclusive = get_start_end_location_from_line(line)
                state_value = chrom_state_mapping[line.strip().split()[-1]]
                
                repetitions = min(end_exclusive, real_end_exclusive) - coord_to_search
                states[collected_basepair_count:collected_basepair_count + repetitions] = state_value
                collected_basepair_count += repetitions
                coord_to_search += repetitions
            
            states_list.append(states)
            
            if line_index % 100 == 1:
                for matrix in states_list:
                    feature_data[serializing_index] = matrix
                    serializing_index += 1
                    
                states_list = []
                print(f'{line_index}/{line_count} = {line_index/line_count:.2%} in {time.time() - stamp:.4f}s'
                      f' current serializing index: {serializing_index}',
                      end='\r')
                
        if states_list:
            for matrix in states_list:
                feature_data[serializing_index] = matrix
                serializing_index += 1
    
            print(f'{serializing_index}/{line_count} = {serializing_index/line_count:.2%} in {time.time() - stamp:.4f}s'
                  f' current serializing index: {serializing_index}')
    
    print(f'\n-> Number of missing states: {num_missing_states}')
    
    close_file_dict(chrom_state_file_dict)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('coord_filename')
    args = parser.parse_args()
    generate_one_hot(args.coord_filename)

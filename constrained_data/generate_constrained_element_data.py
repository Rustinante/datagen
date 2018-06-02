import argparse
import numpy as np
import time
import h5py
import os
from constrained_data_binary_search import search, get_start_end_location_from_line, scan_through_line_for_number

num_chromatin_states = 100


def open_chrom_files():
    file_dict = {}
    
    # chromosome 1 to 22
    for i in range(1, 23):
        filename = f'chr{i}_phast_cons.txt'
        file_dict[f'chr{i}'] = open(filename, 'r'), os.stat(filename).st_size
    
    file_dict['chrX'] = open('chrX_phast_cons.txt', 'r'), os.stat('chrX_phast_cons.txt').st_size
    file_dict['chrY'] = open('chrY_phast_cons.txt', 'r'), os.stat('chrY_phast_cons.txt').st_size
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


def generate(coord_filename):
    chrom_file_dict = open_chrom_files()
    flanking_number = 400
    states_list = []
    num_basepairs = 1000
    
    line_count = get_line_count(coord_filename)
    
    stamp = time.time()
    num_missing_states = 0
    serializing_index = 0
    
    with open(coord_filename, 'r') as coord_file, h5py.File('constrained_elements.hdf5', 'w') as hdf5_file:
        feature_group = hdf5_file.create_group('state')
        # uint8 is enough because there are only 100 states
        feature_data = feature_group.create_dataset('data', (line_count, num_basepairs), dtype='uint8')
        
        for line_index, line in enumerate(coord_file):
            # each line will be of the form
            # 585 chr1 11991 11995 lod=12 240
            # Note that in the actual file the spaces are actually tab characters.
            tokens = line.split()
            chrom, start, end_exclusive = tokens[0], int(tokens[1]), int(tokens[2])
            real_start = start - flanking_number
            real_end_exclusive = end_exclusive + flanking_number
            # assert num_basepairs == real_end_exclusive - real_start
            states = np.zeros(num_basepairs, dtype=np.int8)
            collected_basepair_count = 0
            coord_to_search = real_start
            start_byteoffset_hint = None
            file, file_bytesize = chrom_file_dict[chrom]
            while collected_basepair_count < num_basepairs:
                if start_byteoffset_hint:
                    coord_is_present, line, start_byteoffset_hint = scan_through_line_for_number(file, start_byteoffset_hint,
                                                                                                 coord_to_search)
                else:
                    coord_is_present, line, start_byteoffset_hint = search(file, coord_to_search, file_bytesize)
                
                if coord_is_present:
                    state_value = 1
                    _, current_end_exclusive = get_start_end_location_from_line(line)
                    repetitions = min(current_end_exclusive, real_end_exclusive) - coord_to_search
                else:
                    state_value = 0
                    if line:
                        current_end_exclusive, _ = get_start_end_location_from_line(line)
                        repetitions = min(current_end_exclusive, real_end_exclusive) - coord_to_search
                    else:
                        repetitions = real_end_exclusive - coord_to_search
                
                states[collected_basepair_count:collected_basepair_count + repetitions] = state_value
                collected_basepair_count += repetitions
                coord_to_search += repetitions
            
            states_list.append(states)
            
            if line_index % 100 == 1:
                for vector in states_list:
                    feature_data[serializing_index] = vector
                    serializing_index += 1
                
                states_list = []
                print(f'{line_index}/{line_count} = {line_index/line_count:.2%} in {time.time() - stamp:.4f}s'
                      f' current serializing index: {serializing_index}',
                      end='\r')
    
    print(f'\n-> Number of missing states: {num_missing_states}')
    
    close_file_dict(chrom_file_dict)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('coord_filename')
    args = parser.parse_args()
    generate(args.coord_filename)

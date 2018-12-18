import argparse
import h5py
import numpy as np
import os
import time

from collections import defaultdict

from binary_search import search, scan_through_line_for_number

from line_cache import LineCache


def get_line_count(filename):
    count = 0
    with open(filename, 'r') as file:
        for _ in file:
            count += 1
    return count


def distribute_tokens_to_species_sequence(tokens, species_sequence):
    for species_index, letter in enumerate(tokens):
        species_sequence[species_index] += letter
        
        
def extend_dataset(chr, purpose):
    checkpoint_time_str = time.strftime('%a %b %d %Y %H:%M:%S UTC%z', time.localtime(time.time()))
    print('Current time: {}'.format(checkpoint_time_str))
    
    cache = LineCache()
    
    coordinate_filename = os.path.join('data', '{}_{}'.format(chr, purpose))
    alignment_filename = '{}_maf_sequence.csv'.format(chr)
    
    species_file_dict = {}
    
    print('=> coordinate_filename: {}'.format(coordinate_filename))
    print('=> alignment_filename: {}'.format(alignment_filename))
    file_byte_size = os.stat(alignment_filename).st_size

    total_line_count = get_line_count(coordinate_filename)
    
    dir_name = 'gkm_fasta_{}'.format(chr)
    os.makedirs(dir_name, exist_ok=True)
    
    with open(coordinate_filename, 'r') as file, open(alignment_filename, 'r') as alignment_file:
        
        header = alignment_file.readline().strip().split(',')
        del header[0]
        assert len(header) == 100
        
        for index, species_code in enumerate(header):
            species_filename = '{}_{}_{}_{}.fasta'.format(index, species_code, chr, purpose)
            print('=> Creating {} under {}'.format(species_filename, dir_name))
            species_file_dict[index] = open(os.path.join(dir_name, species_filename), 'w')
        
        processed_line_count = 0
        start_time = time.time()
        flanking_number = 400
        
        for line in file:
            processed_line_count += 1
            start_coordinate = int(line.strip().split(',')[0])
            
            species_sequence = defaultdict(str)
            
            start_line_hint = None
            
            for letter_index, coordinate in enumerate(range(start_coordinate - flanking_number, start_coordinate + 200 + flanking_number)):
                if coordinate in cache:
                    cached_result = cache[coordinate]
                    tokens = cached_result[0]
                    distribute_tokens_to_species_sequence(tokens, species_sequence)
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
                    del tokens[0]
                    # assert len(tokens) == 100
                    
                    distribute_tokens_to_species_sequence(tokens, species_sequence)
                    
                    cache[coordinate] = (tokens, start_line_hint)
                
                else:
                    print('Failed to find the letters at coordinate {}. '
                          'Substituting with N instead'.format(coordinate))
                    tokens = ['N'] * 100
                    distribute_tokens_to_species_sequence(tokens, species_sequence)
                    
            for species_index, _ in enumerate(header):
                species_file_dict[species_index].write('>middle 200 bp start coordinate {}\n'
                                                       '{}\n'
                                                       .format(start_coordinate, species_sequence[species_index]))
                
            if processed_line_count % 1000 == 0:
                elapsed_time = time.time() - start_time
                time_per_line = elapsed_time / processed_line_count
                print('Processed [{}/{}] lines in {:5f}s, averaging: {:5f}s per line, {:3f} done.'
                      .format(processed_line_count, total_line_count,
                              elapsed_time, time_per_line, processed_line_count/total_line_count))
        
    for filename, species_file in species_file_dict.items():
        print('=> Closing {}'.format(filename))
        species_file.close()
        
    print('=> Done!')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('chr')
    parser.add_argument('purpose')
    args = parser.parse_args()
    extend_dataset(args.chr, args.purpose)
